import os
import sys
import logging
import glob
import subprocess
import json
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd


INTERMEDIATE_FILE_PATTERNS = (
    # Per-sample TargetDepth outputs. The coverage files are especially large.
    "*.bam_counts.bed",
    "*.bam_isizes.bed",
    "*.bam_coverage.bed",
    # Cohort matrices used only while normalizing and calling CNVs.
    "*.read.counts.bed",
    "*.per.base.coverage.bed",
    "*.normalized.depth.bed",
    "*.normalized.per.base.bed",
    "*.normalized.offtarget.bed",
    "*.normalized.all.bed",
    "*.ratios.bed",
    # Derived BED inputs and temporary call files.
    "*.splitted.gc.bed",
    "*.splitted.map.bed",
    "*.tmp.rawcalls.bed",
    # Off-target analysis intermediates.
    "GRAPES2.offtarget.raw.counts.bed",
    "ontarget_paded_400bp.bed",
    "ontarget_paded_400bp_centromers_patches.bed",
    "offtarget_unprocessed.bed",
    "offTarget.pseudowindowed.bed",
    "offTarget.pseudowindowed.gc.bed",
    "offTarget.pseudowindowed.gc.map.bed",
)

STAGE_MANIFEST_NAME = ".grapes2.stage_state.json"


def validate_delimited_file(
    file_path,
    required_columns=None,
    min_columns=None,
    min_data_rows=1,
    has_header=True,
    delimiter="\t",
):
    """Validate a delimited stage output without loading the complete file."""
    if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
        return False

    required_columns = set(required_columns or [])
    data_rows = 0
    try:
        with open(file_path, "r", encoding="utf-8") as handle:
            first_line = handle.readline().rstrip("\r\n")
            if not first_line:
                return False

            first_fields = first_line.lstrip("#").split(delimiter)
            if min_columns is not None and len(first_fields) < min_columns:
                return False
            if has_header and not required_columns.issubset(first_fields):
                return False
            if not has_header:
                data_rows = 1

            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\r\n").split(delimiter)
                if min_columns is not None and len(fields) < min_columns:
                    return False
                data_rows += 1
                if data_rows >= min_data_rows:
                    break
    except (OSError, UnicodeError):
        return False

    return data_rows >= min_data_rows


def validate_vcf_file(file_path, min_records=0):
    """Validate that a VCF has a complete column header and enough records."""
    if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
        return False

    found_header = False
    record_count = 0
    try:
        with open(file_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM\t"):
                    if len(line.rstrip("\r\n").split("\t")) < 10:
                        return False
                    found_header = True
                    continue
                if not line.strip():
                    continue
                if not found_header or len(line.rstrip("\r\n").split("\t")) < 10:
                    return False
                record_count += 1
    except (OSError, UnicodeError):
        return False

    return found_header and record_count >= min_records


@contextmanager
def atomic_output_path(final_path):
    """Yield a same-directory temporary path and atomically promote it on success."""
    final_path = Path(final_path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_path = tempfile.mkstemp(
        prefix=f".{final_path.name}.",
        suffix=".tmp",
        dir=str(final_path.parent),
    )
    os.close(fd)
    temporary_path = Path(temporary_path)
    try:
        yield str(temporary_path)
        if not temporary_path.is_file():
            raise RuntimeError(f"Atomic output was not created: {temporary_path}")
        os.replace(temporary_path, final_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def atomic_write_dataframe(dataframe, output_path, **to_csv_kwargs):
    """Write a DataFrame through a temporary file and atomically replace output."""
    with atomic_output_path(output_path) as temporary_path:
        dataframe.to_csv(temporary_path, **to_csv_kwargs)


def _manifest_file_record(file_path, output_root):
    path = Path(file_path)
    try:
        display_path = str(path.resolve().relative_to(output_root))
    except ValueError:
        display_path = str(path.resolve())
    stat = path.stat()
    return {
        "path": display_path,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def stage_manifest_matches(
    output_dir,
    stage_name,
    output_paths,
    input_paths=None,
    metadata=None,
):
    """Return True when a recorded stage still matches its inputs and outputs."""
    manifest_path = Path(output_dir) / STAGE_MANIFEST_NAME
    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
        record = manifest["stages"][stage_name]
    except (OSError, KeyError, TypeError, ValueError):
        return False

    normalized_metadata = json.loads(json.dumps(metadata or {}, sort_keys=True))
    if record.get("metadata", {}) != normalized_metadata:
        return False

    output_root = Path(output_dir).resolve()
    recorded_outputs = {
        item.get("path"): item for item in record.get("outputs", [])
    }
    recorded_inputs = {
        item.get("path"): item for item in record.get("inputs", [])
    }

    for paths, recorded in (
        (output_paths, recorded_outputs),
        (input_paths or [], recorded_inputs),
    ):
        for file_path in paths:
            path = Path(file_path)
            if not path.is_file():
                return False
            current = _manifest_file_record(path, output_root)
            if recorded.get(current["path"]) != current:
                return False

    return True


def update_stage_manifest(
    output_dir,
    stage_name,
    output_paths,
    metadata=None,
    input_paths=None,
):
    """Record a validated stage completion using an atomic JSON manifest update."""
    manifest_path = Path(output_dir) / STAGE_MANIFEST_NAME
    manifest = {"version": 1, "stages": {}}
    if manifest_path.is_file():
        try:
            with open(manifest_path, "r", encoding="utf-8") as handle:
                loaded_manifest = json.load(handle)
            if isinstance(loaded_manifest, dict):
                manifest = loaded_manifest
                manifest.setdefault("version", 1)
                manifest.setdefault("stages", {})
        except (OSError, ValueError, TypeError):
            logging.warning(
                f" WARNING: Replacing invalid GRAPES2 stage manifest {manifest_path}"
            )

    output_root = Path(output_dir).resolve()
    recorded_outputs = [
        _manifest_file_record(output_path, output_root)
        for output_path in output_paths
        if Path(output_path).is_file()
    ]
    recorded_inputs = [
        _manifest_file_record(input_path, output_root)
        for input_path in (input_paths or [])
        if Path(input_path).is_file()
    ]

    stage_record = {
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "outputs": recorded_outputs,
        "inputs": recorded_inputs,
        "metadata": metadata or {},
    }
    manifest["stages"][stage_name] = stage_record

    with atomic_output_path(manifest_path) as temporary_path:
        with open(temporary_path, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")


def remove_tmp_files(input_dir):
    """Remove reproducible intermediates and report files/bytes reclaimed."""
    removed_files = []
    reclaimed_bytes = 0
    matched_files = set()

    for pattern in INTERMEDIATE_FILE_PATTERNS:
        matched_files.update(glob.glob(os.path.join(input_dir, pattern)))

    for file_path in sorted(matched_files):
        if not os.path.isfile(file_path):
            continue

        try:
            file_size = os.path.getsize(file_path)
            os.remove(file_path)
        except OSError as exc:
            logging.warning(
                f" WARNING: Could not remove GRAPES2 intermediate {file_path}: {exc}"
            )
            continue

        removed_files.append(file_path)
        reclaimed_bytes += file_size

    return removed_files, reclaimed_bytes



def sort_bed_file(input_bed):
    """ """

    output_bed = input_bed.replace(".bed", ".tmp.bed")

    # Load the BED file into a DataFrame
    df = pd.read_csv(input_bed, sep='\t', names=['chr', 'start', 'end', 'name'], header=None)

    if df.empty:
        return df

    df['name'] = df['name'].replace(";", "_")
    # Replace 'chrX' and 'chrY' with temporary placeholders
    df['chr'] = df['chr'].replace({'chrM': 'chr0', 'chrX': 'chr23', 'chrY': 'chr24'})

    # Remove 'chr' prefix for sorting but keep it in a separate column
    df['chr_num'] = df['chr'].str.replace('chr', '').astype(int)

    # Sort by chromosomal position
    df = df.sort_values(['chr_num', 'start', 'end'])

    # Drop the temporary column
    df = df.drop(columns=['chr_num'])

    # Replace temporary placeholders with 'chrX' and 'chrY'
    df['chr'] = df['chr'].replace({'chr0':'chrM', 'chr23': 'chrX', 'chr24': 'chrY'})

    # Save the sorted DataFrame back to a BED file
    df.to_csv(output_bed, sep='\t', header=False, index=False)

    os.remove(input_bed)
    os.rename(output_bed, input_bed)


def signal_to_noise(data):
    """ """
    median = np.median(data)
    std = np.std(data)
    if std == 0:
        s2n = 0
    else:
        s2n = abs(median / std)
    return round(s2n, 3)


def remove_bed_header(file, pattern):
    """ """
    no_header_file = file.replace(".bed", ".noheader.bed")
    nh = open(no_header_file, "w")
    with open(file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(pattern):
                continue
            nh.write(line + "\n")
    f.close()
    nh.close()
    return no_header_file


def assign_genotype_based_on_cn(cn):
    """
    Simple genotype assignment based on copy number
    """
    if cn == 0:
        gt = "1/1"
    elif cn == 1:
        gt = "0/1"
    else:
        gt = "./1"
    return gt


def validate_bed(bed):
    """ """
    pass


class MissingInputBamFiles(Exception):
    pass


class NgsUtilsError(Exception):
    pass


def check_executable(program, dump_messages=True):
    """ """
    bashCommand = ("which {}").format(program)
    p1 = subprocess.run(
        bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output = p1.stdout.decode("UTF-8")
    error = p1.stderr.decode("UTF-8")
    executable = False
    if not error:
        if output:
            msg = (" INFO: Found executable for {}").format(os.path.basename(program))
            logging.info(msg)
            pass
        else:
            msg = (" ERROR: Unable to execute {}").format(os.path.basename(program))
            raise NgsUtilsError(msg)
            # print(msg)
    else:
        msg = (" ERROR: Unable to execute {}").format(os.path.basename(program))
        raise NgsUtilsError(msg)


def get_bam_files(input_dir):
    """
    Get bam files from input dir
    """
    bam_list = glob.glob(input_dir + "/*.bam")

    if not bam_list:
        msg = (" ERROR: missing input bam files from {} directory").format(input_dir)
        # logging.error(msg)
        raise MissingInputBamFiles(msg)

    return bam_list
