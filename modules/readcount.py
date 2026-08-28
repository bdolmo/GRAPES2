import os
import sys
from pathlib import Path
import re
import logging
import subprocess
import shutil
import shlex
import tempfile
import pandas as pd
from functools import partial, reduce
from modules.sample import Sample
from modules.params import *
from modules.gc_content import annotate_gc
from modules.mappability import annotate_mappability
from modules.utils import (
    atomic_output_path,
    sort_bed_file,
    stage_manifest_matches,
    update_stage_manifest,
    validate_delimited_file,
)


READ_DEPTH_META_COLUMNS = ["chr", "start", "end", "exon", "gc", "map"]


def _read_depth_outputs_valid(
    unified_raw_depth,
    summary_log,
    sample_list,
    per_base_coverage_file=None,
):
    sample_names = [sample.name for sample in sample_list]
    required_depth_columns = READ_DEPTH_META_COLUMNS + sample_names
    if not validate_delimited_file(
        unified_raw_depth,
        required_columns=required_depth_columns,
        min_columns=len(required_depth_columns),
    ):
        return False
    if not validate_delimited_file(
        summary_log,
        required_columns=[
            "SAMPLE",
            "READS_ON_TARGET",
            "%ROI",
            "MEAN_COVERAGE",
            "MEAN_COVERAGE_X",
        ],
        min_columns=5,
    ):
        return False
    if per_base_coverage_file is not None:
        if not validate_delimited_file(
            per_base_coverage_file,
            required_columns=required_depth_columns,
            min_columns=len(required_depth_columns),
        ):
            return False
    return True


def launch_read_depth(sample_list, analysis_dict, ngs_utils, ann_dict):
    """ """
    sample_list, analysis_dict = extract_read_depth(
        sample_list, analysis_dict, ngs_utils, ann_dict
    )
    return sample_list, analysis_dict


def extract_read_depth(sample_list, analysis_dict, ngs_utils_dict, ann_dict):
    """ """
    analysis_dict = annotate_gc(analysis_dict)
    analysis_dict = annotate_mappability(analysis_dict, ann_dict)
    unified_depth_name = f'{analysis_dict["output_name"]}.read.counts.bed'
    unified_raw_depth = str(Path(analysis_dict["output_dir"]) / unified_depth_name)
    analysis_dict["unified_raw_depth"] = unified_raw_depth


    needs_per_base_coverage = bool(analysis_dict.get("single_exon_cnv", True))
    per_base_coverage_name = f'{analysis_dict["output_name"]}.per.base.coverage.bed'
    per_base_coverage_file = str(
        Path(analysis_dict["output_dir"]) / per_base_coverage_name
    )
    analysis_dict["per_base_coverage"] = (
        per_base_coverage_file if needs_per_base_coverage else None
    )

    summary_log_name = "summary_metrics.log"
    summary_log = str(Path(analysis_dict["output_dir"]) / summary_log_name)
    required_per_base_file = (
        per_base_coverage_file if needs_per_base_coverage else None
    )
    read_depth_inputs = [
        path
        for path in [
            analysis_dict.get("ready_bed"),
            analysis_dict.get("reference"),
            analysis_dict.get("bam_dir"),
            *[getattr(sample, "bam", None) for sample in sample_list],
        ]
        if path and os.path.isfile(path)
    ]
    read_depth_metadata = {
        "per_base_coverage": needs_per_base_coverage,
        "samples": [sample.name for sample in sample_list],
    }
    expected_outputs = [unified_raw_depth, summary_log]
    if needs_per_base_coverage:
        expected_outputs.append(per_base_coverage_file)
    read_depth_current = stage_manifest_matches(
        analysis_dict["output_dir"],
        "read_depth",
        expected_outputs,
        input_paths=read_depth_inputs,
        metadata=read_depth_metadata,
    )

    if (
        analysis_dict.get("force", False)
        or not read_depth_current
        or not _read_depth_outputs_valid(
            unified_raw_depth,
            summary_log,
            sample_list,
            required_per_base_file,
        )
    ):
        msg = f' INFO: Extracting coverage for {analysis_dict["output_name"]}'
        logging.info(msg)
        if not needs_per_base_coverage:
            logging.info(
                " INFO: Per-base coverage disabled because single-exon CNV "
                "analysis is not enabled"
            )

        stage_dir = tempfile.mkdtemp(
            prefix=".targetdepth-",
            dir=analysis_dict["output_dir"],
        )
        targetdepth_log = str(
            Path(analysis_dict["output_dir"]) / "targetdepth.console.log"
        )
        cmd = [
            ngs_utils_dict["targetdepth"],
            "-i",
            analysis_dict["bam_dir"],
            "-o",
            stage_dir,
            "-n",
            analysis_dict["output_name"],
            "-g",
            analysis_dict["reference"],
            "-b",
            analysis_dict["ready_bed"],
            "-t",
            str(analysis_dict["threads"]),
            "-c",
        ]
        if needs_per_base_coverage:
            cmd.append("-d")
        logging.info(f" INFO: TargetDepth command: {shlex.join(cmd)}")

        try:
            with open(targetdepth_log, "wb") as log_handle:
                completed = subprocess.run(
                    cmd,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                )
            staged_unified_depth = str(Path(stage_dir) / unified_depth_name)
            staged_summary_log = str(Path(stage_dir) / summary_log_name)
            staged_per_base = (
                str(Path(stage_dir) / per_base_coverage_name)
                if needs_per_base_coverage
                else None
            )
            if completed.returncode != 0 or not _read_depth_outputs_valid(
                staged_unified_depth,
                staged_summary_log,
                sample_list,
                staged_per_base,
            ):
                raise RuntimeError(
                    "TargetDepth failed or produced invalid outputs; console output "
                    f"is available at {targetdepth_log}"
                )

            for staged_file in Path(stage_dir).iterdir():
                if staged_file.is_file():
                    os.replace(
                        str(staged_file),
                        str(Path(analysis_dict["output_dir"]) / staged_file.name),
                    )
        finally:
            shutil.rmtree(stage_dir, ignore_errors=True)
    else:
        logging.info(" INFO: Reusing validated TargetDepth outputs")

    if needs_per_base_coverage and os.path.isfile(per_base_coverage_file):
        if not check_first_line(per_base_coverage_file):
            with atomic_output_path(per_base_coverage_file) as temporary_path:
                with open(temporary_path, "w", encoding="utf-8") as output_handle:
                    with open(
                        per_base_coverage_file, "r", encoding="utf-8"
                    ) as input_handle:
                        for line in input_handle:
                            if line.startswith("chr\tstart"):
                                output_handle.write("#" + line)
                            else:
                                output_handle.write(line)

    if not _read_depth_outputs_valid(
        unified_raw_depth,
        summary_log,
        sample_list,
        required_per_base_file,
    ):
        raise RuntimeError("TargetDepth outputs failed final validation")

    stage_outputs = [unified_raw_depth, summary_log]
    if needs_per_base_coverage:
        stage_outputs.append(per_base_coverage_file)
    update_stage_manifest(
        analysis_dict["output_dir"],
        "read_depth",
        stage_outputs,
        metadata=read_depth_metadata,
        input_paths=read_depth_inputs,
    )

    with open(summary_log) as f:
        header_line = f.readline().rstrip("\n")
        header_cols = header_line.split("\t")
        header_idx = {name: idx for idx, name in enumerate(header_cols)}

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            tmp = line.split("\t")
            sample_name = tmp[header_idx["SAMPLE"]].replace(".bam", "")

            for sample in sample_list:
                if sample.name != sample_name:
                    continue

                # Robust field extraction by header names.
                enrichment = float(tmp[header_idx["%ROI"]])
                ontarget_reads = int(tmp[header_idx["READS_ON_TARGET"]])
                mean_coverage = float(tmp[header_idx["MEAN_COVERAGE"]])
                mean_coverage_X = float(tmp[header_idx["MEAN_COVERAGE_X"]])

                # Sex inference from chrX/autosomal depth ratio:
                # female ~1.0, male ~0.5
                gender = "Undefined"
                x_ratio = round(mean_coverage_X / mean_coverage, 3) if mean_coverage > 0 else 0
                if x_ratio >= 0.75:
                    gender = "Female"
                elif x_ratio > 0 and x_ratio <= 0.65:
                    gender = "Male"

                sample.add("enrichment", enrichment)
                sample.add("ontarget_reads", ontarget_reads)
                sample.add("mean_coverage", mean_coverage)
                sample.add("mean_coverage_X", mean_coverage_X)
                sample.add("gender", gender)
                msg = f" INFO: {sample.name}\tGender_ratio_X:{x_ratio}\tGender:{gender}"
                logging.info(msg)

    return sample_list, analysis_dict


def unify_read_depths(sample_list, analysis_dict):
    """ """
    df_list = []
    df_ready_bed = pd.read_csv(
        analysis_dict["ready_bed"],
        header=0,
        sep="\t",
        names=["chr", "start", "end", "region", "gc", "map"],
    )

    for sample in sample_list:
        df = pd.read_csv(
            sample.region_coverage,
            compression="gzip",
            header=0,
            sep="\t",
            names=["chr", "start", "end", "region", sample.name],
        )
        df["gc"] = df_ready_bed["gc"]
        df["map"] = df_ready_bed["map"]
        df = df[["chr", "start", "end", "region", "gc", "map", sample.name]]
        df_list.append(df)

    merge = partial(
        pd.merge, on=["chr", "start", "end", "region", "gc", "map"], how="outer"
    )
    merged_df = reduce(merge, df_list)
    merged_df = merged_df[(df["end"] - df["start"]) > 10]
    unified_depth_name = ("{}.region.read.depth.bed").format(
        analysis_dict["output_name"]
    )
    unified_raw_depth = str(Path(analysis_dict["output_dir"]) / unified_depth_name)
    merged_df.to_csv(unified_raw_depth, sep="\t", index=False)

    analysis_dict["unified_raw_depth"] = unified_raw_depth
    return analysis_dict


def extract_read_depth_exome(sample_list, analysis_dict, ngs_utils_dict, ann_dict):
    """ """

    analysis_dict = annotate_gc(analysis_dict)

    analysis_dict = annotate_mappability(analysis_dict, ann_dict)

    all_coverage_files = []
    for sample in sample_list:

        # Mosdepth output files
        mosdepth_summary_name = ("{}{}").format(sample.name, ".mosdepth.summary.txt")
        mosdepth_region_name = ("{}{}").format(sample.name, ".regions.bed.gz")
        mosdepth_per_base_name = ("{}{}").format(sample.name, ".per-base.bed.gz")

        # ROI mean coverage
        sample_region_file = Path(sample.sample_folder) / mosdepth_region_name

        # Per base coverage
        sample_per_base_file = Path(sample.sample_folder) / mosdepth_per_base_name

        sample.add("region_coverage", sample_region_file)
        sample.add("base_coverage", sample_per_base_file)

        mosdepth_output = str(Path(sample.sample_folder) / sample.name)

        cmd = ("{} --fast-mode --by {} {} {}").format(
            ngs_utils_dict["mosdepth"],
            analysis_dict["ready_bed"],
            mosdepth_output,
            sample.bam,
        )

        if not os.path.isfile(sample_per_base_file):
            msg = f" INFO: Extracting coverage from sample {sample.name}"
            logging.info(msg)
            p1 = subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            output = p1.stdout.decode("UTF-8")
            error = p1.stderr.decode("UTF-8")
            total_reads = ""
            if not error:
                msg = " INFO: Coverage extraction ended successfully"
                logging.info(msg)
            else:
                if re.search("error", error):
                    msg = f" INFO: {error.strip()}"
                    logging.error(msg)
                    msg = " INFO: Could not extract coverage"
                    logging.error(msg)
                else:
                    tmp = error.split("\n")
                    for line in tmp:
                        if line.startswith("Read"):
                            m = re.search(r"\d+", line)
                            total_reads = m.group()
        else:
            msg = (" INFO: Skipping coverage extraction from sample {}").format(
                sample.name
            )
            logging.info(msg)

    analysis_dict = unify_read_depths(sample_list, analysis_dict)
    return sample_list, analysis_dict


def check_first_line(filename):
    with open(filename, 'r') as file:
        first_line = file.readline().strip()  # Reads the first line and removes any leading/trailing whitespace
        if first_line.startswith("#"):  # Checks if the first line starts with "#"
            return True
        else:
            return False
