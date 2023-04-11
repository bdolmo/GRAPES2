import os
import sys
from pathlib import Path
import re
import logging
import subprocess
import pandas as pd
from functools import partial, reduce
from modules.sample import Sample
from modules.params import *
from modules.gc_content import annotate_gc
from modules.mappability import annotate_mappability


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
    unified_depth_name = ("{}.read.counts.bed").format(analysis_dict["output_name"])
    unified_raw_depth = str(Path(analysis_dict["output_dir"]) / unified_depth_name)
    analysis_dict["unified_raw_depth"] = unified_raw_depth

    per_base_coverage_name = ("{}.per.base.coverage.bed").format(
        analysis_dict["output_name"]
    )
    per_base_coverage_file = str(
        Path(analysis_dict["output_dir"]) / per_base_coverage_name
    )
    analysis_dict["per_base_coverage"] = per_base_coverage_file

    cmd = ("{} -i {} -o {} -n {} -g {} -b {} -t {} -c -d ").format(
        ngs_utils_dict["targetdepth"],
        analysis_dict["bam_dir"],
        analysis_dict["output_dir"],
        analysis_dict["output_name"],
        analysis_dict["reference"],
        analysis_dict["ready_bed"],
        analysis_dict["threads"],
    )

    if not os.path.isfile(unified_raw_depth) and not os.path.isfile(
        per_base_coverage_file
    ):

        msg = (" INFO: Extracting coverage for {} samples").format(
            analysis_dict["output_name"]
        )
        logging.info(msg)

        p1 = subprocess.run(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        output = p1.stdout.decode("UTF-8")
        error = p1.stderr.decode("UTF-8")

    summary_log_name = "summary_metrics.log"
    summary_log = str(Path(analysis_dict["output_dir"]) / summary_log_name)
    with open(summary_log) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("SAMPLE"):
                continue
            tmp = line.split("\t")
            sample_name = tmp[0]
            for sample in sample_list:
                if sample.name == sample_name:
                    sample.add("ontarget_reads", int(tmp[2]))
                    sample.add("mean_coverage", float(tmp[5]))
                    sample.add("mean_coverage_X", float(tmp[9]))
        f.close()
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
            msg = (" INFO: Extracting coverage from sample {}").format(sample.name)
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
                    msg = " ERROR: Could not extract coverage"
                    logging.error(error)
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
