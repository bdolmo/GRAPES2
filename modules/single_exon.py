#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
from datetime import datetime
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
import pybedtools
import math
from modules.utils import remove_bed_header, signal_to_noise
from modules.plot import plot_single_exon


def get_samples_from_header(coverage_bed) -> dict():
    """
    Function that returns a dict with a sample name as key and an index as value
    """
    samples_dict = dict()
    with open(coverage_bed) as f:
        header = f.readline().rstrip("\n")
        tmp = header.split("\t")
        for i in range(6, len(tmp)):
            sample_name = tmp[i]
            samples_dict[sample_name] = i

    return samples_dict


def evaluate_single_exon_cnvs(sample_list, analysis_dict):
    """ """
    sample_idx_dict = get_samples_from_header(analysis_dict["normalized_per_base"])

    for sample in sample_list:

        references = sample.references

        with open(sample.raw_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                tmp = line.split("\t")
                n_regions = tmp[4]
                exon = tmp[3]
                if n_regions == "1":
                    candidate_name = ("{}.bed").format("_".join(tmp[0:3]))
                    candidate_bed = str(Path(sample.sample_folder) / candidate_name)
                    o = open(candidate_bed, "w")
                    o.write("\t".join(tmp[0:5]) + "\n")
                    o.close()
                    a = pybedtools.BedTool(analysis_dict["normalized_per_base"])
                    b = pybedtools.BedTool(candidate_bed)
                    c = a.intersect(b, wb=True, stream=True)

                    plot_dict = defaultdict()
                    plot_dict["position"] = []
                    plot_dict["sample"] = []
                    plot_dict["log2_ratio"] = []
                    plot_dict["class"] = []

                    log2_ratios_case = []
                    log2_ratios_controls = []

                    for line in iter(c):
                        line = str(line)
                        line = line.rstrip()
                        tmp = line.split("\t")
                        start = int(tmp[1])

                        baseline_cov_list = []
                        sample_cov_list = []
                        sample_cov = float(tmp[sample_idx_dict[sample.name]])
                        sample_cov_list.append(sample_cov)

                        for control in references:
                            idx = sample_idx_dict[control[0]]
                            cov_control = float(tmp[idx])
                            control_cov_list = []
                            for control2 in references:
                                idx2 = sample_idx_dict[control2[0]]
                                cov2 = float(tmp[idx2])
                                if control == control2:
                                    continue
                                control_cov_list.append(cov2)
                            median_controls = np.median(control_cov_list)
                            if median_controls == 0:
                                ratio_controls = 0.1
                            else:
                                ratio_controls = cov_control / median_controls
                            control_log2_ratio = round(math.log2(ratio_controls), 3)
                            log2_ratios_controls.append(control_log2_ratio)

                            baseline_cov_list.append(cov_control)
                            plot_dict["position"].append(start)
                            plot_dict["sample"].append(control[0])
                            plot_dict["log2_ratio"].append(control_log2_ratio)
                            plot_dict["class"].append("control")
                        median_cov_baseline = np.median(baseline_cov_list)
                        std_cov_baseline = np.std(baseline_cov_list)
                        signal_to_noise_baseline = signal_to_noise(baseline_cov_list)

                        median_cov_sample = np.median(sample_cov_list)
                        std_cov_sample = np.std(sample_cov_list)
                        signal_to_noise_sample = signal_to_noise(sample_cov_list)

                        ratio = median_cov_sample / median_cov_baseline
                        if ratio == 0:
                            ratio = 0.1
                        log2_ratio = round(math.log2(ratio), 3)
                        log2_ratios_case.append(log2_ratio)

                        plot_dict["position"].append(start)
                        plot_dict["sample"].append(sample.name)
                        plot_dict["log2_ratio"].append(log2_ratio)
                        plot_dict["class"].append("case")

                    median_log2ratio_case = np.median(log2_ratios_case)
                    s2n_case = signal_to_noise(log2_ratios_case)
                    median_log2ratio_controls = np.median(log2_ratios_controls)
                    s2n_controls = signal_to_noise(log2_ratios_controls)

                    print(
                        sample.name
                        + " "
                        + exon
                        + " "
                        + str(median_log2ratio_case)
                        + " "
                        + str(s2n_case)
                        + " "
                        + str(median_log2ratio_controls)
                        + " "
                        + str(s2n_controls)
                    )

                    output_png_name = ("{}_{}.png").format(
                        sample.name, "_".join(tmp[0:4])
                    )
                    output_png = str(Path(sample.sample_folder) / output_png_name)
                    # plot_single_exon(plot_dict, sample.name, exon, output_png)

                    os.remove(candidate_bed)
        f.close()
