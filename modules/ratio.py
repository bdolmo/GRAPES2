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
from functools import reduce


def calculate_coverage_ratios(sample_list, analysis_dict, log2=True):
    """
    Calculate coverage ratios
    """
    df_list = []
    ratio_fields = []

    all_ratios_name = ("{}.ratios.bed").format(analysis_dict["output_name"])
    all_ratios = str(Path(analysis_dict["output_dir"]) / all_ratios_name)

    merged_df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")
    for sample in sample_list:

        msg = f" INFO: Calculating coverage ratios for sample {sample.name}"
        logging.info(msg)

        baseline_samples = []
        num = 0
        for oth in sample.references:
            if num > 10:
                break
            normalized_depth_tag = ("{}_normalized_final").format(oth[0])
            baseline_samples.append(normalized_depth_tag)
            num += 1

        ratio_file_name = ("{}.ratios.bed").format(sample.name)
        ratio_file = str(Path(sample.sample_folder) / ratio_file_name)
        sample.add("ratio_file", ratio_file)

        # if not os.path.isfile(ratio_file):
        if not os.path.isfile(ratio_file):

            sample_ratio = ("{}_ratio").format(sample.name)
            normalized_depth_tag = ("{}_normalized_final").format(sample.name)
            new_df = merged_df
            new_df[sample_ratio] = merged_df.apply(
                do_ratio_ref,
                baseline=baseline_samples,
                sample=normalized_depth_tag,
                axis=1,
            )

            # Now substract sample ratios
            new_df = new_df[["chr", "start", "end", "exon", "gc", "map", sample_ratio]]
            df_list.append(new_df)
            # Write dataframe as bed
            new_df.to_csv(ratio_file, sep="\t", mode="w", index=None)
        else:
            new_df = pd.read_csv(ratio_file, sep="\t")
            df_list.append(new_df)

    result = reduce(
        lambda df1, df2: pd.merge(
            df1, df2, on=["chr", "start", "end", "exon", "gc", "map"]
        ),
        df_list,
    )
    result.to_csv(all_ratios, sep="\t", mode="w", index=None)
    analysis_dict["all_ratios"] = all_ratios

    return sample_list, analysis_dict


def do_ratio_ref(row, baseline, sample, log2=True):
    """
    log2 ratio calculation
    """
    a = row[baseline].to_numpy()
    median_baseline = np.median(a)
    if median_baseline == 0:
        median_baseline = 0.01
    ratio = row[sample] / float(median_baseline)
    if ratio == 0:
        ratio = 0.01
    log2_ratio = round(math.log2(ratio), 3)
    return log2_ratio
