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
import subprocess
from natsort import natsorted, index_natsorted, order_by_index
import pybedtools
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from adjustText import adjust_text

pd.options.mode.chained_assignment = None  # default='warn'


def launch_sample_clustering(sample_list, analysis_dict):
    """ """
    sample_list, analysis_dict = calculate_depth_correlation(sample_list, analysis_dict)

    sample_list = cluster_samples(analysis_dict["correlation_tsv"], sample_list)

    analysis_dict = create_heatmap(sample_list, analysis_dict)

    return sample_list, analysis_dict


def cluster_samples(corr_tsv, sample_list, min_correlation=0.91, min_refs=1):
    """ """
    n_line = 0
    header = []
    corr_dict = defaultdict(dict)
    with open(corr_tsv) as f:
        for line in f:
            n_line += 1
            line = line.rstrip("\n")
            tmp = line.split("\t")
            # header
            if n_line == 1:
                for sample in tmp:
                    if sample == "":
                        continue
                    header.append(sample)
            # data
            else:
                current_sample = tmp[0]
                corr_dict[current_sample] = defaultdict(dict)
                corr_dict[current_sample]["correlations"] = defaultdict(dict)
                idx = 0
                tmp_corr_dict = defaultdict(dict)
                tmp_corr_list = []
                n_refs = 0
                for corr in tmp[1:]:
                    sample = header[idx]
                    if sample != current_sample:
                        if float(corr) >= min_correlation:
                            tmp_corr_dict[sample] = str(round(float(corr), 3))
                            tmp_corr_list.append(round(float(corr), 3))
                            n_refs += 1
                    idx += 1
                corr_dict[current_sample]["n_references"] = n_refs
                corr_dict[current_sample]["correlations"] = tmp_corr_dict
                corr_dict[current_sample]["mean_correlation"] = round(
                    np.mean(tmp_corr_list), 3
                )

    for sample in sample_list:
        sample_name = str(sample.name)
        nrefs = corr_dict[sample_name]["n_references"]
        if nrefs >= min_refs:
            sample.add("analyzable", "True")
        else:
            sample.add("analyzable", "False")
        sample.add("mean_correlation", corr_dict[sample_name]["mean_correlation"])
        ref_dict = sorted(
            corr_dict[sample_name]["correlations"].items(),
            key=lambda item: item[1],
            reverse=True,
        )
        sample.add("references", ref_dict)

    sys.exit()

    return sample_list


def calculate_depth_correlation(sample_list, analysis_dict):
    """ """
    df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")
    newdf = df
    names_list = []
    for sample in sample_list:
        sample_tag = ("{}").format(sample.name)
        sample_tag = ("{}_normalized_final").format(sample.name)

        newdf[sample.name] = df[sample_tag]
        names_list.append(sample.name)
    data = newdf[names_list]
    dat_corr = data.corr(method="pearson")

    correaltion_tsv = str(Path(analysis_dict["output_dir"]) / "correlation.tsv")
    dat_corr.to_csv(correaltion_tsv, sep="\t", mode="w")
    analysis_dict["correlation_tsv"] = correaltion_tsv
    return sample_list, analysis_dict


def create_heatmap(sample_list, analysis_dict):
    """ """
    df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")
    names_list = []
    for sample in sample_list:
        sample_tag = ("{}_normalized_final").format(sample.name)
        df[sample.name] = df[sample_tag]
        names_list.append(sample.name)
    data = df[names_list]
    dat_corr = data.corr()

    sns.set_context("talk")
    sns.set(font_scale=1.4)
    heatmap = sns.clustermap(dat_corr, metric="correlation", cmap="coolwarm")
    correlation_plot = str(
        Path(analysis_dict["output_dir"]) / "correlation.heatmap.png"
    )
    heatmap.figure.savefig(correlation_plot, dpi=400)

    return analysis_dict
