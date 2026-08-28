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
import pandas as pd

from .baseline_db import init_db, calculate_bed_md5, calculate_baseline_median_depth
from .utils import atomic_write_dataframe, update_stage_manifest

pd.options.mode.chained_assignment = None  # default='warn'


def launch_sample_clustering(sample_list, analysis_dict):
    """ """

    if analysis_dict["use_baseline_db"]:
        init_db(analysis_dict["baseline_db"])

    sample_list, analysis_dict = calculate_depth_correlation(sample_list, analysis_dict)

    sample_list = cluster_samples(
        analysis_dict["correlation_tsv"],
        sample_list,
        analysis_dict,
        min_correlation=float(analysis_dict.get("min_reference_correlation", 0.85)),
        min_refs=int(analysis_dict.get("min_reference_samples", 3)),
        max_refs=int(analysis_dict.get("max_reference_samples", 10)),
    )

    analysis_dict = create_heatmap(sample_list, analysis_dict)

    return sample_list, analysis_dict


def cluster_samples(corr_tsv, sample_list, analysis_dict, min_correlation=0.85, min_refs=3, max_refs=10):
    """Select the best references using a raw Spearman correlation threshold."""
    if not -1.0 <= min_correlation <= 1.0:
        raise ValueError("min_reference_correlation must be between -1 and 1")
    if min_refs < 1:
        raise ValueError("min_reference_samples must be at least 1")
    if max_refs < min_refs:
        raise ValueError("max_reference_samples must be >= min_reference_samples")

    similarity_df = pd.read_csv(corr_tsv, sep="\t", index_col=0)
    corr_dict = defaultdict(dict)
    for current_sample in similarity_df.index:
        candidates = []
        for reference_name, similarity in similarity_df.loc[current_sample].items():
            if reference_name == current_sample or not np.isfinite(similarity):
                continue
            raw_correlation = (2.0 * float(similarity)) - 1.0
            if raw_correlation < min_correlation:
                continue
            candidates.append(
                (reference_name, float(similarity), raw_correlation)
            )

        candidates.sort(key=lambda item: item[2], reverse=True)
        selected = candidates[:max_refs]
        corr_dict[current_sample]["n_references"] = len(selected)
        corr_dict[current_sample]["correlations"] = {
            name: round(similarity, 6) for name, similarity, _ in selected
        }
        corr_dict[current_sample]["mean_correlation"] = (
            round(float(np.mean([item[1] for item in selected])), 6)
            if selected
            else 0.0
        )
        corr_dict[current_sample]["mean_raw_correlation"] = (
            round(float(np.mean([item[2] for item in selected])), 6)
            if selected
            else 0.0
        )

    # Get all available non-overlapping baselines (group of samples with high correlation)
    seen_samples = []
    baselines = []
    sorted_samples = sorted(corr_dict.items(), key=lambda x: x[1]["n_references"], reverse=True)
    for sample, data in sorted_samples:
        baseline_samples = []
        seen_samples.append(sample)

        for baseline_member in data["correlations"]:
            if baseline_member in seen_samples:
                continue
            baseline_samples.append(baseline_member)
            seen_samples.append(baseline_member)
        if baseline_samples:
            baseline_samples.append(sample)
            baselines.append(baseline_samples)
    if analysis_dict["use_baseline_db"]:
        calculate_baseline_median_depth(analysis_dict, sample_list, baselines)

    # Gather data for reporting later
    for sample in sample_list:
        sample_name = str(sample.name)

        # print(sample_name, corr_dict[sample_name])
        nrefs = corr_dict[sample_name]["n_references"]
        if nrefs >= min_refs:
            sample.add("analyzable", "True")
            sample.analysis_json["analyzable"] = "True"
        else:
            sample.add("analyzable", "False")
            sample.analysis_json["analyzable"] = "False"
            logging.warning(
                " WARNING: Sample %s has only %d references with raw Spearman "
                ">= %.3f; at least %d are required",
                sample_name,
                nrefs,
                min_correlation,
                min_refs,
            )
        sample.add("mean_correlation", corr_dict[sample_name]["mean_correlation"])
        sample.analysis_json["mean_correlation"] = corr_dict[sample_name]["mean_correlation"]
        sample.analysis_json["mean_raw_reference_correlation"] = corr_dict[
            sample_name
        ]["mean_raw_correlation"]
        sample.add(
            "mean_raw_reference_correlation",
            corr_dict[sample_name]["mean_raw_correlation"],
        )
        sample.analysis_json["reference_count"] = nrefs
        ref_dict = list(corr_dict[sample_name]["correlations"].items())
        sample.add("references", ref_dict)

        if nrefs >= min_refs:
            logging.info(
                " INFO: Selected %d references for %s (mean raw Spearman %.3f)",
                nrefs,
                sample_name,
                corr_dict[sample_name]["mean_raw_correlation"],
            )

    return sample_list


def calculate_depth_correlation(sample_list, analysis_dict):
    """ """
    df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")
    newdf = df
    names_list = []

    normalized_final_cols = [col for col in df.columns if '_normalized_final' in col]
    for sample_tag in normalized_final_cols:
        sample_name = sample_tag.replace("_normalized_final", "")
        
        newdf[sample_name] = df[sample_tag]
        names_list.append(sample_name)
    data = newdf[names_list]
    raw_correlation = data.corr(method="spearman")
    similarity_correlation = (raw_correlation + 1) / 2

    for sample in sample_list:
        sample.analysis_json["correlation_matrix"] = similarity_correlation.to_json()

    correlation_tsv = str(Path(analysis_dict["output_dir"]) / "correlation.tsv")
    raw_correlation_tsv = str(
        Path(analysis_dict["output_dir"]) / "correlation.raw.tsv"
    )
    atomic_write_dataframe(
        similarity_correlation,
        correlation_tsv,
        sep="\t",
        index=True,
    )
    atomic_write_dataframe(
        raw_correlation,
        raw_correlation_tsv,
        sep="\t",
        index=True,
    )
    analysis_dict["correlation_tsv"] = correlation_tsv
    analysis_dict["raw_correlation_tsv"] = raw_correlation_tsv
    update_stage_manifest(
        analysis_dict["output_dir"],
        "sample_correlation",
        [correlation_tsv, raw_correlation_tsv],
        metadata={"sample_and_baseline_count": len(names_list)},
    )

    return sample_list, analysis_dict


def create_heatmap(sample_list, analysis_dict):
    """ """
    df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")
    names_list = []
    normalized_final_cols = [col for col in df.columns if '_normalized_final' in col]
    for sample_tag in normalized_final_cols:
        sample_name = sample_tag.replace("_normalized_final", "")
        df[sample_name] = df[sample_tag]
        names_list.append(sample_name)
    data = df[names_list]
    dat_corr = data.corr(method="spearman")

    sns.set_context("talk")
    sns.set(font_scale=1.4)
    try:
        heatmap = sns.clustermap(dat_corr, metric="correlation", cmap="coolwarm")
        correlation_plot = str(
            Path(analysis_dict["output_dir"]) / "correlation.heatmap.png"
        )
        heatmap.figure.savefig(correlation_plot, dpi=400)
    except:
        pass

    return analysis_dict
