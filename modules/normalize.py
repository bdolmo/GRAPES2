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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.extmath import randomized_svd
from statsmodels.nonparametric.smoothers_lowess import lowess


def normalize_read_depth(bed_file, window_size=100):
    # Read the BED file into a pandas DataFrame
    df = pd.read_csv(bed_file, sep='\t')

    # Split DataFrame based on 'chr' column
    autosomal_df = df[df['chr'].str.lower().isin(['chr'+str(i) for i in range(1, 23)])]
    sex_df = df[df['chr'].str.lower().isin(['chrx', 'chry'])]


    # Define the samples columns
    samples = df.columns[6:]  # Assumes sample columns are from index 6 onwards

    # Define a normalization function
    def normalize_df(df):
        # Normalize for total sequencing depth
        df[samples] = df[samples].div(df[samples].sum(axis=1), axis=0)

        # Normalize for GC content using rolling median
        for sample in samples:
            # Create a DataFrame with 'gc' and the current sample
            gc_df = pd.DataFrame({'gc': df['gc'], 'depth': df[sample]})
            # Sort by 'gc' and then apply rolling median to 'depth'
            gc_df = gc_df.sort_values('gc')
            gc_df['depth'] = gc_df['depth'].rolling(window_size, min_periods=30).median()
            # Assign the normalized depths back to the original DataFrame in their original order
            df[sample] = gc_df.sort_values('gc', ascending=False)['depth'].values
            
        return df


    # Apply the normalization function to both DataFrames
    autosomal_df = normalize_df(autosomal_df)
    sex_df = normalize_df(sex_df)

    normalized_bed = bed_file.replace(".bed", ".test.norm.bed")
    normalized_df = pd.concat([autosomal_df, sex_df])
    normalized_df.to_csv(normalized_bed, sep='\t', index=False)
    return normalized_bed


def loess_normalization(df, sample_name, field):

    # Adding a small constant to avoid log(0)
    y = df[field] + 1e-10

    # Log-transform y
    y_log = np.log(y)

    # Perform LOESS smoothing on log-transformed y
    lowess_smoothed = lowess(y_log, df['gc'], frac=0.05)

    # Interpolate the LOESS-smoothed values at the points in x
    y_smoothed_interpolated = np.interp(df['gc'], lowess_smoothed[:, 0], lowess_smoothed[:, 1])
    df[f"{sample_name}_normalized_final"] = y_log - y_smoothed_interpolated

    # lowess_smoothed = lowess(df[field], df['gc'], frac=0.3)
    # df[f"{sample_name}_normalized_final"] = df[f"{sample_name}_normalized_final"] - lowess_smoothed[:, 1]
    shift = abs(df[f"{sample_name}_normalized_final"].min())
    df[f"{sample_name}_normalized_final"] += shift

    return df


def svd_noise_reduction(df, var_cutoff=0.9, max_components=25):
    # Iterate over each sample

    start_index = df.columns.get_loc('map') + 1
    end_index = df.columns.get_loc('gc_bin')
    depth_cols = [col for col in df.columns if col.endswith("_normalized_final")]
    # depth_cols = df.columns[start_index:end_index]
    # Center the data
    # print(depth_cols)

    colmeans = np.median(df[depth_cols], axis=0)
    centered = df[depth_cols]-colmeans

    # Compute the SVD
    U, S, V = randomized_svd(centered, n_components=max_components, n_oversamples=max_components)

    # Determine the number of components to keep
    prop_v = np.cumsum(S**2 / np.sum(S**2))


    n_components = np.nonzero(prop_v > var_cutoff)[0][0]
    n_components = 15

    X_projected = np.dot(centered, V[n_components:, :].T)
    X_reconstructed = np.dot(X_projected, V[n_components:, :])
    residuals = centered - X_reconstructed
    shift_value = abs(np.min(residuals))
    residuals += shift_value

    # Add the residuals to the dataframe
    #residual_cols = [col + '_normalized_final' for col in depth_cols]
    residual_cols = [col + '_normalized_svd' for col in depth_cols]


    df[residual_cols] = pd.DataFrame(residuals, index=df.index)

    return df


def launch_normalization(sample_list, analysis_dict, ann_dict):
    """ """
    for sample in sample_list:

        # Load a dataframe of raw coverage data
        df = pd.read_csv(analysis_dict["unified_raw_depth"], sep="\t")

        # Calculate median coverage
        median_depth = round(df[sample.name].median(), 6)
        sample.add("median_depth", median_depth)

    # Normalize exon-level coverage by GC-content
    normalized_depth_name = f"{analysis_dict['output_name']}.normalized.depth.bed"
    normalized_depth = str(Path(analysis_dict["output_dir"]) / normalized_depth_name)
    analysis_dict["normalized_depth"] = normalized_depth

    norm_factors = ["gc"]
    if not os.path.isfile(normalized_depth):       
        # normalize ontarget exon-level normalized depth
        df = normalize_exon_level(analysis_dict["unified_raw_depth"], sample_list, norm_factors)

        # Export ontarget exon-level normalized depth
        df.to_csv(normalized_depth, sep="\t", mode="w", index=None)


    normalized_offtarget_name = f"{analysis_dict['output_name']}.normalized.offtarget.bed"
    normalized_offtarget = str(Path(analysis_dict["output_dir"]) / normalized_offtarget_name)
    analysis_dict["normalized_offtarget"] = normalized_depth

    if analysis_dict["offtarget"]:

        if not os.path.isfile(normalized_offtarget):
            # normalize off-target normalized depth
            df = normalize_exon_level(analysis_dict["offtarget_raw_counts"], sample_list, norm_factors)

            # Export offtarget normalized depth
            df.to_csv(normalized_offtarget, sep="\t", mode="w", index=None)
    

    if os.path.isfile(normalized_depth) and os.path.isfile(normalized_offtarget):

        # Load the two BED files into pandas DataFrames
        bed1 = pd.read_csv(normalized_depth, sep="\t" )
        bed2 = pd.read_csv(normalized_offtarget, sep="\t")

        # Get column names from the first DataFrame
        col_names = bed1.columns.tolist()

        # Concatenate the DataFrames
        combined = pd.concat([bed1, bed2])
        # Reset index to ensure no duplicate indices
        combined.reset_index(drop=True, inplace=True)
        # print(combined)

        # Set column names for the combined DataFrame
        combined.columns = col_names
 
        # # Use natsorted to get indices that would sort by chromosome in a natural way, and sort start position within each chromosome
        combined['chr'] = combined['chr'].astype(str)
        combined['start'] = combined['start'].astype(int)

        combined = combined.reindex(index=order_by_index(combined.index, index_natsorted(combined['chr'])))
        combined = combined.sort_values(['chr', 'start'])


        # Write the sorted DataFrame to a new BED file
        normalized_all_name = f"{analysis_dict['output_name']}.normalized.all.bed"
        normalized_all = str(Path(analysis_dict["output_dir"]) / normalized_all_name)
        analysis_dict["normalized_all"] = normalized_all
        combined.to_csv(normalized_all, sep="\t", header=True, index=False)


    normalized_per_base_name = f"{analysis_dict['output_name']}.normalized.per.base.bed"
    normalized_per_base_file = str(
        Path(analysis_dict["output_dir"]) / normalized_per_base_name
    )
    normalized_per_base_pca_file = normalized_per_base_file.replace(".bed", ".pca.bed")

    analysis_dict["normalized_per_base"] = normalized_per_base_file
    analysis_dict["normalized_per_base_pca"] = normalized_per_base_pca_file


    # analysis_dict["normalized_depth"] = normalize_read_depth(analysis_dict["unified_raw_depth"])

    if not os.path.isfile(normalized_per_base_file):
        sample_list, analysis_dict = normalize_per_base(
            sample_list, analysis_dict, norm_factors
        )
    

    return sample_list, analysis_dict


def normalize_per_base(sample_list, analysis_dict, fields):
    """ """
    fields_str = ",".join(fields)
    msg = f" INFO: Normalizing per base coverage by {fields_str}"
    logging.info(msg)

    normalized_per_base_name = f"{analysis_dict['output_name']}.normalized.per.base.bed"
    normalized_per_base_file = str(
        Path(analysis_dict["output_dir"]) / normalized_per_base_name
    )
    analysis_dict["normalized_per_base"] = normalized_per_base_file

    sample_stats = defaultdict(dict)
    for sample in sample_list:
        sample_stats[sample.name]["MEAN_COVERAGE"] = sample.mean_coverage
        sample_stats[sample.name]["MEAN_COVERAGEX"] = sample.mean_coverage_X

    sample_idx = {}
    o = open(normalized_per_base_file, "w")
    with open(analysis_dict["per_base_coverage"]) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            chromosome = tmp[0]
            if line.startswith("#chr\tstart") or line.startswith("chr\tstart"):
                o.write(line + "\n")
                for i in range(6, len(tmp)):
                    sample_name = tmp[i]
                    sample_idx[i] = sample_name
            else:
                gc = int(float(tmp[4]))
                norm_list = []
                for i in range(6, len(tmp)):
                    sample_name = sample_idx[i]
                    raw_coverage = int(tmp[i])
                    if "X" in chromosome:
                        normalized_lib = str(
                            round(
                                raw_coverage
                                / sample_stats[sample_name]["MEAN_COVERAGEX"],
                                3,
                            )
                        )
                    else:
                        normalized_lib = str(
                            round(
                                raw_coverage
                                / sample_stats[sample_name]["MEAN_COVERAGE"],
                                3,
                            )
                        )
                    norm_list.append(normalized_lib)
                o.write("\t".join(tmp[0:6]) + "\t" + "\t".join(norm_list) + "\n")
    o.close()

    return sample_list, analysis_dict


def norm_by_lib(row, sample_median, chrX_median, sample_name, total_reads):
    """ """
    length = row["end"] - row["start"]
    if sample_median == 0:
        sample_median = 0.01

    if row["chr"] == "chrX" or row["chr"] == 23 or row["chr"] == "X":
        norm_lib = row[sample_name]/chrX_median
    else:
        norm_lib = row[sample_name]/sample_median

    return norm_lib


def normalize_exon_level(input_bed, sample_list, fields):
    """ """

    df = pd.read_csv(input_bed, sep="\t")

    df['length'] = df['end'] - df['start']
    df = df[df['length'] > 10]

    sample_names = [sample.name for sample in sample_list ]

    # Calculate median coverage row-wise
    df['median_coverage'] = df[sample_names].median(axis=1)

    # Filter rows where median coverage is less than the threshold
    df = df[df['median_coverage'] >= 30]

    # Drop the 'median_coverage' column as it's not needed anymore
    df = df.drop(columns=['median_coverage'])

    fields_str = ",".join(fields)
    msg = f" INFO: Normalizing exon level coverage by {fields_str}"
    logging.info(msg)

    bins = np.arange(0, 105, 5) # assuming GC content is a fraction

    # Create a new column 'gc_bin' using pd.cut
    df['gc_bin'] = pd.cut(df['gc'], bins)

    for sample in sample_list:
        sample_lib_tag = f"{sample.name}_normalized_library"

        sample_median = df[sample.name].median()
        median_cov_chrX = round(df[(df["chr"] != "chrX")][sample.name].median(), 6)

        df[sample_lib_tag] = df.apply(
            norm_by_lib,
            sample_median=sample_median,
            chrX_median=median_cov_chrX,
            sample_name=sample.name,
            total_reads=sample.ontarget_reads,
            axis=1,
        )

        cov_target = sample_lib_tag
        idx = 0     

        # df = loess_normalization(df, sample.name, sample_lib_tag)

        for field in fields:
            idx += 1
            field_int = f"{field}_integer"

            # Get integer field (gc or map) value
            df[field_int] = df[field].apply(int)

            median_field_cov_autosomes = f"median_{field}_cov_autosomes"
            median_field_cov_chrX = f"median_{field}_cov_chrX"
            median_field_cov_chrY = f"median_{field}_cov_chrY"

            median_cov_autosomes = round(
                df[(df["chr"] != "chrX") & (df["chr"] != "chrY")][cov_target].median(),
                6,
            )
            median_cov_chrX = round(df[(df["chr"] != "chrX")][cov_target].median(), 6)
            median_cov_chrY = round(df[(df["chr"] == "chrY")][cov_target].median(), 6)

            stats_dict = {
                "median_cov_autosomes": median_cov_autosomes,
                "median_cov_chrX": median_cov_chrX,
                "median_cov_chrY": median_cov_chrY
            }

            # Group coverage by field and calculate the median
            df[median_field_cov_autosomes] = (
                df[(df["chr"] != "chrX") & (df["chr"] != "chrY")]
                .groupby("gc_bin")[cov_target]
                .transform(
                    median_by_interval,
                    cov_target=cov_target,
                    median=median_cov_autosomes,
                )
            )
            df[median_field_cov_chrX] = (
                df[(df["chr"] == "chrX")]
                .groupby("gc_bin")[cov_target]
                .transform(
                    median_by_interval, cov_target=cov_target, median=median_cov_chrX
                )
            )
            df[median_field_cov_chrY] = (
                df[(df["chr"] == "chrY")]
                .groupby("gc_bin")[cov_target]
                .transform(
                    median_by_interval, cov_target=cov_target, median=median_cov_chrY
                )
            )   

            #print(df.groupby("gc_bin").size().to_dict())
            gc_bin_dict = df.groupby("gc_bin").size().to_dict()

            # Apply rolling median
            normalized_field = f"{sample.name}_normalized_{field}"
            df[normalized_field] = df.apply(
                apply_normalization,
                cov_target=cov_target,
                stats_dict=stats_dict,
                field=field,
                gc_bin_dict=gc_bin_dict,
                axis=1,
            )
            if idx == len(fields):
                normalized_final = f"{sample.name}_normalized_final"
                df[normalized_final] = df[normalized_field]

    # df = svd_noise_reduction(df)
     
    return df


def median_by_interval(x, cov_target, median):
    """ """
    median_interval = round(x.median(), 3)
    if not median:
        median_interval = median

    return median_interval


def apply_normalization(row, cov_target, stats_dict, field, gc_bin_dict):
    """
    Normalization
    """
    median_field_chrX = f"median_{field}_cov_chrX"
    median_field_chrY = f"median_{field}_cov_chrY"
    median_field_autosomes = f"median_{field}_cov_autosomes"
    
    gc_bin = row["gc_bin"]
    gc_bin_size = 0
    
    try:
        gc_bin_dict[gc_bin]
    except:
        pass
    else:
        gc_bin_size = gc_bin_dict[gc_bin]

    norm_cov = row[cov_target]

    if row[median_field_chrY] == 0:
        row[median_field_chrY] = 0.01

    if row["chr"] == "chrX":
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrX"]) / row[median_field_chrX],
                6,
            )
    elif row["chr"] == "chrY":
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrY"]) / row[median_field_chrY],
                6,
            )
    else:
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_autosomes"])
                / row[median_field_autosomes],
                6,
            )

    return norm_cov
