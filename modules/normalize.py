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
from .baseline_db import calculate_bed_md5, import_baselines_to_df
import statsmodels.api as sm
from modules.loopca import simulate_artificial_cnvs, loo_pca, find_optimal_pc_removal



def normalize_read_depth(bed_file, gc_col='gc', num_gc_bins=50, scale_factor=1e6, meta_columns=None):
    """
    Normalize read depth by library size and correct for GC bias.
    The BED file is split into autosomal and sex chromosomes. Sample columns are auto-detected
    by excluding meta columns.

    Parameters:
      bed_file      : Path to the input BED file.
      gc_col        : Column name for GC content (should be numeric).
      num_gc_bins   : Number of bins to use when grouping by GC content.
      scale_factor  : Scaling factor for library size normalization (e.g., 1e6 for counts-per-million).
      meta_columns  : List of columns that are metadata (e.g. ['chr', 'start', 'end', 'name', 'gc']).
                      If None, a default list is used.

    Returns:
      autosomal_df: DataFrame with normalized data for autosomal chromosomes (chr1–chr22).
      sex_df      : DataFrame with normalized data for sex chromosomes (chrX, chrY).
    """
    # Define default meta columns if none are provided.
    if meta_columns is None:
        meta_columns = ['chr', 'start', 'end', 'name', gc_col]
    
    # Load the BED file into a DataFrame.
    df = pd.read_csv(bed_file, sep='\t')
    
    # Ensure the 'chr' column is a string and in lowercase (to standardize comparisons).
    df['chr'] = df['chr'].astype(str).str.lower()
    
    # Automatically determine sample columns: those not in meta_columns.
    sample_columns = df[6:]
    
    # Split the DataFrame into autosomal and sex chromosomes.
    autosomal_chr = ['chr' + str(i) for i in range(1, 23)]
    sex_chr = ['chrx', 'chry']
    autosomal_df = df[df['chr'].isin(autosomal_chr)].copy()
    sex_df = df[df['chr'].isin(sex_chr)].copy()

    def normalize_df(df_sub):
        # ---- Step 1: Library Size Normalization ----
        # Compute the total reads per sample in this subset.
        library_sizes = {sample: df_sub[sample].sum() for sample in sample_columns}
        # Scale each sample to counts per the given scale_factor.
        for sample in sample_columns:
            df_sub[sample] = (df_sub[sample] / library_sizes[sample]) * scale_factor

        # ---- Step 2: GC Bias Correction ----
        # Create bins for GC content.
        # df_sub['gc_bin'] = pd.cut(df_sub[gc_col], bins=num_gc_bins)
        # # For each sample, compute the median normalized depth per GC bin and map it back.
        # for sample in sample_columns:
        #     expected_depth = df_sub.groupby('gc_bin')[sample].median()
        #     expected = df_sub['gc_bin'].map(expected_depth)
        #     # To avoid division by zero, clip the expected values.
        #     epsilon = 1e-6
        #     expected = expected.clip(lower=epsilon)
        #     df_sub[sample] = df_sub[sample] / expected
        # # Remove the temporary gc_bin column.
        # df_sub.drop(columns='gc_bin', inplace=True)
        return df_sub

    autosomal_df = normalize_df(autosomal_df)
    sex_df = normalize_df(sex_df)

    return autosomal_df, sex_df


def launch_normalization(sample_list, analysis_dict, ann_dict):
    """ """

    analysis_dict["bed_md5"] = calculate_bed_md5(analysis_dict["bed"])

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
    else:
        df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")

    if analysis_dict["use_baseline_db"]:
        df = import_baselines_to_df(analysis_dict, df)

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
        if sample.mean_coverage_X == 0:
            sample_stats[sample.name]["MEAN_COVERAGEX"] = 0.01

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


def loess_normalization(df, sample_name, cov_target, frac=0.3):
    """
    Normalize exon-level coverage for a given sample by fitting a LOWESS model
    to the relationship between GC content and library-normalized coverage.
    
    Parameters:
      df: DataFrame with at least columns 'gc' and the library-normalized coverage.
      sample_name: The sample name.
      cov_target: The column name (e.g. "{sample}_normalized_library") holding the library-normalized coverage.
      frac: The fraction of data used when estimating each y-value in LOWESS.
    
    Returns:
      A pandas Series of the GC-normalized coverage.
    """
    # Get GC content (assumed to be percent values; if fraction, you may need to adjust)
    gc_values = df['gc'].values
    coverage = df[cov_target].values
    # Fit LOWESS model
    lowess_fit = sm.nonparametric.lowess(coverage, gc_values, frac=frac, return_sorted=False)
    # Avoid division by zero
    lowess_fit[lowess_fit == 0] = 1e-8
    # Calculate a global median to keep overall scale
    global_median = np.median(coverage)
    normalized = (df[cov_target] / lowess_fit) * global_median
    return normalized

# def normalize_exon_level(input_bed, sample_list, fields):
#     """
#     Normalize exon-level coverage for each sample.
#     Steps:
#       1. Read the BED file and compute exon length.
#       2. (Optional) Compute row-wise median if needed.
#       3. Create GC bins (here we still use pd.cut, but then we apply LOESS on the full data).
#       4. For each sample:
#            - Compute library-normalized coverage using norm_by_lib (assumed to be defined).
#            - Use loess_normalization to adjust the library-normalized coverage by GC content.
#            - For each additional field (e.g. GC, map), you may further adjust normalization using
#              grouping by GC bin (if desired).
#       5. Return the DataFrame with a new column, e.g. "{sample}_normalized_final".
#     """
#     df = pd.read_csv(input_bed, sep="\t")
    
#     df['length'] = df['end'] - df['start']
#     # Optionally, filter very short exons
#     # df = df[df['length'] > 10]
    
#     sample_names = [sample.name for sample in sample_list]
    
#     # Calculate row-wise median coverage (if needed)
#     df['median_coverage'] = df[sample_names].median(axis=1)
#     # Optionally filter rows with low coverage, then drop the column
#     # df = df[df['median_coverage'] >= 30]
#     df = df.drop(columns=['median_coverage'])
    
#     fields_str = ",".join(fields)
#     msg = f" INFO: Normalizing exon level coverage by {fields_str}"
#     logging.info(msg)
    
#     # Create GC bins (if needed for additional normalization)
#     bins = np.arange(0, 110, 10)  # assuming GC content in percent
#     df['gc_bin'] = pd.cut(df['gc'], bins)
    
#     # For each sample, perform library normalization and then GC normalization
#     for sample in sample_list:
#         sample_lib_tag = f"{sample.name}_normalized_library"
#         # Compute the sample median from the raw coverage column
#         sample_median = df[sample.name].median()
#         median_cov_chrX = round(df[df["chr"] != "chrX"][sample.name].median(), 6)
    
#         # norm_by_lib should adjust for library size, etc. (assumed defined elsewhere)
#         df[sample_lib_tag] = df.apply(
#             norm_by_lib,
#             sample_median=sample_median,
#             chrX_median=median_cov_chrX,
#             sample_name=sample.name,
#             total_reads=sample.ontarget_reads,
#             axis=1,
#         )
    
#         # Now apply LOESS normalization for GC content.
#         normalized_final = f"{sample.name}_normalized_final"
#         df[normalized_final] = loess_normalization(df, sample.name, sample_lib_tag, frac=0.3)
    
#         # If additional normalization by other fields is needed, you could loop over fields here.
#         # For example, if you want to further adjust for GC or mapping score using group medians,
#         # you can perform that calculation here.
    
#     # (Optionally apply PCA/SVD noise reduction, etc.)
#     # df_svd = pca_noise_reduction(df)
#     # df = df_svd
#     # normalized_svd_bed = input_bed.replace("read.counts.bed", "normalized.svd.bed")
#     # df.to_csv(normalized_svd_bed, sep='\t', index=False)
    
#     return df

def normalize_exon_level(input_bed, sample_list, fields):
    """ """

    df = pd.read_csv(input_bed, sep="\t")

    df['length'] = df['end'] - df['start']
    # df = df[df['length'] > 10]

    sample_names = [sample.name for sample in sample_list ]

    # Calculate median coverage row-wise
    df['median_coverage'] = df[sample_names].median(axis=1)

    # Filter rows where median coverage is less than the threshold
    # df = df[df['median_coverage'] >= 30]

    # Drop the 'median_coverage' column as it's not needed anymore
    df = df.drop(columns=['median_coverage'])

    fields_str = ",".join(fields)
    msg = f" INFO: Normalizing exon level coverage by {fields_str}"
    logging.info(msg)

    bins = np.arange(0, 110, 10) # assuming GC content is a fraction

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

    # df_svd = pca_noise_reduction(df)
    # df = df_svd

    # normalized_svd_bed = input_bed.replace("read.counts.bed", "normalized.svd.bed")
    # df_svd.to_csv(normalized_svd_bed, sep='\t', index=False)
     
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
        if gc_bin_size > 30:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrX"]) / row[median_field_chrX],
                6,
            )
    elif row["chr"] == "chrY":
        if gc_bin_size > 30:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrY"]) / row[median_field_chrY],
                6,
            )
    else:
        if gc_bin_size > 30:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_autosomes"])
                / row[median_field_autosomes],
                6,
            )

    return norm_cov
