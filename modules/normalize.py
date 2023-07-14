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



def apply_pca(bed_file_path, explained_variance_ratio=0.20):


    # Open the file and read the headers
    with open(bed_file_path, 'r') as bed_file:
        headers = [line for line in bed_file if line.startswith('#')]
    
    # Read the BED file into a pandas DataFrame
    df = pd.read_csv(bed_file_path, sep='\t', comment='#')

    # Keep only the columns corresponding to the samples
    samples_columns = df.columns[6:]
    df_samples = df[samples_columns]

    # Standardizing the features
    x = StandardScaler().fit_transform(df_samples)

    # Apply PCA
    pca = PCA()
    x_pca = pca.fit_transform(x)

    # Calculate the cumulative explained variance and find the number of components for the threshold
    cumulative_explained_variance = np.cumsum(pca.explained_variance_ratio_)
    num_components = np.argmax(cumulative_explained_variance >= explained_variance_ratio) + 1
    print(f"REMOVING {num_components}")

    # Create an array of zeros with the same shape as x_pca
    x_pca_zeroed = np.zeros_like(x_pca)

    # Copy the components that you want to remove to the zeroed array
    x_pca_zeroed[:, :num_components] = x_pca[:, :num_components]


    # Inverse transform the zeroed array and subtract the result from the original data
    x_reconstructed = x - pca.inverse_transform(x_pca_zeroed)


    # Scale the data to the range [0, 1]
    scaler = MinMaxScaler()
    x_scaled = scaler.fit_transform(x_reconstructed)

    # Convert scaled data into a DataFrame
    df_scaled = pd.DataFrame(data=x_scaled, columns=samples_columns)

    # Convert scaled data into a DataFrame
    df_scaled = pd.DataFrame(data=x_scaled, columns=samples_columns)

    # Replace the original sample columns with the deflated data
    df_final = df.drop(columns=samples_columns).join(df_scaled)


    normalized_pca_per_base = bed_file_path.replace(".bed", ".pca.bed")

    # Write the DataFrame to a new BED file
    with open(normalized_pca_per_base, 'w') as out_file:
        out_file.writelines(headers)
        df_final.to_csv(out_file, sep='\t', index=False)


def launch_normalization(sample_list, analysis_dict, ann_dict):
    """ """
    for sample in sample_list:

        # Load a dataframe of raw coverage data
        df = pd.read_csv(analysis_dict["unified_raw_depth"], sep="\t")

        # Calculate median coverage
        median_depth = round(df[sample.name].median(), 2)
        sample.add("median_depth", median_depth)

    # Normalize exon-level coverage by GC-content
    normalized_depth_name = f"{analysis_dict['output_name']}.normalized.depth.bed"
    normalized_depth = str(Path(analysis_dict["output_dir"]) / normalized_depth_name)
    analysis_dict["normalized_depth"] = normalized_depth

    norm_factors = ["gc"]
    if not os.path.isfile(normalized_depth):
        df = normalize_exon_level(analysis_dict, sample_list, norm_factors)
        df.to_csv(normalized_depth, sep="\t", mode="w", index=None)

    # Export exon level normalized coverage
    normalized_per_base_name = f"{analysis_dict['output_name']}.normalized.per.base.bed"
    normalized_per_base_file = str(
        Path(analysis_dict["output_dir"]) / normalized_per_base_name
    )
    normalized_per_base_pca_file = normalized_per_base_file.replace(".bed", ".pca.bed")

    analysis_dict["normalized_per_base"] = normalized_per_base_file
    analysis_dict["normalized_per_base_pca"] = normalized_per_base_pca_file


    # analysis_dict["normalized_depth"] = normalize_read_depth(analysis_dict["unified_raw_depth"])

    if not os.path.isfile(normalized_per_base_file) or not os.path.isfile(normalized_per_base_pca_file):
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

    # apply_pca(analysis_dict["normalized_per_base"])

    return sample_list, analysis_dict


def norm_by_lib(row, sample_name, total_reads):
    """ """
    length = row["end"] - row["start"]
    norm_by_len = round(10e4 * (row[sample_name] / length), 3)
    norm_lib = round(10e2 * (norm_by_len / total_reads), 3)
    return norm_lib


def normalize_exon_level(analysis_dict, sample_list, fields):
    """ """

    df = pd.read_csv(analysis_dict["unified_raw_depth"], sep="\t")

    fields_str = ",".join(fields)
    msg = f" INFO: Normalizing exon level coverage by {fields_str}"
    logging.info(msg)


    bins = np.arange(0, 105, 5) # assuming GC content is a fraction

    # Create a new column 'gc_bin' using pd.cut
    df['gc_bin'] = pd.cut(df['gc'], bins)


    for sample in sample_list:
        sample_lib_tag = f"{sample.name}_normalized_library"
        df[sample_lib_tag] = df.apply(
            norm_by_lib,
            sample_name=sample.name,
            total_reads=sample.ontarget_reads,
            axis=1,
        )

        cov_target = sample_lib_tag
        idx = 0     


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
                3,
            )
            median_cov_chrX = round(df[(df["chr"] != "chrX")][cov_target].median(), 3)
            median_cov_chrY = round(df[(df["chr"] == "chrY")][cov_target].median(), 3)

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

    print(gc_bin, "bin_size", gc_bin_size)

    norm_cov = row[cov_target]

    if row["chr"] == "chrX":
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrX"]) / row[median_field_chrX],
                2,
            )
    elif row["chr"] == "chrY":
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_chrY"]) / row[median_field_chrY],
                2,
            )
    else:
        if gc_bin_size > 10:
            norm_cov = round(
                (row[cov_target] * stats_dict["median_cov_autosomes"])
                / row[median_field_autosomes],
                2,
            )

    return norm_cov
