#!/usr/bin/env python3
"""
offpeak_master.py

This module simulates artificial CNVs in a master BED file (with normalized exon‐level coverage)
and uses a leave-one-out PCA (LOO-PCA) approach to determine the optimal number of principal components
to remove for noise reduction. The user specifies a total number of ROIs to modify.

Usage:
    python offpeak_master.py --bed master_normalized.bed --test_sample sample1_normalized_final --n_rois 100 --max_pc 10 --threshold 5

Author: Your Name
Date: YYYY-MM-DD
"""

import numpy as np
import pandas as pd
import random
import logging
from sklearn.decomposition import PCA
import argparse
import math
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


##############################################
# Normalization function (provided by user)
##############################################

# def normalize_exon_level(input_bed, sample_list, fields):
#     """
#     Read a master BED file with normalized exon-level coverage and perform additional processing.
#     (This function is as provided by your pipeline.)
#     """
#     df = pd.read_csv(input_bed, sep="\t")
#     df['length'] = df['end'] - df['start']
#     sample_names = [sample.name for sample in sample_list]

#     # Calculate median coverage row-wise and drop it.
#     df['median_coverage'] = df[sample_names].median(axis=1)
#     df = df.drop(columns=['median_coverage'])

#     fields_str = ",".join(fields)
#     logging.info(f"INFO: Normalizing exon level coverage by {fields_str}")

#     bins = np.arange(0, 110, 10)  # assuming GC content is a fraction
#     df['gc_bin'] = pd.cut(df['gc'], bins)

#     # Apply library normalization per sample (norm_by_lib and median_by_interval, apply_normalization are assumed to be defined)
#     for sample in sample_list:
#         sample_lib_tag = f"{sample.name}_normalized_library"
#         sample_median = df[sample.name].median()
#         median_cov_chrX = round(df[(df["chr"] != "chrX")][sample.name].median(), 6)
#         df[sample_lib_tag] = df.apply(
#             norm_by_lib,
#             sample_median=sample_median,
#             chrX_median=median_cov_chrX,
#             sample_name=sample.name,
#             total_reads=sample.ontarget_reads,
#             axis=1,
#         )
#         # Here you might also call loess_normalization, etc.
#         normalized_final = f"{sample.name}_normalized_final"
#         df[normalized_final] = df[sample_lib_tag]
#     return df


##############################################
# Artificial CNV Simulation
##############################################

def simulate_artificial_cnvs(df, test_sample, total_rois):
    """
    Randomly select 'total_rois' from the master BED DataFrame for the test sample
    and modify them:
      - Half receive 50% of their original coverage (simulate deletion)
      - Half receive 150% of their original coverage (simulate duplication)
      
    Parameters:
      df: pandas DataFrame from normalize_exon_level()
      test_sample: string, column name (e.g., "sample1_normalized_final")
      total_rois: integer, total number of ROIs to modify (must be even)
      
    Returns:
      modified_df: DataFrame with modifications applied to test_sample.
      modifications: dict mapping ROI index (DataFrame index) to 'DEL' or 'DUP'
    """
    # print(df.columns)

    # sys.exit()


    modified_df = df.copy()
    if total_rois % 2 != 0:
        raise ValueError("total_rois must be even.")
    indices = random.sample(list(df.index), total_rois)
    half = total_rois // 2
    del_indices = indices[:half]
    dup_indices = indices[half:]
    modifications = {}
    for idx in del_indices:
        original = modified_df.loc[idx, test_sample]
        modified_df.loc[idx, test_sample] = original * 0.5
        modifications[idx] = 'DEL'
    for idx in dup_indices:
        original = modified_df.loc[idx, test_sample]
        modified_df.loc[idx, test_sample] = original * 1.5
        modifications[idx] = 'DUP'
    return modified_df, modifications

##############################################
# Leave-One-Out PCA (LOO-PCA)
##############################################

def loo_pca(df, test_sample, n_components_removed):
    """
    Perform leave-one-out PCA on the master DataFrame.
    1. Fit PCA on all samples except the test_sample.
    2. Project the test sample onto the PCA space.
    3. Zero out the first n_components_removed scores.
    4. Reconstruct the test sample.
    
    Parameters:
      df: DataFrame with ROIs as rows and sample columns.
      test_sample: string, column name of the test sample.
      n_components_removed: int, number of top PCs to remove.
    
    Returns:
      reconstructed: numpy array of reconstructed coverage for the test sample.
    """
    training_df = df.drop(columns=[test_sample])
    pca = PCA()
    pca.fit(training_df)
    # Project test sample onto PC space.
    test_values = df[[test_sample]].values
    scores = pca.transform(test_values)
    # Zero out the first n_components_removed principal components.
    scores[:, :n_components_removed] = 0
    reconstructed = pca.inverse_transform(scores)
    return reconstructed.flatten()

##############################################
# Performance Evaluation
##############################################

def evaluate_performance(original, corrected, modifications, threshold):
    """
    Evaluate how many of the modified ROIs are recovered.
    For each modified ROI, if the absolute difference between the corrected and original
    coverage is at least 'threshold', count it as detected.
    
    Parameters:
      original: numpy array of original coverage for the test sample.
      corrected: numpy array of reconstructed coverage after LOO-PCA.
      modifications: dict mapping ROI indices to modification label.
      threshold: float, minimum absolute difference to call a change detected.
    
    Returns:
      fraction_detected: float, fraction of modified ROIs detected.
    """
    detected = 0
    total = len(modifications)
    for idx, mod in modifications.items():
        if abs(corrected[idx] - original[idx]) >= threshold:
            detected += 1
    return detected / total if total > 0 else 0

##############################################
# Optimize PC Removal
##############################################

def find_optimal_pc_removal(df, test_sample, modifications, max_pc_removal=10, threshold=5):
    """
    Iterate over candidate numbers of PCs to remove and evaluate performance.
    
    Parameters:
      df: master DataFrame with normalized coverage.
      test_sample: string, column name for the test sample.
      modifications: dict of ROI indices to modification type.
      max_pc_removal: int, maximum number of PCs to try.
      threshold: float, threshold for detecting a modified ROI.
      
    Returns:
      best_pc: int, the number of PCs to remove that gives the best performance.
      performance_dict: dict mapping number of PCs removed to performance.
    """
    original = df[test_sample].values
    performance_dict = {}
    best_perf = -1
    best_pc = 0
    for k in range(0, max_pc_removal+1):
        corrected = loo_pca(df, test_sample, n_components_removed=k)
        perf = evaluate_performance(original, corrected, modifications, threshold)
        performance_dict[k] = perf
        if perf > best_perf:
            best_perf = perf
            best_pc = k
    return best_pc, performance_dict

##############################################
# Main Module Execution
##############################################

def main():
    parser = argparse.ArgumentParser(description="Simulate artificial CNVs and optimize PC removal via LOO-PCA.")
    parser.add_argument("--bed", required=True, help="Master BED file with normalized exon-level coverage.")
    parser.add_argument("--test_sample", required=True, help="Column name for the test sample (e.g., sample1_normalized_final).")
    parser.add_argument("--n_rois", type=int, required=True, help="Total number of ROIs to modify (even number).")
    parser.add_argument("--max_pc", type=int, default=10, help="Maximum number of PCs to try removing.")
    parser.add_argument("--threshold", type=float, default=5, help="Minimum absolute difference threshold for detection.")
    args = parser.parse_args()

    # logging.info("Reading master BED file.")
    df = pd.read_csv(args.bed, sep="\t")
    
    # For demonstration, assume the master bed file has been produced by normalize_exon_level() 
    # and already contains columns such as sample1_normalized_final, etc.
    test_sample = args.test_sample
    # logging.info(f"Simulating artificial CNVs in {args.n_rois} ROIs for test sample {test_sample}.")
    modified_df, modifications = simulate_artificial_cnvs(df, test_sample, args.n_rois)
    
    best_pc, perf_dict = find_optimal_pc_removal(modified_df, test_sample, modifications, max_pc_removal=args.max_pc, threshold=args.threshold)
    # logging.info(f"Optimal number of PCs to remove: {best_pc}")
    # logging.info(f"Performance by number of PCs removed: {perf_dict}")
    
    # Optionally, write out the modified and corrected coverage.
    corrected = loo_pca(modified_df, test_sample, best_pc)
    modified_df[f"{test_sample}_corrected"] = corrected
    out_bed = args.bed.replace(".bed", ".offpeak_normalized.bed")
    modified_df.to_csv(out_bed, sep="\t", index=False)
    # logging.info(f"Output written to {out_bed}")

if __name__ == "__main__":
    main()
