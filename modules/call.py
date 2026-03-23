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
from modules.plot import plot_single_exon_cnv, plot_gene
from scipy.spatial import distance
from scipy.stats import nbinom, poisson, norm
from scipy.special import logsumexp, gammaln
import math 
import shutil
import seaborn as sns
import matplotlib.pyplot as plt


def compute_zscore_and_cv(log2_ratio, sample):
    sample_mean = getattr(sample, "mean_log2_ratio", np.nan)
    sample_std = getattr(sample, "std_log2_ratio", np.nan)

    if (
        np.isfinite(log2_ratio)
        and np.isfinite(sample_mean)
        and np.isfinite(sample_std)
        and sample_std > 0
    ):
        zscore = round((log2_ratio - sample_mean) / sample_std, 3)
    else:
        zscore = np.nan

    denom = log2_ratio + 0.01
    if (
        not np.isfinite(sample_std)
        or sample_std <= 0
        or not np.isfinite(denom)
        or abs(denom) < 1e-8
    ):
        cv = np.nan
    else:
        cv = sample_std / denom

    return zscore, cv


def export_all_calls(sample_list, analysis_dict):
    min_cnv_score = 0.5
    all_calls_name = f'{analysis_dict["output_name"]}.all.calls.bed'
    all_calls_bed = str(Path(analysis_dict["output_dir"]) / all_calls_name)
    sample_calls = defaultdict(list)

    def _get_cnv_score_from_info(info_field):
        for token in info_field.split(";"):
            if token.startswith("CNV_SCORE="):
                try:
                    return float(token.split("=", 1)[1])
                except ValueError:
                    return None
        return None

    def _parse_info_field(info_field):
        info_dict = {}
        for token in info_field.split(";"):
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            info_dict[key] = value
        return info_dict

    def _format_value(value):
        try:
            return f"{float(value):.4f}"
        except (TypeError, ValueError):
            return str(value)

    with open(all_calls_bed, "w") as o:
        o.write("sample\tchr\tstart\tend\tinfo\n")
        for sample in sample_list:
            if sample.analyzable == "False":
                continue
            with open(sample.calls_bed) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line.startswith("chr\tstart"):
                        continue
                    tmp = line.split("\t")
                    if len(tmp) < 4:
                        continue
                    cnv_score = _get_cnv_score_from_info(tmp[3])
                    if cnv_score is None or cnv_score < min_cnv_score:
                        continue
                    o.write(sample.name + "\t" + line + "\n")
                    info_dict = _parse_info_field(tmp[3])
                    sample_calls[sample.name].append({
                        "chr": tmp[0],
                        "start": tmp[1],
                        "end": tmp[2],
                        "svtype": info_dict.get("SVTYPE", "NA"),
                        "cn": info_dict.get("CN", "NA"),
                        "log2ratio": info_dict.get("LOG2RATIO", "NA"),
                        "cnv_score": cnv_score,
                    })

    msg = f" INFO: Final call summary (CNV_SCORE >= {min_cnv_score:.1f})"
    print(msg)
    total_calls = 0
    for sample in sample_list:
        if sample.analyzable == "False":
            continue
        calls = sample_calls.get(sample.name, [])
        n_calls = len(calls)
        if n_calls == 0:
            continue
        total_calls += n_calls
        msg = f" INFO: {sample.name}\t{n_calls} calls"
        print(msg)
        for call in calls:
            msg = (
                f" INFO: {call['chr']}:{call['start']}-{call['end']} "
                f"SVTYPE={call['svtype']} "
                f"CN={call['cn']} "
                f"LOG2RATIO={_format_value(call['log2ratio'])} "
                f"CNV_SCORE={call['cnv_score']:.4f}"
            )
            print(msg)
    if total_calls == 0:
        msg = " INFO: No calls passed the CNV_SCORE threshold"
        print(msg)
    msg = f" INFO: Total calls\t{total_calls}"
    print(msg)


def logarithmic_cv_score(coeff_variation, lower=0.05, upper=0.2, low_score=1.0, high_score=0):
    """
    Compute a CV score using logarithmic interpolation between low_score and high_score.
    
    Parameters:
      coeff_variation: the coefficient of variation (a positive number).
      lower: CV value corresponding to low_score (default 0.05).
      upper: CV value corresponding to high_score (default 0.2).
      low_score: score when CV is at or below lower (default 1.0, high quality).
      high_score: score when CV is at or above upper (default 0, low quality).
      
    Returns:
      A score (float) that decreases as CV increases, using a log-scale interpolation.
    """
    if coeff_variation <= lower:
        return low_score
    elif coeff_variation >= upper:
        return high_score
    else:
        # Compute the logarithms of the bounds and the observed CV.
        log_lower = np.log(lower)
        log_upper = np.log(upper)
        log_cv = np.log(coeff_variation)
        # Compute the fraction (0 when coeff_variation == lower, 1 when coeff_variation == upper)
        fraction = (log_cv - log_lower) / (log_upper - log_lower)
        # Linear interpolation between low_score and high_score in score-space:
        return low_score + fraction * (high_score - low_score)


def _single_exon_dispersion_penalty(sample_std, coeff_variation):
    """
    Extra penalty for highly dispersed single-exon calls.
    Values near/below the start thresholds keep penalty close to 1.0.
    Higher dispersion is down-weighted exponentially.
    """
    # Slightly stronger than the previous version.
    std_penalty = np.exp(-9.0 * max(0.0, sample_std - 0.10))
    cv_penalty = np.exp(-12.0 * max(0.0, coeff_variation - 0.10))
    # Use geometric mean to avoid over-penalizing when only one metric is mildly noisy.
    combined_penalty = np.sqrt(std_penalty * cv_penalty)
    return max(min(combined_penalty, 1.0), 0.18)


def _multiple_exon_cv_penalty(coeff_variation):
    """
    Additional CV penalty for multi-exon calls.
    This is milder than single-exon, but still down-weights noisy events.
    """
    cv_penalty = np.exp(-5.0 * max(0.0, coeff_variation - 0.12))
    return max(cv_penalty, 0.20)


def _clamp01(value):
    return min(max(value, 0.0), 1.0)


def _scaled_noise(value, low, high):
    """
    Scale a metric into [0,1], where 0 is low-noise and 1 is high-noise.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return 1.0
    if not np.isfinite(v):
        return 1.0
    if high <= low:
        return 1.0 if v > low else 0.0
    return _clamp01((v - low) / (high - low))


def compute_sample_noise(sample_std, pct_calls, mean_corr):
    """
    Build a sample-level noise index from global metrics.
    Returns a value in [0,1], where higher means noisier sample.
    """
    std_noise = _scaled_noise(sample_std, 0.10, 0.35)
    pct_calls_noise = _scaled_noise(pct_calls, 0.8, 6.0)
    try:
        corr_val = float(mean_corr)
    except (TypeError, ValueError):
        corr_val = 0.95
    if not np.isfinite(corr_val):
        corr_val = 0.95
    corr_deficit = max(0.0, 1.0 - corr_val)
    corr_noise = _scaled_noise(corr_deficit, 0.01, 0.08)
    return _clamp01(0.45 * std_noise + 0.25 * pct_calls_noise + 0.30 * corr_noise)


def _sample_noise_penalty(sample_noise, n_rois):
    """
    Convert sample-level noise to multiplicative score penalty.
    Penalize single-exon calls more strongly than multi-exon calls.
    """
    try:
        noise = float(sample_noise)
    except (TypeError, ValueError):
        noise = 0.5
    if not np.isfinite(noise):
        noise = 0.5
    noise = _clamp01(noise)
    is_single_exon = int(n_rois) == 1
    k = 2.2 if is_single_exon else 1.3
    floor = 0.18 if is_single_exon else 0.30
    return max(min(np.exp(-k * noise), 1.0), floor)


def _kl_resolution_score(kl_divergence):
    """
    Convert KL divergence to a bounded [0, 1] score.
    """
    if kl_divergence is None:
        return 0.5
    try:
        kl_value = float(kl_divergence)
    except (TypeError, ValueError):
        return 0.5
    if not np.isfinite(kl_value):
        return 0.5
    kl_value = max(0.0, kl_value)
    return min(1.0 - np.exp(-kl_value / 0.7), 1.0)


def _single_exon_log2_snr_score(per_base_log2_snr):
    """
    Normalize per-base log2-ratio SNR score into [0,1].
    Accepts either a pre-normalized metric or a raw SNR value.
    """
    if per_base_log2_snr is None:
        return 0.5
    try:
        snr_val = float(per_base_log2_snr)
    except (TypeError, ValueError):
        return 0.5
    if not np.isfinite(snr_val):
        return 0.5
    if 0.0 <= snr_val <= 1.0:
        return snr_val
    snr_val = max(0.0, snr_val)
    return min(1.0 - np.exp(-snr_val / 2.5), 1.0)


def score_single_exon(
    svtype,
    posterior_prob,
    sampleStd,
    pctCalls,
    corr,
    log2ratio,
    coeff_variation,
    cn,
    nRois,
    kl_divergence=None,
    per_base_log2_snr=None,
    sample_noise=0.0,
):
    """
    Compute a CNV score for a single-exon call based on multiple evidence,
    converting the observed log2 ratio to a linear ratio.
    
    Parameters:
      svtype: CNV type (unused here, but available for future adjustments)
      sampleStd: standard deviation metric for the sample (float)
      pctCalls: percent of calls metric (float)
      corr: correlation metric (float)
      log2ratio: observed log₂ ratio (float); e.g., 0 for normal, positive for duplications, negative for deletions.
      s2n: signal-to-noise metric for the sample (float)
      s2nc: signal-to-noise metric for controls (float)
      cn: copy number (float)
      nRois: number of regions (should be 1 for single-exon calls)
    
    Returns:
      A score (float) between 0 and 1.
    """
    # Sample-Level Metrics. Here we are heavy penalyzing samples with high standard deviation and low correlation
    stdMetric = 1 if sampleStd < 0.2 else 0
    pctCallsMetric = 1 if pctCalls < 1.5 else 0
    corrMetric = 1 if corr > 0.97 else 0
    SLM = 0.45 * stdMetric + 0.45 * pctCallsMetric + 0.1 * corrMetric

    # Major quality contributor for single-exon calls:
    # distance between observed and expected log2 ratio.
    expected_ratio = cn / 2.0
    if expected_ratio <= 0:
        expected_log2ratio = -2.5
    else:
        expected_log2ratio = math.log2(expected_ratio)
    log2_diff = abs(log2ratio - expected_log2ratio)
    absDiffMetric = 1.0 / (1.0 + (log2_diff / 0.2) ** 2)

    cvScore = logarithmic_cv_score(coeff_variation)
    try:
        posteriorMetric = float(posterior_prob)
    except (TypeError, ValueError):
        posteriorMetric = 0.5
    posteriorMetric = min(max(posteriorMetric, 0.0), 1.0)
    klMetric = _kl_resolution_score(kl_divergence)
    snrMetric = _single_exon_log2_snr_score(per_base_log2_snr)

    RLM = (
        0.52 * absDiffMetric
        + 0.15 * cvScore
        + 0.13 * posteriorMetric
        + 0.10 * klMetric
        + 0.10 * snrMetric
    )

    score = 0.15 * SLM + 0.85 * RLM
    score *= _single_exon_dispersion_penalty(sampleStd, coeff_variation)
    score *= _sample_noise_penalty(sample_noise, nRois)
    score = min(max(score, 0.0), 1.0)
    return score


def score_multiple_exon(svtype, posterio_prob, sampleStd, pctCalls, corr, log2ratio, coeff_variation, cn, nRois, sample_noise=0.0):
    """
    Compute a CNV score for a multiple-exon call based on multiple evidence.
    The observed log2 ratio is converted to a linear ratio.
    
    Parameters are similar to score_single_exon, but nRois should be >1.
    
    Returns:
      A score (float) between 0 and 1.
    """
    # Sample-Level Metrics:
    stdMetric = 1 if sampleStd < 0.2 else 0
    pctCallsMetric = 1 if pctCalls < 1.5 else 0
    corrMetric = 1 if corr > 0.97 else 0
    SLM = 0.45 * stdMetric + 0.45 * pctCallsMetric + 0.1 * corrMetric

    # Convert observed log2 ratio to linear ratio.
    observed_ratio = 2 ** log2ratio
    expected_ratio = cn / 2.0

    if observed_ratio > expected_ratio:
        absDiff = 2.5 * abs(observed_ratio - expected_ratio)
    else:
        absDiff = 2.5 * abs(expected_ratio - observed_ratio)
    absDiffMetric = 1 - absDiff if (1 - absDiff) > 0 else 0

    # Compute a z-score metric as before.
    cvScore = logarithmic_cv_score(coeff_variation)

    # For multiple exon calls, we incorporate a weight based on the number of regions.
    nRoisMetric = 1 if nRois > 1 else 0
    nRoiWeight = 0.2 + (nRois * 0.1)
    signalWeight = 1 - nRoiWeight

    # Combine the ROI-level evidence: here we blend the nRoiMetric with a weighted combination
    # of the absolute difference and zscore metrics.
    combinedROI = 0.5 * absDiffMetric + 0.5 * cvScore
    RLM = (nRoiWeight * nRoisMetric) + (signalWeight * combinedROI)
    if RLM > 1:
        RLM = 1

    score = 0.4 * SLM + 0.6 * RLM
    score *= _multiple_exon_cv_penalty(coeff_variation)
    score *= _sample_noise_penalty(sample_noise, nRois)
    score = min(max(score, 0.0), 1.0)
    return score

#########################################
def get_total_lines_bed(input_bed):
    """ """
    n_rois = 0
    with open(input_bed) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if line.startswith("chrom\tstart"):
                continue
            if line.startswith("chr\tstart"):
                continue
            n_rois+=1
    return n_rois

########################################################################
# export_cnv_calls_to_bed function with new CNV_SCORE calculation
########################################################################

def export_cnv_calls_to_bed(sample_list, analysis_dict):
    """
    Join CNV calls and export a single file with an updated CNV_SCORE computed from multiple evidence.
    For each sample, it reads the CNV calls bed file (sample.cnv_calls_bed), and for each call,
    it computes a new CNV_SCORE using the evidence metrics. The new score is added to the INFO field.
    
    The expected input file is tab-separated with columns:
      [0] sample, [1] chr, [2] start, [3] end, [4] GC, [5] MAP, [6] ZSCORE, [7] NREGIONS, 
      [8] LOG2RATIO, [9] CN, [10] stdev, [11] ... , [second-to-last] SVTYPE
    (Adjust column indices as needed.)
    """

    total_rois = get_total_lines_bed(analysis_dict["bed"])
    min_cnv_score = 0.5
    

    for sample in sample_list:
        if sample.analyzable == "False":
            continue

        total_calls = get_total_lines_bed(sample.cnv_calls_bed)

        pctCalls = 100*(total_calls/total_rois)
        sample_noise = compute_sample_noise(sample.std_log2_ratio, pctCalls, sample.mean_correlation)

        # Output file name
        target_name = f'{sample.name}.GRAPES2.cnv.bed'
        target = os.path.join(analysis_dict["output_dir"], sample.name, target_name)
        with open(target, "w") as g, open(sample.cnv_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                tmp = line.split("\t")
                # Here, we assume the last column is stdev and we skip calls with high stdev.
                try:
                    stdev = float(tmp[-1])
                except:
                    stdev = 0
                if stdev >= 0.3:
                    continue

                # Extract fields from the line.
                # Adjust indices as needed for your file format.
                # For our scoring, we assume:
                #   - SVTYPE: column -2
                #   - REGION: column 3 (end coordinate string)
                #   - GC: column 4
                #   - MAP: column 5
                #   - ZSCORE: column 6 (we use this as sampleStd)
                #   - NREGIONS: column 7 (we use as nRois)
                #   - LOG2RATIO: column 8 (we convert this to ratio by 2^(value))
                #   - CN: column 9
                #   - For the other metrics, we use defaults:
                #         pctCalls: default 1.0
                #         corr: default 0.98
                #         s2n: default 12
                #         s2nc: default 12
                svtype = tmp[-3]
                region_str = tmp[3].replace(";", "_")
                gc_val = tmp[4]
                map_val = tmp[5]
                # Use ZSCORE as a proxy for sampleStd
                zscore = float(tmp[6])
                nRois = int(tmp[7])
                log2ratio = float(tmp[8])
                posterior_prob = float(tmp[10])
                kl_divergence = None
                per_base_log2_snr = None
                if nRois == 1 and len(tmp) > 11:
                    aux_metric = tmp[11]
                    if "|" in aux_metric:
                        kl_str, snr_str = aux_metric.split("|", 1)
                        try:
                            kl_divergence = float(kl_str)
                        except (TypeError, ValueError):
                            kl_divergence = None
                        try:
                            per_base_log2_snr = float(snr_str)
                        except (TypeError, ValueError):
                            per_base_log2_snr = None
                    else:
                        try:
                            kl_divergence = float(aux_metric)
                        except (TypeError, ValueError):
                            kl_divergence = None

                # Convert log2ratio to ratio: 2^(log2ratio)
                # ratio = 2 ** log2ratio
                cn = float(tmp[9])
                cv = round(float(tmp[-2]), 3)

                #sample11.simulated DUP 0.6973133908081831 0.0010204081632653062 0.992 0.573 7.255 3.0 12
                # print(sample.name, svtype, posterior_prob, sample.std_log2_ratio, pctCalls, sample.mean_correlation, log2ratio, zscore, cn, nRois)
                # Choose scoring function based on nRois (if >1, treat as multiple exon call)
                if nRois > 1:
                    cnv_score = score_multiple_exon(svtype, posterior_prob, sample.std_log2_ratio, 
                        pctCalls, sample.mean_correlation, log2ratio, cv, cn, nRois, sample_noise=sample_noise)
                else:
                    cnv_score = score_single_exon(svtype, posterior_prob, sample.std_log2_ratio, 
                        pctCalls, sample.mean_correlation, log2ratio, cv, cn, nRois,
                        kl_divergence=kl_divergence,
                        per_base_log2_snr=per_base_log2_snr,
                        sample_noise=sample_noise)

                if cnv_score < min_cnv_score:
                    continue

                # Prepare the INFO field
                info = {
                    "SVTYPE": svtype,
                    "REGION": region_str,
                    "GC": gc_val,
                    "MAP": map_val,
                    "ZSCORE": tmp[6],
                    "CV": cv,
                    "NREGIONS": tmp[7],
                    "LOG2RATIO": tmp[8],
                    "CN": tmp[9],
                    "CNV_SCORE": str(round(cnv_score, 3))
                }
                info_str = "IMPRECISE;" + ";".join(f"{k}={v}" for k, v in info.items())
               
                coordinates = "\t".join([tmp[0], tmp[1], tmp[2]])
                g.write(coordinates + "\t" + info_str + "\n")
    return sample_list


def calculate_z_score(case_coverage, background_coverage_list):
    # cnv_case is the log2 ratio for the CNV in your case sample (a single value)
    # cnv_background is a list/array of log2 ratios for the CNV in your background samples

    mean_background = np.median(background_coverage_list)
    std_dev_background = np.std(background_coverage_list)

    if case_coverage == 0:
        case_coverage = 0.01

    if std_dev_background == 0:
        std_dev_background = 0.01

    z_score = (case_coverage - mean_background) / std_dev_background

    return z_score


def calculate_mad_z_score(case_coverage, background_coverage_list):
    # cnv_case is the log2 ratio for the CNV in your case sample (a single value)
    # cnv_background is a list/array of log2 ratios for the CNV in your background samples

    median_background = np.median(background_coverage_list)
    mad_background = np.median(np.abs(background_coverage_list - median_background))
    mad_z_score = 0.6745 * (case_coverage - median_background) / mad_background

    return mad_z_score


def ratio_to_dict(ratio_file):

    ratios_dict = {}
    with open(ratio_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            tmp_line = line.split("\t")
            region = "\t".join(tmp_line[0:3])
            ratios_dict[region] = {}
            #chr1	26380309	26380529	TRIM63	43.639999	100.0	0.249
            ratios_dict[region]["gc"] = tmp_line[4]
            ratios_dict[region]["map"] = tmp_line[5]
    f.close()

    return ratios_dict


def compute_exon_cv(coverage_array):
    """
    Compute the coefficient of variation (CV) for an array of per-base coverage values.
    """
    mean_cov = np.mean(coverage_array)
    std_cov = np.std(coverage_array)
    if mean_cov == 0:
        return 0
    return std_cov / mean_cov


def kl_gaussian(mu1, sigma1, mu2, sigma2, epsilon=1e-10):
    """
    Compute the KL divergence between two univariate Gaussians:
         N(mu1, sigma1^2) and N(mu2, sigma2^2).
    """
    return np.log((sigma2+epsilon)/(sigma1+epsilon)) + (sigma1**2 + (mu1 - mu2)**2)/(2*(sigma2**2+epsilon)) - 0.5

def filter_single_exon_cnv(sample_list, upper_del_threshold, dup_threshold, analysis_dict):
    """
    Filter single-exon CNVs using a Gaussian-based likelihood model and include a KL divergence measure.
    We compute the KL divergence between the heterozygous deletion state (effective ratio = 0.5)
    and the diploid state (effective ratio = 1.0) as a proxy for target resolution.
    """
    # Build header dictionary from normalized per-base coverage file.
    header_dict = {}
    with open(analysis_dict["normalized_per_base"]) as f:
        header_line = f.readline().strip('\n')
        tmp = header_line.split("\t")
        for i in range(6, len(tmp)):
            sample_name = tmp[i]
            header_dict[sample_name] = i
            header_dict[i] = sample_name

    for sample in sample_list:
        if sample.analyzable == "False":
            continue
        msg = f" INFO: Calling single-exon CNVs on sample {sample.name}"
        print(msg)
        
        # Intersect candidate calls with normalized per-base coverage.
        a = pybedtools.BedTool(sample.raw_single_exon_calls)
        b = pybedtools.BedTool(analysis_dict["normalized_per_base"])
        c = a.intersect(b, wa=True, wb=True, stream=True)
        
        filtered_file = os.path.join(sample.sample_folder, f"{sample.name}.filtered.single.exon.calls.bed")
        sample.add("filtered_single_exon_calls", filtered_file)
        o = open(filtered_file, "w")
        
        candidate_cnvs = {}
        sample_names = [sample.name]
        for idx, s in enumerate(sample.references):
            if s[0] == sample.name:
                continue
            if idx > 10:
                break
            sample_names.append(s[0])
        columns_list = ["chr", "start", "end", "info"]
        columns_list.extend(sample_names)
        
        for line in iter(c):
            line = str(line).rstrip()
            tmp_line = line.split("\t")
            cnv_call = "\t".join(tmp_line[0:12])
            if cnv_call not in candidate_cnvs:
                candidate_cnvs[cnv_call] = {
                    "list_rows": [],
                    "case_coverage": [],
                    "controls_coverage": [],
                    "case_ratios": [],
                    "case_log2_ratios": [],
                    "control_ratios": [],
                    "case_median_ratio": "",
                    "cv": None
                }
            # Extract case coverage.
            case_coverage = float(tmp_line[header_dict[sample.name] + 13])
            candidate_cnvs[cnv_call]["case_coverage"].append(case_coverage)
            
            samples_cov_dict = {}
            background_cov_list = []
            for idx, s in enumerate(sample.references):
                if "baseline" in s[0]:
                    continue
                if s[0] == sample.name:
                    continue
                if idx > 10:
                    break
                control_coverage = float(tmp_line[header_dict[s[0]] + 13])
                candidate_cnvs[cnv_call]["controls_coverage"].append(control_coverage)
                samples_cov_dict[s[0]] = control_coverage
                background_cov_list.append(control_coverage)
            
            # Keep only finite control coverage values to avoid invalid mean/std operations.
            valid_background_cov = [x for x in background_cov_list if np.isfinite(x)]

            # Compute coefficient of variation (CV) for controls.
            if valid_background_cov:
                median_bg_cov = np.median(valid_background_cov)
                if np.isfinite(median_bg_cov) and median_bg_cov > 0:
                    cv = np.std(valid_background_cov) / median_bg_cov
                else:
                    cv = 1.0
            else:
                cv = 1.0
            candidate_cnvs[cnv_call]["cv"] = cv
            
            row_dict = {
                "chr": tmp_line[0],
                "start": tmp_line[1],
                "end": tmp_line[2],
                "info": tmp_line[3]
            }
            mean_bg_coverage = np.median(valid_background_cov) if valid_background_cov else 0.001
            if not np.isfinite(mean_bg_coverage) or mean_bg_coverage <= 0:
                mean_bg_coverage = 0.001
            if not np.isfinite(case_coverage) or case_coverage <= 0:
                case_sample_ratio = -3
                case_bg_ratio = 0.0
            else:
                case_bg_ratio = case_coverage / mean_bg_coverage
                case_sample_ratio = math.log2(case_bg_ratio)
            row_dict[sample.name] = case_sample_ratio
            candidate_cnvs[cnv_call]["case_ratios"].append(case_bg_ratio)
            candidate_cnvs[cnv_call]["case_log2_ratios"].append(case_sample_ratio)
            candidate_cnvs[cnv_call]["case_median_ratio"] = case_bg_ratio
            
            for control_sample in samples_cov_dict:
                bg_for_control = [
                    samples_cov_dict[other]
                    for other in samples_cov_dict
                    if other != control_sample and np.isfinite(samples_cov_dict[other])
                ]
                median_bg = np.median(bg_for_control) if len(bg_for_control) > 0 else 0.001
                if (
                    not np.isfinite(samples_cov_dict[control_sample])
                    or samples_cov_dict[control_sample] <= 0
                    or not np.isfinite(median_bg)
                    or median_bg <= 0
                ):
                    control_ratio = -3
                else:
                    control_ratio = math.log2(samples_cov_dict[control_sample] / median_bg)
                row_dict[control_sample] = control_ratio
                candidate_cnvs[cnv_call]["control_ratios"].append(control_ratio)
            candidate_cnvs[cnv_call]["list_rows"].append(row_dict)
        
        # Process candidate CNVs.
        for cnv_call in candidate_cnvs:
            tmp_cnv_call = cnv_call.split("\t")
            valid_case_cov = [x for x in candidate_cnvs[cnv_call]["case_coverage"] if np.isfinite(x)]
            valid_control_cov = [x for x in candidate_cnvs[cnv_call]["controls_coverage"] if np.isfinite(x)]
            if not valid_case_cov or not valid_control_cov:
                continue

            median_cov_case = np.median(valid_case_cov)
            median_cov_controls = np.median(valid_control_cov)
            std_cov_controls = np.std(valid_control_cov)
            if not np.isfinite(median_cov_controls) or median_cov_controls <= 0:
                median_cov_controls = 0.001
            if not np.isfinite(std_cov_controls) or std_cov_controls <= 0:
                std_cov_controls = 0.001

            if not np.isfinite(median_cov_case) or median_cov_case <= 0:
                signal_ratio = -3
            else:
                signal_ratio = round(math.log2(median_cov_case / median_cov_controls), 3)
            z_score = calculate_z_score(signal_ratio, candidate_cnvs[cnv_call]["control_ratios"])
            
            # --- Original Gaussian-based likelihood (for comparison) ---
            probs_list = []
            for state in [0,1,2,3,4]:
                effective_ratio = 0.01 if state == 0 else state / 2.0
                expected_depth = median_cov_controls * effective_ratio
                prob = norm.pdf(median_cov_case, loc=expected_depth, scale=std_cov_controls)
                probs_list.append(prob)
            max_log_prob = np.max(probs_list)
            probs = np.exp(np.array(probs_list) - max_log_prob)
            probs /= np.sum(probs)
            most_state = np.argmax(probs)
            most_prob = probs[most_state]
            cn = tmp_cnv_call[-3]            

            single_prob = probs[int(cn)]

            # Now compute KL divergence between heterozygous deletion and diploid states distributions
            # THe idea is to calculate the overlap between these distributions to account for how much uncertainity may be
            # For heterozygous deletion (CN=1): effective ratio = 0.5; for diploid (CN=2): ratio = 1.
            expected_depth_het = median_cov_controls * 0.5
            expected_depth_dip = median_cov_controls * 1.0
            # Assume the standard deviation is std_cov_controls for both.
            kl_val = kl_gaussian(expected_depth_het, std_cov_controls, expected_depth_dip, std_cov_controls)

            valid_case_log2 = [
                x for x in candidate_cnvs[cnv_call]["case_log2_ratios"]
                if np.isfinite(x) and x > -2.95
            ]
            valid_controls_log2 = [
                x for x in candidate_cnvs[cnv_call]["control_ratios"]
                if np.isfinite(x) and x > -2.95
            ]
            case_log2_snr = signal_to_noise(valid_case_log2) if valid_case_log2 else 0.0
            controls_log2_snr = signal_to_noise(valid_controls_log2) if valid_controls_log2 else 0.0
            case_snr_component = 1.0 - np.exp(-case_log2_snr / 2.5)
            contrast_component = case_log2_snr / (case_log2_snr + controls_log2_snr + 0.25)
            snr_metric = min(max(0.7 * case_snr_component + 0.3 * contrast_component, 0.0), 1.0)
            
            cv_value = candidate_cnvs[cnv_call].get("cv", 1.0)
            if signal_ratio <= upper_del_threshold or signal_ratio >= dup_threshold:
                svtype_final = "DEL" if signal_ratio <= upper_del_threshold else "DUP"
                # Apply all filtering criteria: robust z-score, low CV, and adequate KL divergence.
                if abs(z_score) >= 2.5 and cv_value <= 0.25:
                    # Here, you can either output the KL divergence as part of the call or use it to adjust a score.
                    # We'll output the call with the KL value.
                    tmp_cnv_call[-2] = str(single_prob)
                    tmp_cnv_call[-1] = f"{kl_val:.4f}|{snr_metric:.4f}"
                    final_call = '\t'.join(tmp_cnv_call)
                    o.write(final_call + "\t" + svtype_final + "\t" + str(cv_value) + "\n")
                    msg = (
                        f" INFO: {sample.name} median_cov_controls={median_cov_controls:.4f} "
                        f"std_dev_controls={std_cov_controls:.4f} CV={cv_value:.4f} "
                        f"KL_divergence={kl_val:.4f} cnv_call={cnv_call} "
                        f"signal_ratio={signal_ratio:.4f} zscore={z_score:.4f}"
                    )
                    print(msg)
        o.close()
    return sample_list


# def filter_single_exon_cnv(sample_list, upper_del_threshold, dup_threshold, analysis_dict):
#     """ """
#     # get sample indices from per base coverage bed
#     header_dict = {}
#     with open(analysis_dict["normalized_per_base"]) as f:
#         header_line = f.readline().strip('\n')
#         tmp = header_line.split("\t")
#         for i in range(6, len(tmp)):
#             sample_name = tmp[i]
#             header_dict[sample_name] = i
#             header_dict[i] = sample_name

#     for sample in sample_list:

#         if sample.analyzable == "False":
#             continue
#         msg = f" INFO: Calling single-exon CNVs on sample {sample.name}"
#         logging.info(msg)       

#         a = pybedtools.BedTool(sample.raw_single_exon_calls)
#         b = pybedtools.BedTool(analysis_dict["normalized_per_base"])
#         c = a.intersect(b, wa=True, wb=True, stream=True)

#         filtered_single_cnv_name = f"{sample.name}.filtered.single.exon.calls.bed"
#         filtered_single_cnv_file = os.path.join(sample.sample_folder, filtered_single_cnv_name)
#         sample.add("filtered_single_exon_calls", filtered_single_cnv_file)
#         o = open(filtered_single_cnv_file, "w")

#         candidate_cnvs = {}
#         coverage_list = []

#         sample_names = [sample.name]
#         for idx,s in enumerate(sample.references):
#             if s[0] == sample.name:
#                 continue
#             if idx > 10:
#                 break
#             sample_names.append(s[0])
#         columns_list = ["chr", "start", "end", "info"]
#         columns_list.extend(sample_names)
        
#         case_log_ratios = []
#         controls_log_ratios = []

#         for line in iter(c):
#             line = str(line)
#             line = line.rstrip()

#             tmp_line = line.split("\t")
#             cnv_call = "\t".join(tmp_line[0:12])

#             #CNVCALL chrX	119590505	119590624	NM_002294_1_2;LAMP2	36.130001	100.0	1	0.802	1	60
#             if not cnv_call in candidate_cnvs:
#                 df = pd.DataFrame(columns=columns_list)
#                 candidate_cnvs[cnv_call] = {}
#                 candidate_cnvs[cnv_call]["list_rows"] = []
#                 candidate_cnvs[cnv_call]["dataframe"] = df
#                 candidate_cnvs[cnv_call]["case_coverage"] = []
#                 candidate_cnvs[cnv_call]["controls_coverage"] = []
#                 candidate_cnvs[cnv_call]["case_ratios"] =  []
#                 candidate_cnvs[cnv_call]["control_ratios"] =  []
#                 candidate_cnvs[cnv_call]["case_median_ratio"] = ""
#                 # Save the CV with the candidate call for later filtering/penalty:
#                 candidate_cnvs[cnv_call]["cv"] = ""

#             case_coverage = float(tmp_line[header_dict[sample.name]+13])
#             candidate_cnvs[cnv_call]["case_coverage"].append(case_coverage)

#             samples_cov_dict = {}
#             background_cov_list = []

#             for idx,s in enumerate(sample.references):
#                 if "baseline" in s[0]:
#                     continue
#                 if s[0] == sample.name:
#                     continue
#                 if idx > 10:
#                     break
#                 control_coverage = float(tmp_line[header_dict[s[0]]+13])
#                 candidate_cnvs[cnv_call]["controls_coverage"].append(control_coverage)
#                 samples_cov_dict[s[0]] = control_coverage
#                 background_cov_list.append(control_coverage)

#             cv = compute_exon_cv(background_cov_list)
#             candidate_cnvs[cnv_call]["cv"] = cv
#             row_dict = {
#                 "chr": tmp_line[0], 
#                 "start": tmp_line[1], 
#                 "end": tmp_line[2], 
#                 "info": tmp_line[3]
#             }

#             mean_bg_coverage = np.median(background_cov_list)
#             if mean_bg_coverage == 0:
#                 mean_bg_coverage = 0.001

#             if case_coverage == 0:
#                 case_sample_ratio = -3
#             else:
#                 case_sample_ratio  = math.log2(case_coverage/mean_bg_coverage)

#             row_dict[sample.name] = case_sample_ratio
#             case_log_ratios.append(case_coverage/mean_bg_coverage)
#             candidate_cnvs[cnv_call]["case_ratios"].append(case_coverage/mean_bg_coverage)
#             candidate_cnvs[cnv_call]["case_median_ratio"] = case_coverage/mean_bg_coverage

#             log_out_list = []

#             for control_sample in samples_cov_dict:
#                 background_cov_list = []
#                 for other in samples_cov_dict:
#                     if control_sample == other:
#                         continue
#                     background_cov_list.append(samples_cov_dict[other])
#                 mean_bg_coverage = np.median(background_cov_list)

#                 if mean_bg_coverage == 0:
#                     mean_bg_coverage = 0.001

#                 if samples_cov_dict[control_sample] == 0:
#                     control_ratio = -3
#                 else:
#                     control_ratio = math.log2(samples_cov_dict[control_sample]/mean_bg_coverage)

#                 log_out_list.append(f"{control_sample}, {control_ratio}")
#                 row_dict[control_sample] = control_ratio
#                 candidate_cnvs[cnv_call]["control_ratios"].append(control_ratio)
#             candidate_cnvs[cnv_call]["list_rows"].append(row_dict)

#         for cnv_call in candidate_cnvs:

#             tmp_cnv_call = cnv_call.split("\t")
#             coordinates = f"{tmp_cnv_call[0]}:{tmp_cnv_call[1]}-{tmp_cnv_call[2]}"
#             variant_title = f"{coordinates} {tmp_cnv_call[3]} {tmp_cnv_call[-1]}"

#             candidate_cnvs[cnv_call]["dataframe"] = \
#                 pd.DataFrame.from_records(candidate_cnvs[cnv_call]["list_rows"])
#             # plot_single_exon_cnv(candidate_cnvs[cnv_call]["dataframe"], sample, variant_title)

#             median_cov_case =  np.median(candidate_cnvs[cnv_call]["case_coverage"])
#             median_cov_controls = np.median(candidate_cnvs[cnv_call]["controls_coverage"])
#             std_cov_controls = np.std(candidate_cnvs[cnv_call]["controls_coverage"])
#             s2n_case = median_cov_case/std_cov_controls
#             s2n_controls = median_cov_controls/std_cov_controls

#             cv_value = candidate_cnvs[cnv_call].get("cv", 0)

#             if median_cov_controls == 0:
#                 median_cov_controls = 0.001
#             if median_cov_case == 0:
#                 signal_ratio = -3
#             else:
#                 signal_ratio = round(math.log2(median_cov_case/median_cov_controls), 3)

#             z_score = calculate_z_score(signal_ratio, 
#                 candidate_cnvs[cnv_call]["control_ratios"])

#             states = [0, 1, 2, 3, 4]
#             probs_list = []
#             for state in states:
#                 if state == 0:
#                     state = 0.01
#                 ratio = state / 2
#                 expected_depth = median_cov_controls * ratio
#                 prob = norm.pdf(median_cov_case, loc=expected_depth, scale=std_cov_controls)
#                 probs_list.append(prob)

#             max_log_prob = np.max(probs_list)

#             probs = np.exp(probs_list - max_log_prob)
#             probs /= np.sum(probs)
#             error_probs = 1 - probs
#             cn = tmp_cnv_call[-3]            

#             epsilon = 1e-10
#             Q = -10 * np.log10(np.clip(error_probs, epsilon, None))

#             Q_rounded = np.round(Q)
#             Q_capped = np.clip(Q_rounded, 0, 60)
#             if int(cn) >= len(probs):
#                 cn = str(len(probs)-1)

#             single_prob = probs[int(cn)]
#             # print(cv_value)
#             if signal_ratio <= upper_del_threshold or signal_ratio >= dup_threshold:
#                 if signal_ratio <= upper_del_threshold:
#                     svtype = "DEL"
#                 else:
#                     svtype = "DUP"

#                 if abs(z_score) > 2.5 and float(cv_value) <= 0.15:
#                     print(sample.name, "Coeff.Variation:", cv_value, cnv_call, "s2n_controls:", median_cov_controls/std_cov_controls,"std_cov_controls:", std_cov_controls, "s2n_case:",s2n_case, "s2n_controls:", s2n_controls, "case_ratio:", signal_ratio, "median_cov_case:", median_cov_case,"median_cov_controls:", median_cov_controls, "zscore:", z_score)
#                     tmp_cnv_call[-2] = str(single_prob)
#                     cnv_call = '\t'.join(tmp_cnv_call)
#                     o.write(cnv_call+"\t"+svtype+"\t"+str(cv_value)+"\n")
#         o.close()


def get_gc_map_from_segment(chr, start, end, ratios_bed):
    """
    """
    tmp_bed = ratios_bed.replace(".bed", ".tmp.intersect.bed")
    o = open(tmp_bed, "w")
    o.write(f"{chr}\t{start}\t{end}\n")
    o.close()

    bed_out = tmp_bed.replace(".bed", ".bedout.bed")

    a = pybedtools.BedTool(ratios_bed)
    b = pybedtools.BedTool(tmp_bed)
    c = a.intersect(b)
    c.saveas(bed_out)

    df = pd.read_csv(bed_out, sep="\t", header=None, names=["chr", "start", "end", "info", "gc", "map", "ratio"])

    gc = round(df["gc"].mean(), 3)
    map = round(df["map"].mean(), 3)

    os.remove(bed_out)
    os.remove(tmp_bed)

    return gc, map


def get_single_rois_from_segment(chr, start, end, ratios_bed):
    """
    """
    tmp_bed = ratios_bed.replace(".bed", ".tmp.intersect.bed")
    o = open(tmp_bed, "w")
    o.write(f"{chr}\t{start}\t{end}\n")
    o.close()

    rois_list = []

    bed_out = tmp_bed.replace(".bed", ".bedout.bed")

    a = pybedtools.BedTool(ratios_bed)
    b = pybedtools.BedTool(tmp_bed)
    c = a.intersect(b)
    for roi in c:
        roi = str(roi).rstrip("\n")
        tmp_roi = roi.split("\t")
        coordinate = '\t'.join(tmp_roi[0:3])
        if not coordinate in rois_list:
            rois_list.append(coordinate)

    os.remove(tmp_bed)

    return rois_list


def call_raw_cnvs(sample_list, analysis_dict, upper_del_threshold, dup_threshold):
    """
    Release a list of raw segmented calls and also single-exon cnvs.
    Only keep single-exon CNVs that do NOT overlap any multi-exon CNVs from the same sample.
    """
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)

    for sample in sample_list:
        if sample.analyzable == "False":
            continue

        msg = f" INFO: Calling segmented CNVs on sample {sample.name}"
        logging.info(msg)

        # A place to record multi-exon intervals so we can exclude overlapping single-exon calls later
        multi_exon_intervals = {}  # dict of lists, keyed by chromosome

        seg_calls_name = f"{sample.name}.seg.calls.bed"
        seg_calls_bed = str(Path(sample.sample_folder) / seg_calls_name)
        sample.add("seg_calls_bed", seg_calls_bed)

        raw_seg_calls = f"{sample.name}.raw.seg.calls.bed"
        raw_seg_calls_bed = str(Path(sample.sample_folder) / raw_seg_calls)
        sample.add("raw_seg_calls", raw_seg_calls_bed)

        raw_single_cnv_name = f"{sample.name}.raw.single.exon.calls.bed"
        raw_single_cnv_file = str(Path(sample.sample_folder) / raw_single_cnv_name)
        sample.add("raw_single_exon_calls", raw_single_cnv_file)

        # We will use plain dicts and bedtools stuff, instead of pandas
        o = open(raw_seg_calls_bed, "w")
        p = open(seg_calls_bed, "w")
        q = open(raw_single_cnv_file, "w")
        ratio_dict = ratio_to_dict(sample.ratio_file)

        # Write a no-header version of ratio_file
        ratio_file_no_header = sample.ratio_file.replace(".bed", ".noheader.bed")
        with open(ratio_file_no_header, "w") as fh:
            with open(sample.ratio_file, "r") as rh:
                for line in rh:
                    if line.startswith("chr\tstart"):
                        continue
                    fh.write(line)

        # Prepare headers
        o.write("chr\tstart\tend\tregions\tn_regions\tlog2_ratio\tcn\tphred\tzscore\tcnvtype\tcv\n")
        
        # ------------------------------
        # Multi-exon calls
        # ------------------------------
        with open(sample.segment_file) as seg:
            for line in seg:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"): 
                    continue

                tmp = line.split("\t")
                chrom = tmp[0]
                start = int(tmp[1])
                end   = int(tmp[2])
                regions   = tmp[3]
                n_regions = int(tmp[4])

                log2_ratio = float(tmp[5])
                cn = int(tmp[6])
                prob_score = float(tmp[7])

                mean_gc, mean_map = get_gc_map_from_segment(chrom, start, end, ratio_file_no_header)
                zscore, cv = compute_zscore_and_cv(log2_ratio, sample)

                # Filters:
                if mean_gc < 20 or mean_gc > 80:
                    continue
                if (end - start) < 10:
                    continue
                if cn == 2:
                    continue
                if prob_score < 0.5:
                    continue

                # Threshold-based type determination
                cnvtype = ""
                if log2_ratio <= upper_del_threshold or log2_ratio >= dup_threshold:
                    if cn > 2:
                        if log2_ratio >= dup_threshold and abs(zscore) > 2:
                            cnvtype = "DUP"
                    else:
                        if log2_ratio <= upper_del_threshold and abs(zscore) > 2:
                            cnvtype = "DEL"

                # If we have a valid CNV
                if cnvtype != "":
                    # We consider multi-exon CNVs only if n_regions > 1
                    outline = f"{line}\t1\t{cnvtype}\t{str(cv)}\n"
                    if n_regions > 1:
                        # Check if it's on/off target
                        if analysis_dict["offtarget"] == False and "pwindow" in regions:
                            continue

                        # Insert mean_gc, mean_map, zscore in the output
                        tmp_outline = outline.split("\t")
                        tmp_outline.insert(4, str(zscore))
                        tmp_outline.insert(4, str(mean_map))
                        tmp_outline.insert(4, str(mean_gc))
                        outline = "\t".join(tmp_outline)

                        o.write(outline)
                        p.write(outline)

                        # Record interval in a dictionary for overlap checks
                        if chrom not in multi_exon_intervals:
                            multi_exon_intervals[chrom] = []
                        multi_exon_intervals[chrom].append((start, end))

                    # If you wanted to handle single-exon CNVs from this file, you could do so,
                    # but presumably we want only from segment_file_map below.
        
        o.close()
        p.close()

        # -------------------
        # Single-exon calls
        # -------------------
        with open(sample.segment_file_map) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")

                chrom      = tmp[0]
                start      = int(tmp[1])
                end        = int(tmp[2])
                regions    = tmp[3]
                gc_content = float(tmp[4])
                mapval     = float(tmp[5])
                log2_ratio = float(tmp[6])
                state      = tmp[-2]
                prob       = tmp[-1] 

                # If off-target is False and 'pwindow' in regions, skip
                if not analysis_dict["offtarget"] and "pwindow" in regions:
                    continue
                
                cn = int(state)
                zscore, cv = compute_zscore_and_cv(log2_ratio, sample)

                # cv = sample.std_log2_ratio/log2_ratio

                if gc_content < 20 or gc_content > 80:
                    continue

                # Require the single-exon log2 ratio to pass DEL or DUP thresholds.
                if not (log2_ratio <= upper_del_threshold or log2_ratio >= dup_threshold):
                    continue

                cnvtype = ""
                # Check duplication
                if cn > 2:
                    if log2_ratio >= dup_threshold and abs(zscore) > 2:
                        cnvtype = "DUP"
                # Check deletion
                elif cn < 2:
                    if log2_ratio <= upper_del_threshold and abs(zscore) > 2:
                        cnvtype = "DEL"
                # Only proceed if not state == 2
                if float(prob) < 0.7:
                    continue

                if cn != 2 and cnvtype:
                    # --------------------
                    # Overlap check
                    # --------------------
                    # If there's no multi-exon calls on this chromosome, or no intervals stored, 
                    # we don’t worry about overlap
                    no_overlap = True
                    if chrom in multi_exon_intervals:
                        for (mstart, mend) in multi_exon_intervals[chrom]:
                            # Simple overlap check: (start < mend) and (end > mstart)
                            if start < mend and end > mstart:
                                no_overlap = False
                                break
                    
                    if no_overlap:
                        # Write single-exon CNV if it does NOT overlap
                        # Format the line as needed:
                        tmp_out = [
                            chrom,
                            str(start),
                            str(end),
                            regions,
                            tmp[4],  # gc_content
                            tmp[5],  # mapval
                            str(zscore),
                            "1",
                            str(log2_ratio),
                            str(cn),
                            prob,
                            "1",
                            cnvtype
                        ]
                        q.write("\t".join(tmp_out) + "\n")

        q.close()

    return sample_list


def unify_raw_calls(sample_list):
    """
    Combine raw single-exon and segmented calls
    """

    for sample in sample_list:

        if sample.analyzable == "False":
            continue

        df_list = []
        bed_list = [sample.raw_seg_calls, sample.filtered_single_exon_calls]

        raw_calls_name = f"{sample.name}.raw.calls.bed"
        raw_calls_bed = str(Path(sample.sample_folder) / raw_calls_name)
        sample.add("raw_calls_bed", raw_calls_bed)
        o = open(raw_calls_bed, "w")
        for file in bed_list:
            with open(file) as f:
                for line in f:
                    o.write(line)
            f.close()
        o.close()
    return sample_list


def call_cnvs(sample_list, upper_del_threshold, dup_threshold, z_score):
    """
    Calling CNVs
    """
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)

    for sample in sample_list:
        
        if sample.analyzable == "False":
            continue

        cnv_calls_name = f"{sample.name}.calls.bed"
        cnv_calls_bed = str(Path(sample.sample_folder) / cnv_calls_name)
        sample.add("cnv_calls_bed", cnv_calls_bed)

        ratio_no_header = remove_bed_header(sample.ratio_file, "chr\tstart")
        seg_no_header = remove_bed_header(sample.segment_file, "chr\tstart")

        # We will use plain dicts instead of pandas
        tmp_calls = cnv_calls_bed.replace(".bed", ".tmp.bed")
        o = open(tmp_calls, "w")

        with open(sample.raw_calls_bed) as seg:
            for line in seg:
                line = line.rstrip("\n")
                # chr11	19200000	26550000 GENE_ANNOT	133	0.5013	DEL
                if line.startswith("chr\tstart"):
                    continue
                o.write(line+"\n")
        seg.close()
        o.close()

        segmented_cnvs_dict = get_segmented_cnvs(ratio_no_header, tmp_calls)

        o = open(cnv_calls_bed, "w")
        o.write("chr\tstart\tend\tregions\tgc\tmap\tz_score\tn_regions\tlog2_ratio\tcopy_number\tscore\tkl_snr\tcnvtype\tstd\tcv\n")
        for variant in segmented_cnvs_dict:
            arr = np.array(segmented_cnvs_dict[variant]["ratios"])
            std = round(np.std(arr), 3)
            tmp_variant = variant.split("\t")

            # Calculate coefficient of variation
            # cv = float(tmp_variant[8])/std
            # tmp_variant.insert(4, str(map_mean))
            # tmp_variant.insert(4, str(gc_mean))

            variant = '\t'.join(tmp_variant)

            outline = f"{variant}\t{str(std)}\n"
            o.write(outline)
        o.close()
        sample.add("ready_cnv_bed", cnv_calls_bed)

    return sample_list


def get_segmented_cnvs(ratio_no_header, tmp_calls):
    """ """
    calls_dict = defaultdict(list)
    a = pybedtools.BedTool(ratio_no_header)
    b = pybedtools.BedTool(tmp_calls)

    c = a.intersect(b, wo=True, stream=True)

    for line in iter(c):
        line = str(line)
        #chr3	57882599	57882659	NM_007159_13_14;SLMAP	45.0	100.0	-1.037	chr3	57850274	57882659	NM_007159_9_10;SLMAP,NM_007159_10_11;	5	-0.844	1	DEL	    60
        #chr11	2905233	    2905365	    CDKN1C	                68.18	50.0	-0.783	chr11	2905233	    2906720	    CDKN1C,CDKN1C	                        2	-0.663	1	32.0	DEL	132
        line = line.rstrip() 
        tmp = line.split("\t")
        variant = "\t".join(tmp[7:-1])
        if not variant in calls_dict:
            calls_dict[variant] = {
                "ratios": [],
                "gc": [],
                "map": []
            }
            calls_dict[variant]["ratios"].append(float(tmp[6]))
            calls_dict[variant]["gc"].append(float(tmp[4]))
            calls_dict[variant]["map"].append(float(tmp[5]))
        else:
            calls_dict[variant]["ratios"].append(float(tmp[6]))
            calls_dict[variant]["gc"].append(float(tmp[4]))
            calls_dict[variant]["map"].append(float(tmp[5]))
        
    return calls_dict
