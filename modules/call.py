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
from modules.utils import (
    atomic_output_path,
    remove_bed_header,
    signal_to_noise,
    update_stage_manifest,
    validate_delimited_file,
)
from modules.plot import plot_single_exon_cnv, plot_gene
from scipy.spatial import distance
from scipy.stats import nbinom, poisson, norm
from scipy.special import logsumexp, gammaln
import math 
import shutil
import seaborn as sns
import matplotlib.pyplot as plt


def compute_zscore(log2_ratio, sample):
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

    return zscore


def export_all_calls(sample_list, analysis_dict):
    min_cnv_quality = float(analysis_dict.get("min_cnv_quality", 0.5))
    all_calls_name = f'{analysis_dict["output_name"]}.all.calls.bed'
    all_calls_bed = str(Path(analysis_dict["output_dir"]) / all_calls_name)
    sample_calls = defaultdict(list)

    def _get_cnv_quality_from_info(info_field):
        fallback_score = None
        for token in info_field.split(";"):
            if token.startswith("CNV_QUALITY="):
                try:
                    return float(token.split("=", 1)[1])
                except ValueError:
                    return None
            if token.startswith("CNV_SCORE="):
                try:
                    fallback_score = float(token.split("=", 1)[1])
                except ValueError:
                    fallback_score = None
        return fallback_score

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

    with atomic_output_path(all_calls_bed) as temporary_calls_bed:
        with open(temporary_calls_bed, "w", encoding="utf-8") as o:
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
                        cnv_quality = _get_cnv_quality_from_info(tmp[3])
                        if cnv_quality is None or cnv_quality < min_cnv_quality:
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
                            "cnv_quality": cnv_quality,
                        })

        if not validate_delimited_file(
            temporary_calls_bed,
            required_columns=["sample", "chr", "start", "end", "info"],
            min_columns=5,
            min_data_rows=0,
        ):
            raise RuntimeError("Final cohort call table failed validation")

    update_stage_manifest(
        analysis_dict["output_dir"],
        "final_cohort_calls",
        [all_calls_bed],
    )

    msg = f" INFO: Final call summary (CNV_QUALITY >= {min_cnv_quality:.1f})"
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
                f"CNV_QUALITY={call['cnv_quality']:.4f}"
            )
            print(msg)
    if total_calls == 0:
        msg = " INFO: No calls passed the CNV_QUALITY threshold"
        print(msg)
    msg = f" INFO: Total calls\t{total_calls}"
    print(msg)


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


def _low_is_good_quality(value, good, poor, default=0.5):
    """Scale a non-negative dispersion metric to quality in [0, 1]."""
    try:
        metric = float(value)
    except (TypeError, ValueError):
        return default
    if not np.isfinite(metric) or metric < 0:
        return default
    return 1.0 - _scaled_noise(metric, good, poor)


def _reference_quality(raw_correlation):
    """Scale raw Spearman correlation to a bounded reference-quality value."""
    try:
        correlation = float(raw_correlation)
    except (TypeError, ValueError):
        return 0.5
    if not np.isfinite(correlation):
        return 0.5
    return _clamp01((correlation - 0.80) / (0.98 - 0.80))


def compute_sample_quality(log2_mad, pct_calls, raw_correlation):
    """Combine sample-level evidence once into an uncalibrated [0, 1] score."""
    dispersion_quality = _low_is_good_quality(log2_mad, 0.10, 0.35)
    call_burden_quality = _low_is_good_quality(pct_calls, 0.8, 6.0)
    reference_quality = _reference_quality(raw_correlation)
    return _clamp01(
        0.45 * dispersion_quality
        + 0.25 * call_burden_quality
        + 0.30 * reference_quality
    )


def _signal_fit_score(log2ratio, copy_number, tolerance):
    """Score agreement between observed log2 ratio and assigned copy number."""
    try:
        observed = float(log2ratio)
        copy_number = float(copy_number)
    except (TypeError, ValueError):
        return 0.0
    if not np.isfinite(observed) or not np.isfinite(copy_number):
        return 0.0
    expected = -2.5 if copy_number <= 0 else math.log2(copy_number / 2.0)
    difference = abs(observed - expected)
    return 1.0 / (1.0 + (difference / tolerance) ** 2)


def _posterior_quality(posterior):
    try:
        value = float(posterior)
    except (TypeError, ValueError):
        return 0.5
    if not np.isfinite(value):
        return 0.5
    return _clamp01(value)


def _roi_support_score(n_rois):
    """Smoothly reward multi-target support without forcing a perfect score."""
    try:
        count = max(int(n_rois), 1)
    except (TypeError, ValueError):
        count = 1
    return _clamp01(1.0 - np.exp(-count / 3.0))


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
    control_cv,
    cn,
    nRois,
    kl_divergence=None,
    per_base_log2_snr=None,
    per_base_state_prob=None,
    return_components=False,
):
    """Return an uncalibrated single-exon CNV quality score in [0, 1]."""
    signal_fit = _signal_fit_score(log2ratio, cn, tolerance=0.25)
    hmm_posterior = _posterior_quality(posterior_prob)
    control_dispersion = _low_is_good_quality(
        control_cv,
        good=0.05,
        poor=0.25,
    )
    sample_quality = compute_sample_quality(sampleStd, pctCalls, corr)
    per_base_support = _clamp01(
        0.50 * _posterior_quality(per_base_state_prob)
        + 0.25 * _kl_resolution_score(kl_divergence)
        + 0.25 * _single_exon_log2_snr_score(per_base_log2_snr)
    )
    quality = _clamp01(
        0.32 * signal_fit
        + 0.20 * hmm_posterior
        + 0.23 * per_base_support
        + 0.15 * control_dispersion
        + 0.10 * sample_quality
    )
    components = {
        "CNV_QUALITY": quality,
        "HMM_POSTERIOR": hmm_posterior,
        "SIGNAL_FIT": signal_fit,
        "DISPERSION_SCORE": control_dispersion,
        "SAMPLE_QUALITY": sample_quality,
        "ROI_SUPPORT": _roi_support_score(nRois),
        "PERBASE_SUPPORT": per_base_support,
    }
    return components if return_components else quality


def score_multiple_exon(
    svtype,
    posterior_prob,
    sampleStd,
    pctCalls,
    corr,
    log2ratio,
    event_std,
    cn,
    nRois,
    return_components=False,
):
    """Return an uncalibrated multi-exon CNV quality score in [0, 1]."""
    signal_fit = _signal_fit_score(log2ratio, cn, tolerance=0.35)
    hmm_posterior = _posterior_quality(posterior_prob)
    event_dispersion = _low_is_good_quality(
        event_std,
        good=0.05,
        poor=0.30,
    )
    sample_quality = compute_sample_quality(sampleStd, pctCalls, corr)
    roi_support = _roi_support_score(nRois)
    quality = _clamp01(
        0.30 * signal_fit
        + 0.25 * hmm_posterior
        + 0.20 * event_dispersion
        + 0.15 * roi_support
        + 0.10 * sample_quality
    )
    components = {
        "CNV_QUALITY": quality,
        "HMM_POSTERIOR": hmm_posterior,
        "SIGNAL_FIT": signal_fit,
        "DISPERSION_SCORE": event_dispersion,
        "SAMPLE_QUALITY": sample_quality,
        "ROI_SUPPORT": roi_support,
        "PERBASE_SUPPORT": None,
    }
    return components if return_components else quality

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
# Export final CNV calls with explicit quality components.
########################################################################

def export_cnv_calls_to_bed(sample_list, analysis_dict):
    """
    Join CNV calls and export a file with CNV_QUALITY and its components.
    For each sample, it reads the CNV calls bed file (sample.cnv_calls_bed), and for each call,
    CNV_SCORE is retained only as a compatibility alias of CNV_QUALITY.
    
    The expected input file is tab-separated with columns:
      [0] sample, [1] chr, [2] start, [3] end, [4] GC, [5] MAP, [6] ZSCORE, [7] NREGIONS, 
      [8] LOG2RATIO, [9] CN, [10] HMM_POSTERIOR, [11] auxiliary metrics,
      [12] SVTYPE, [13] CONTROL_CV, [14] EVENT_STD.
    """

    total_rois = get_total_lines_bed(analysis_dict["bed"])
    min_cnv_quality = float(analysis_dict.get("min_cnv_quality", 0.5))
    

    for sample in sample_list:
        if sample.analyzable == "False":
            continue

        total_calls = get_total_lines_bed(sample.cnv_calls_bed)

        pctCalls = 100*(total_calls/total_rois)
        sample_log2_mad = getattr(sample, "log2_mad", sample.std_log2_ratio)
        raw_reference_correlation = getattr(
            sample,
            "mean_raw_reference_correlation",
            (2.0 * float(sample.mean_correlation)) - 1.0,
        )

        # Output file name
        target_name = f'{sample.name}.GRAPES2.cnv.bed'
        target = os.path.join(analysis_dict["output_dir"], sample.name, target_name)
        with open(target, "w") as g, open(sample.cnv_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                tmp = line.split("\t")
                # The last column is within-event log2-ratio standard deviation.
                try:
                    event_std = abs(float(tmp[-1]))
                except (TypeError, ValueError):
                    event_std = 0.3
                if event_std >= 0.3:
                    continue

                svtype = tmp[-3]
                region_str = tmp[3].replace(";", "_")
                gc_val = tmp[4]
                map_val = tmp[5]
                zscore = float(tmp[6])
                nRois = int(tmp[7])
                log2ratio = float(tmp[8])
                hmm_posterior = float(tmp[10])
                kl_divergence = None
                per_base_log2_snr = None
                per_base_state_prob = None
                if nRois == 1 and len(tmp) > 11:
                    aux_metric = tmp[11]
                    if "|" in aux_metric:
                        aux_parts = aux_metric.split("|")
                        if len(aux_parts) == 3:
                            probability_str, kl_str, snr_str = aux_parts
                            try:
                                per_base_state_prob = float(probability_str)
                            except (TypeError, ValueError):
                                per_base_state_prob = None
                        else:
                            # Backward-compatible parsing of legacy KL|SNR files.
                            kl_str, snr_str = aux_parts[0:2]
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

                cn = float(tmp[9])
                try:
                    control_cv = abs(float(tmp[-2]))
                except (TypeError, ValueError):
                    control_cv = None

                if nRois > 1:
                    score_components = score_multiple_exon(
                        svtype,
                        hmm_posterior,
                        sample_log2_mad,
                        pctCalls,
                        raw_reference_correlation,
                        log2ratio,
                        event_std,
                        cn,
                        nRois,
                        return_components=True,
                    )
                else:
                    score_components = score_single_exon(
                        svtype,
                        hmm_posterior,
                        sample_log2_mad,
                        pctCalls,
                        raw_reference_correlation,
                        log2ratio,
                        control_cv,
                        cn,
                        nRois,
                        kl_divergence=kl_divergence,
                        per_base_log2_snr=per_base_log2_snr,
                        per_base_state_prob=per_base_state_prob,
                        return_components=True,
                    )
                cnv_quality = score_components["CNV_QUALITY"]
                reported_dispersion = (
                    control_cv if nRois == 1 and control_cv is not None else event_std
                )

                if cnv_quality < min_cnv_quality:
                    continue

                # Prepare the INFO field
                info = {
                    "QUALITY_MODEL": "GRAPES2_HEURISTIC_V2",
                    "SVTYPE": svtype,
                    "REGION": region_str,
                    "GC": gc_val,
                    "MAP": map_val,
                    "ZSCORE": tmp[6],
                    "CV": round(reported_dispersion, 4),
                    "NREGIONS": tmp[7],
                    "LOG2RATIO": tmp[8],
                    "CN": tmp[9],
                    "HMM_POSTERIOR": round(
                        score_components["HMM_POSTERIOR"], 4
                    ),
                    "SIGNAL_FIT": round(score_components["SIGNAL_FIT"], 4),
                    "DISPERSION_SCORE": round(
                        score_components["DISPERSION_SCORE"], 4
                    ),
                    "SAMPLE_QUALITY": round(
                        score_components["SAMPLE_QUALITY"], 4
                    ),
                    "ROI_SUPPORT": round(score_components["ROI_SUPPORT"], 4),
                    "PERBASE_SUPPORT": (
                        round(score_components["PERBASE_SUPPORT"], 4)
                        if score_components["PERBASE_SUPPORT"] is not None
                        else "."
                    ),
                    "EVENT_STD": round(event_std, 4),
                    "CONTROL_CV": (
                        round(control_cv, 4) if control_cv is not None else "."
                    ),
                    # CNV_SCORE is retained as a backward-compatible alias.
                    "CNV_QUALITY": round(cnv_quality, 4),
                    "CNV_SCORE": round(cnv_quality, 4),
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
            
            # Normalize Gaussian state likelihoods on the log scale.
            log_likelihoods = []
            for state in [0,1,2,3,4]:
                effective_ratio = 0.01 if state == 0 else state / 2.0
                expected_depth = median_cov_controls * effective_ratio
                log_likelihoods.append(
                    norm.logpdf(
                        median_cov_case,
                        loc=expected_depth,
                        scale=std_cov_controls,
                    )
                )
            probs = np.exp(np.asarray(log_likelihoods) - logsumexp(log_likelihoods))
            probs /= np.sum(probs)
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
                    # Preserve the HMM posterior and record per-base evidence
                    # separately instead of replacing one probability with another.
                    tmp_cnv_call[-1] = (
                        f"{single_prob:.6f}|{kl_val:.4f}|{snr_metric:.4f}"
                    )
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


def call_raw_cnvs(
    sample_list,
    analysis_dict,
    upper_del_threshold,
    dup_threshold,
    analyze_single_exon_cnv=True,
):
    """
    Release a list of raw segmented calls and optionally single-exon cnvs.
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

        filtered_single_cnv_name = f"{sample.name}.filtered.single.exon.calls.bed"
        filtered_single_cnv_file = str(Path(sample.sample_folder) / filtered_single_cnv_name)
        sample.add("filtered_single_exon_calls", filtered_single_cnv_file)

        # We will use plain dicts and bedtools stuff, instead of pandas
        o = open(raw_seg_calls_bed, "w")
        p = open(seg_calls_bed, "w")
        q = open(raw_single_cnv_file, "w") if analyze_single_exon_cnv else None
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
        o.write(
            "chr\tstart\tend\tregions\tn_regions\tlog2_ratio\tcn\t"
            "hmm_posterior\tzscore\tcnvtype\tdispersion\n"
        )
        
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
                zscore = compute_zscore(log2_ratio, sample)

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
                    outline = f"{line}\t.\t{cnvtype}\t.\n"
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

        if not analyze_single_exon_cnv:
            open(raw_single_cnv_file, "w").close()
            open(filtered_single_cnv_file, "w").close()
            logging.info(f" INFO: Skipping single-exon CNVs on sample {sample.name}")
            continue

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
                zscore = compute_zscore(log2_ratio, sample)

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
                            ".",
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
        o.write(
            "chr\tstart\tend\tregions\tgc\tmap\tz_score\tn_regions\t"
            "log2_ratio\tcopy_number\thmm_posterior\tperbase_metrics\t"
            "cnvtype\tcontrol_cv\tevent_std\n"
        )
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
