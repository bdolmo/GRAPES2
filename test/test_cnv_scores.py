def score_single_exon(svtype, sampleStd, pctCalls, corr, log2ratio, zscore, cn, nRois):
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
    # Sample-Level Metrics:
    stdMetric = 1 if sampleStd < 0.2 else 0
    pctCallsMetric = 1 if pctCalls < 1.5 else 0
    corrMetric = 1 if corr > 0.97 else 0
    SLM = 0.45 * stdMetric + 0.45 * pctCallsMetric + 0.1 * corrMetric

    # Convert observed log2 ratio to linear ratio.
    observed_ratio = 2 ** log2ratio
    # Expected ratio on a linear scale: For CN=2, expected=1; for CN=3, expected=1.5; etc.
    expected_ratio = cn / 2.0

    # Compute the absolute difference metric (scaled by 2.5)
    absDiff = 2.5 * abs(observed_ratio - expected_ratio)
    # print(observed_ratio, expected_ratio, svtype, sampleStd, pctCalls, corr, log2ratio, cn, nRois, s2n , s2nc)

    absDiffMetric = 1 - absDiff if (1 - absDiff) > 0 else 0

    signal2noiseMetric = 1 if s2n > 10 else 0
    signal2noiseControlsMetric = 1 if s2nc > 10 else 0

    if zscore > 3:
        zscoreMetric = 1
    elif zscore > 2 and zscore < 3:
        zscoreMetric = 0.5
    else:
        zscoreMetric = 0

    # RLM = 0.4 * absDiffMetric + 0.3 * signal2noiseMetric + 0.3 * signal2noiseControlsMetric
    RLM = 0.4 * absDiffMetric + 0.6*zscoreMetric

    score = 0.4 * SLM + 0.6 * RLM
    return score


def score_multiple_exon(svtype, sampleStd, pctCalls, corr, log2ratio, zscore, cn, nRois):
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
    if zscore > 3:
        zscoreMetric = 1
    elif zscore > 2:
        zscoreMetric = 0.5
    else:
        zscoreMetric = 0

    # For multiple exon calls, we incorporate a weight based on the number of regions.
    nRoisMetric = 1 if nRois > 1 else 0
    nRoiWeight = 0.2 + (nRois * 0.1)
    signalWeight = 1 - nRoiWeight

    # Combine the ROI-level evidence: here we blend the nRoiMetric with a weighted combination
    # of the absolute difference and zscore metrics.
    combinedROI = 0.4 * absDiffMetric + 0.6 * zscoreMetric
    RLM = (nRoiWeight * nRoisMetric) + (signalWeight * combinedROI)
    if RLM > 1:
        RLM = 1

    score = 0.4 * SLM + 0.6 * RLM
    return score
###############################################
# Example testing of the scoring functions:
###############################################
if __name__ == '__main__':
    # Example 1: Single exon deletion call:
    # For a deletion (CN = 1), expected linear ratio is 0.5, so expected log₂ ratio is -1.
    svtype = "DEL"
    sampleStd = 0.15    # low std, good quality
    pctCalls = 1.0      # good
    corr = 0.98         # good
    log2ratio = -1.0    # as expected for CN = 1
    zscore = 3.5        # high, so zscoreMetric=1
    s2n = 12            # defaults (not used now)
    s2nc = 12
    cn = 1
    nRois = 1          # single exon
    score1 = score_single_exon(svtype, sampleStd, pctCalls, corr, log2ratio, zscore, cn, nRois)
    print("Single exon score (deletion, CN=1):", score1)

    # Example 2: Single exon duplication call:
    # For duplication (CN = 3), expected linear ratio is 1.5 (log₂ ~0.585).
    svtype = "DUP"
    sampleStd = 0.15
    pctCalls = 1.0
    corr = 0.98
    log2ratio = 0.7    # slightly higher than expected; some noise
    zscore = 2.5       # in between, so zscoreMetric=0.5
    s2n = 12
    s2nc = 12
    cn = 3
    nRois = 1
    score2 = score_single_exon(svtype, sampleStd, pctCalls, corr, log2ratio, zscore, cn, nRois)
    print("Single exon score (duplication, CN=3):", score2)

    # Example 3: Multiple exon call:
    # For a multi-exon call (nRois = 5) with CN = 3.
    svtype = "DUP"
    sampleStd = 0.15
    pctCalls = 1.0
    corr = 0.98
    log2ratio = 0.7
    zscore = 2.5
    s2n = 12
    s2nc = 12
    cn = 3
    nRois = 5
    score3 = score_multiple_exon(svtype, sampleStd, pctCalls, corr, log2ratio, zscore, cn, nRois)
    print("Multiple exon score (duplication, CN=3, nRois=5):", score3)
