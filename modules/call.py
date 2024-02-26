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
from scipy.special import logsumexp
import math 
import shutil
import seaborn as sns
import matplotlib.pyplot as plt


def export_cnv_calls_to_bed(sample_list, analysis_dict):
    """
    Join CNV calls and export a single file
    """
    all_calls_name = f'{analysis_dict["output_name"]}.all.calls.bed'
    all_calls_bed = str(Path(analysis_dict["output_dir"]) / all_calls_name)
    o = open(all_calls_bed, "w")
    o.write("sample\tchr\tstart\tend\tregion\tgc\tmap\tzscore\tn_regions\tlog2_ratio\tcopy_number\tprob_score\tx\tcnvtype\tstdev\n")

    for sample in sample_list:
        if sample.analyzable == "False":
            continue

        tags = ['SVTYPE', 'REGION', 'GC', 'MAP', 'NREGIONS', 'LOG2RATIO', 'CN', 'CNV_SCORE']
        original = sample.cnv_calls_bed
        target_name = f'{sample.name}.GRAPES2.cnv.bed'
        target = os.path.join(analysis_dict["output_dir"], sample.name, target_name)
        g = open(target, "w")
        with open(sample.cnv_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                # sample	chr	start	end	regions	n_regions	log2_ratio	copy_number	cnvtype						
                # sample1.simulated	chr3	38591811	38598775	exon	50.252	100	-7.113	5	-0.892	1	0.990636622781486	1	DEL	0.077

                tmp = line.split("\t")
                info = {
                    "SVTYPE": tmp[-2],
                    "REGION": tmp[3].replace(";", "_"),
                    "GC": tmp[4],
                    "MAP": tmp[5],
                    "ZSCORE": tmp[6],
                    "NREGIONS": tmp[7],
                    "LOG2RATIO": tmp[8],
                    "CN": tmp[9],
                    "CNV_SCORE": str(round(float(tmp[10]),3))
                }
                info_str = "IMPRECISE;" + ";".join(f"{k}={v}" for k, v in info.items())
                coordinates = f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}"
                g.write(coordinates +"\t" + info_str +"\n")
        f.close()
        g.close()

        with open(sample.cnv_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                # print(line)
                o.write(sample.name + "\t" + line + "\n")
        f.close()
    o.close()

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


def filter_single_exon_cnv(sample_list, upper_del_threshold, dup_threshold, analysis_dict):
    """ """
    # get sample indices from per base coverage bed
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
        logging.info(msg)       

        a = pybedtools.BedTool(sample.raw_single_exon_calls)
        b = pybedtools.BedTool(analysis_dict["normalized_per_base"])
        c = a.intersect(b, wa=True, wb=True, stream=True)

        filtered_single_cnv_name = f"{sample.name}.filtered.single.exon.calls.bed"
        filtered_single_cnv_file = os.path.join(sample.sample_folder, filtered_single_cnv_name)
        sample.add("filtered_single_exon_calls", filtered_single_cnv_file)
        o = open(filtered_single_cnv_file, "w")

        candidate_cnvs = {}
        coverage_list = []

        sample_names = [sample.name]
        for idx,s in enumerate(sample.references):
            if s[0] == sample.name:
                continue
            if idx > 10:
                break
            sample_names.append(s[0])
        columns_list = ["chr", "start", "end", "info"]
        columns_list.extend(sample_names)
        
        case_log_ratios = []
        controls_log_ratios = []

        for line in iter(c):
            line = str(line)
            line = line.rstrip()

            tmp_line = line.split("\t")
            cnv_call = "\t".join(tmp_line[0:12])

            #CNVCALL chrX	119590505	119590624	NM_002294_1_2;LAMP2	36.130001	100.0	1	0.802	1	60
            if not cnv_call in candidate_cnvs:
                df = pd.DataFrame(columns=columns_list)
                candidate_cnvs[cnv_call] = {}
                candidate_cnvs[cnv_call]["list_rows"] = []
                candidate_cnvs[cnv_call]["dataframe"] = df
                candidate_cnvs[cnv_call]["case_coverage"] = []
                candidate_cnvs[cnv_call]["controls_coverage"] = []
                candidate_cnvs[cnv_call]["case_ratios"] =  []
                candidate_cnvs[cnv_call]["control_ratios"] =  []
                candidate_cnvs[cnv_call]["case_median_ratio"] = ""

            case_coverage = float(tmp_line[header_dict[sample.name]+13])
            candidate_cnvs[cnv_call]["case_coverage"].append(case_coverage)

            samples_cov_dict = {}
            background_cov_list = []

            for idx,s in enumerate(sample.references):
                if "baseline" in s[0]:
                    continue
                if s[0] == sample.name:
                    continue
                if idx > 10:
                    break
                control_coverage = float(tmp_line[header_dict[s[0]]+13])
                candidate_cnvs[cnv_call]["controls_coverage"].append(control_coverage)
                samples_cov_dict[s[0]] = control_coverage
                background_cov_list.append(control_coverage)

            row_dict = {
                "chr": tmp_line[0], 
                "start": tmp_line[1], 
                "end": tmp_line[2], 
                "info": tmp_line[3]
            }

            mean_bg_coverage = np.median(background_cov_list)
            if mean_bg_coverage == 0:
                mean_bg_coverage = 0.001

            if case_coverage == 0:
                case_sample_ratio = -3
            else:
                case_sample_ratio  = math.log2(case_coverage/mean_bg_coverage)

            row_dict[sample.name] = case_sample_ratio
            case_log_ratios.append(case_coverage/mean_bg_coverage)
            candidate_cnvs[cnv_call]["case_ratios"].append(case_coverage/mean_bg_coverage)
            candidate_cnvs[cnv_call]["case_median_ratio"] = case_coverage/mean_bg_coverage

            log_out_list = []

            for control_sample in samples_cov_dict:
                background_cov_list = []
                for other in samples_cov_dict:
                    if control_sample == other:
                        continue
                    background_cov_list.append(samples_cov_dict[other])
                mean_bg_coverage = np.median(background_cov_list)

                if mean_bg_coverage == 0:
                    mean_bg_coverage = 0.001

                if samples_cov_dict[control_sample] == 0:
                    control_ratio = -3
                else:
                    control_ratio = math.log2(samples_cov_dict[control_sample]/mean_bg_coverage)

                log_out_list.append(f"{control_sample}, {control_ratio}")
                row_dict[control_sample] = control_ratio
                candidate_cnvs[cnv_call]["control_ratios"].append(control_ratio)
            candidate_cnvs[cnv_call]["list_rows"].append(row_dict)

        for cnv_call in candidate_cnvs:

            tmp_cnv_call = cnv_call.split("\t")
            coordinates = f"{tmp_cnv_call[0]}:{tmp_cnv_call[1]}-{tmp_cnv_call[2]}"
            variant_title = f"{coordinates} {tmp_cnv_call[3]} {tmp_cnv_call[-1]}"

            candidate_cnvs[cnv_call]["dataframe"] = \
                pd.DataFrame.from_records(candidate_cnvs[cnv_call]["list_rows"])
            plot_single_exon_cnv(candidate_cnvs[cnv_call]["dataframe"], sample, variant_title)

            median_cov_case =  np.median(candidate_cnvs[cnv_call]["case_coverage"])
            median_cov_controls = np.median(candidate_cnvs[cnv_call]["controls_coverage"])
            std_cov_controls = np.std(candidate_cnvs[cnv_call]["controls_coverage"])
            s2n_case = median_cov_case/std_cov_controls
            s2n_controls = median_cov_controls/std_cov_controls

            if median_cov_controls == 0:
                median_cov_controls = 0.001
            if median_cov_case == 0:
                signal_ratio = -3
            else:
                signal_ratio = round(math.log2(median_cov_case/median_cov_controls), 3)

            z_score = calculate_z_score(signal_ratio, 
                candidate_cnvs[cnv_call]["control_ratios"])

            states = [0, 1, 2, 3, 4]
            probs_list = []
            for state in states:
                if state == 0:
                    state = 0.01
                ratio = state / 2
                expected_depth = median_cov_controls * ratio
                # print(median_cov_case, median_cov_controls, state, expected_depth, std_cov_controls, mad_cov_controls)
                #prob = np.log(norm.pdf(median_cov_case, loc=expected_depth, scale=std_cov_controls))
                prob = norm.pdf(median_cov_case, loc=expected_depth, scale=std_cov_controls)
                probs_list.append(prob)

            max_log_prob = np.max(probs_list)

            probs = np.exp(probs_list - max_log_prob)
            probs /= np.sum(probs)
            error_probs = 1 - probs
            cn = tmp_cnv_call[-3]            

            Q = -10 * np.log10(error_probs)
            Q_rounded = np.round(Q)
            Q_capped = np.clip(Q_rounded, 0, 60)
            if int(cn) > len(probs):
                cn = str(len(probs))

            single_prob = probs[int(cn)]

            if signal_ratio <= upper_del_threshold or signal_ratio >= dup_threshold:
                if signal_ratio <= upper_del_threshold:
                    svtype = "DEL"
                else:
                    svtype = "DUP"
                # print(tmp_cnv_call, Q_capped, probs, error_probs, s2n_case, s2n_controls)
                # print(sample.name, cnv_call, "s2n_controls:", median_cov_controls/std_cov_controls,"std_cov_controls:", std_cov_controls, "s2n_case:",s2n_case, "s2n_controls:", s2n_controls, "case_ratio:", signal_ratio, "median_cov_case:", median_cov_case,"median_cov_controls:", median_cov_controls, "zscore:", z_score)
                if s2n_controls >= 5 and abs(z_score) > 2.5:
                    tmp_cnv_call[-2] = str(single_prob)
                    cnv_call = '\t'.join(tmp_cnv_call)
                    o.write(cnv_call+"\t"+svtype+"\n")
        o.close()


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
    Release a list of raw segmented calls and also single-exon cnvs
    """
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)
    for sample in sample_list:

        if sample.analyzable == "False":
            continue

        msg = f" INFO: Calling segmented CNVs on sample {sample.name}"
        logging.info(msg)

        seen_roi_dict = {}

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


        ratio_file_no_header = sample.ratio_file.replace(".bed", ".noheader.bed")
        with open(ratio_file_no_header, "w") as fh:
            with open(sample.ratio_file, "r") as rh:
                for line in rh:
                    if line.startswith("chr\tstart"):
                        continue
                    fh.write(line)
            rh.close()
        fh.close()

        rois_list = []
        o.write("chr\tstart\tend\tregions\tn_regions\tlog2_ratio\tcn\tphred\tzscore\tcnvtype\n")
        with open(sample.segment_file) as seg:
            for line in seg:
                line = line.rstrip("\n")
                # chr11	19200000	26550000	133	0.5013	DEL
                if line.startswith("chr\tstart"):
                    continue

                tmp = line.split("\t")
                chr = tmp[0]
                start = tmp[1]
                end = tmp[2]
                regions = tmp[3]
                n_regions = tmp[4]

                log2_ratio = float(tmp[5])
                cn = int(tmp[6])
                prob_score = float(tmp[7])
                mean_gc, mean_map = get_gc_map_from_segment(chr, start, end, ratio_file_no_header)

                zscore = round((log2_ratio-sample.mean_log2_ratio)/sample.std_log2_ratio, 3)
                cnvtype = ""

                if cn == 2:
                    continue

                if prob_score < 0.5:
                    continue

                # if prob_score < 15:
                #     continue
                if log2_ratio <= upper_del_threshold or log2_ratio >= dup_threshold:
                    pass
                else:
                    continue
                if cn > 2:
                    if log2_ratio >= dup_threshold and abs(zscore) > 2.5:
                        cnvtype = "DUP"
                else:
                    if log2_ratio <= upper_del_threshold and abs(zscore) > 2.5:
                        cnvtype = "DEL"
                if cnvtype != "":
                    outline = f"{line}\t1\t{cnvtype}\n"
                    if int(n_regions) > 1:
                        if analysis_dict["offtarget"] == False and  "pwindow" in regions:
                            continue

                        tmp_outline = outline.split("\t")
                        tmp_outline.insert(4, str(zscore))
                        tmp_outline.insert(4, str(mean_map))
                        tmp_outline.insert(4, str(mean_gc))
                        outline = "\t".join(tmp_outline)
                        o.write(outline)
                        p.write(outline)


                        tmp_rois_list = get_single_rois_from_segment(chr, start, end, ratio_file_no_header)
                        for coord in tmp_rois_list:
                            seen_roi_dict[coord] = coord

                        for roi in tmp_rois_list:
                            rois_list.append(roi)

                        # Gene is assumed to be located at offset 3
                        gene = regions                          
                        tmp_regions = re.split('[, ;]', regions)
                        if len(tmp_regions) > 1:
                            gene = tmp_regions[-1]
                        gene_plot, sample = plot_gene(sample.name, sample_list, gene, analysis_dict)
                    else:
                        coordinate = f"{chr}\t{start}\t{end}"
                        seen_roi_dict[coordinate] = coordinate
                        outline = f"{line}\t1\t{cnvtype}\n"

                        variant = '\t'.join(tmp[0:3])
                        gc = ratio_dict[variant]["gc"]
                        map = ratio_dict[variant]["map"]

                        if analysis_dict["offtarget"] == False and  "pwindow" in regions:
                            continue
                        get_gc_map_from_segment(chr, start, end, ratio_file_no_header)

                        tmp_outline = outline.split("\t")
                        tmp_outline.insert(4, str(zscore))
                        tmp_outline.insert(4, str(mean_map))
                        tmp_outline.insert(4, str(mean_gc))

                        outline = "\t".join(tmp_outline)

                        q.write(outline)
                        #o.write(outline)
                        #p.write(outline)

        seg.close()
        o.close()
        p.close()

        with open(sample.ratio_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\tstart"):
                    continue
                tmp = line.split("\t")
                chr = tmp[0]
                start = tmp[1]
                end = tmp[2]
                region = tmp[3]
                gc = tmp[4]
                map = tmp[5]
                coordinate = f"{chr}\t{start}\t{end}"
                if coordinate in seen_roi_dict:
                    continue
                if analysis_dict["offtarget"] and  "pwindow" in region:
                    continue

                #chr1	26378363	26378374	NM_032588_8_9;TRIM63	45.450001	100.0	0.013
                log2_ratio = float(tmp[-1])
                fold_change = 2 ** (log2_ratio)

                zscore = round((log2_ratio-sample.mean_log2_ratio)/sample.std_log2_ratio,3)

                cn = str(int( (fold_change * 2)+.5))
                tmp_list = [chr, start, end, region, gc, map, str(zscore), "1", str(log2_ratio), cn, "60"]

                if log2_ratio <= upper_del_threshold and abs(zscore) > 2.5:
                    cnvtype = "DEL"
                    tmp_list.append(cn)
                    tmp_list.append(cnvtype)
                    q.write("\t".join(tmp_list) + "\n")
                if log2_ratio >= dup_threshold and abs(zscore) > 2.5:
                    cnvtype = "DUP"
                    tmp_list.append(cn)
                    tmp_list.append(cnvtype)
                    q.write("\t".join(tmp_list) + "\n")
        f.close()
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

        msg = f" INFO: Calling CNVs on sample {sample.name}"
        logging.info(msg)

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
        o.write("chr\tstart\tend\tregions\tgc\tmap\tz_score\tn_regions\tlog2_ratio\tcopy_number\tscore\tx\tcnvtype\tstd\n")
        for variant in segmented_cnvs_dict:
            arr = np.array(segmented_cnvs_dict[variant]["ratios"])
            std = round(np.std(arr), 3)
            # if std >= 0.4:
            #     continue
            # arr = np.array(segmented_cnvs_dict[variant]["gc"])
            # gc_mean = round(np.mean(arr), 3)

            # arr = np.array(segmented_cnvs_dict[variant]["map"])
            # map_mean = round(np.mean(arr), 3)

            tmp_variant = variant.split("\t")
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
