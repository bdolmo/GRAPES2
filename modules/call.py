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
from modules.utils import remove_bed_header

def export_all_calls(sample_list, analysis_dict):
    '''
        Join CNV calls and export a single file
    '''
    all_calls_name = ("{}.all.calls.bed").format(analysis_dict['output_name'])
    all_calls_bed  = str(Path(analysis_dict['output_dir']) / all_calls_name)
    o = open(all_calls_bed, 'w')
    o.write("sample\tchr\tstart\tend\tregions\tlog2_ratio\tcn\tcnvtype\n")
    for sample in sample_list:
        with open(sample.cnv_calls_bed) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith('chr\tstart'):
                    continue
                o.write(sample.name+"\t"+line+"\n")
        f.close()
    o.close()
    return sample_list

def call_raw_single_exon_cnv(sample_list, upper_del_threshold, dup_threshold):
    '''
        Calling additional cnvs that do not overlap with segmented calls
    '''
    for sample in sample_list:
        msg = (" INFO: Calling single-exon CNVs on sample {}").format(sample.name)
        logging.info(msg)

        raw_single_cnv_name = ("{}.raw.single.exon.calls.bed").format(sample.name)
        raw_single_cnv_file = str(Path(sample.sample_folder) / raw_single_cnv_name)
        sample.add("raw_single_exon_calls", raw_single_cnv_file)
        o = open(raw_single_cnv_file, 'w')
        ratio_no_header = remove_bed_header(sample.ratio_file, "chr\tstart")
        a = pybedtools.BedTool(ratio_no_header)
        b = pybedtools.BedTool(sample.seg_calls_bed)
        c = a.intersect(b, wo=True, stream=True, v=True)
        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp  = line.split('\t')
            chr      = tmp[0]
            start    = tmp[1]
            end      = tmp[2]
            region   = tmp[3]
            n_region = 1
            log2_ratio = float(tmp[-1])
            fold_change = 2**(log2_ratio)
            cn = str(int(fold_change*2))
            tmp_list = [chr, start, end, region, str(n_region), str(log2_ratio)]
            if log2_ratio <= upper_del_threshold:
                cnvtype = "DEL"
                tmp_list.append(cn)
                tmp_list.append(cnvtype)
                o.write('\t'.join(tmp_list)+"\n")
            if log2_ratio >= dup_threshold:
                cnvtype = "DUP"
                tmp_list.append(cn)
                tmp_list.append(cnvtype)
                o.write('\t'.join(tmp_list)+"\n")
        o.close()
    return sample_list

def call_raw_segmented_cnvs(sample_list, upper_del_threshold, dup_threshold):
    '''
        Calling CNVs2
    '''
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)
    for sample in sample_list:
        msg = (" INFO: Calling segmented CNVs on sample {}").format(sample.name)
        logging.info(msg)

        seg_calls_name = ("{}.seg.calls.bed").format(sample.name)
        seg_calls_bed  = str(Path(sample.sample_folder)/seg_calls_name)
        sample.add("seg_calls_bed", seg_calls_bed)

        raw_seg_calls = ("{}.raw.seg.calls.bed").format(sample.name)
        raw_seg_calls_bed  = str(Path(sample.sample_folder)/raw_seg_calls)
        sample.add("raw_seg_calls", raw_seg_calls_bed)

        # We will use plain dicts and bedtools stuff, instead of pandas
        o = open(raw_seg_calls_bed, 'w')
        p = open(seg_calls_bed, 'w')

        o.write("chr\tstart\tend\tregions\tn_regions\tlog2_ratio\tcn\tcnvtype\n")
        with open (sample.segment_file) as seg:
            for line in seg:
                line = line.rstrip("\n")
                #chr11	19200000	26550000	133	0.5013	DEL
                if line.startswith('chr\tstart'):
                    continue
                tmp = line.split('\t')
                chr   = tmp[0]
                start = tmp[1]
                end   = tmp[2]
                regions = tmp[3]
                n_regions = tmp[4]
                log2_ratio= tmp[5]
                cn = int(tmp[6])
                cnvtype = ""
                if cn != 2:
                    if cn > 2:
                        if float(log2_ratio) >= dup_threshold:
                            cnvtype = "DUP"
                    else:
                        if float(log2_ratio) <= upper_del_threshold:
                            cnvtype = "DEL"
                    if cnvtype != "":
                        outline = ("{}\t{}\n").format(line, cnvtype)
                        o.write(outline)
                        p.write(outline)
        seg.close()
        o.close()
        p.close()
    return sample_list

def unify_raw_calls(sample_list):

    for sample in sample_list:
        df_list = []
        bed_list = [sample.raw_seg_calls, sample.raw_single_exon_calls]

        raw_calls_name = ("{}.raw.calls.bed").format(sample.name)
        raw_calls_bed = str(Path(sample.sample_folder)/raw_calls_name)
        sample.add("raw_calls_bed", raw_calls_bed)
        for bed in bed_list:
            if os.path.getsize(bed):
                df = pd.read_csv(bed, sep="\t")
                df_list.append(df)

        o = open(raw_calls_bed, 'w')
        for file in bed_list:
            with open(file) as f:
                for line in f:
                    o.write(line)
        o.close()
        # p.close()
    return sample_list

def call_cnvs(sample_list, upper_del_threshold, dup_threshold, z_score):
    '''
        Calling CNVs
    '''
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)
    for sample in sample_list:
        msg = (" INFO: Calling CNVs on sample {}").format(sample.name)
        logging.info(msg)

        cnv_calls_name = ("{}.calls.bed").format(sample.name)
        cnv_calls_bed  = str(Path(sample.sample_folder)/cnv_calls_name)
        sample.add("cnv_calls_bed", cnv_calls_bed)

        ratio_no_header = remove_bed_header(sample.ratio_file, "chr\tstart")
        seg_no_header   = remove_bed_header(sample.segment_file, "chr\tstart")

        # We will use plain dicts and bedtools stuff, instead of pandas
        tmp_calls = cnv_calls_bed.replace(".bed", ".tmp.bed")
        o = open(tmp_calls, 'w')

        with open (sample.segment_file) as seg:
            for line in seg:
                line = line.rstrip("\n")
                #chr11	19200000	26550000	133	0.5013	DEL
                if line.startswith('chr\tstart'):
                    continue
                tmp = line.split('\t')
                chr   = tmp[0]
                start = tmp[1]
                end   = tmp[2]
                n_bins= tmp[3]
                ratio = tmp[4]
                cnvtype = ""
                if float(ratio) >= dup_threshold:
                    cnvtype = "DUP"
                elif float(ratio) <= upper_del_threshold:
                    cnvtype = "DEL"
                if cnvtype != "":
                    outline = ("{}\t{}\n").format(line, cnvtype)
                    o.write(outline)
        seg.close()
        o.close()

        segmented_cnvs_dict = get_segmented_cnvs(ratio_no_header, tmp_calls)

        o = open(cnv_calls_bed, 'w')
        o.write("chr\tstart\tend\tregions\tratio\tcnvtype\tstd\n")

        for variant in calls_dict:
            arr = np.array(segmented_cnvs_dict[variant])
            std = round(np.std(arr), 3)
            outline = ("{}\t{}\n").format(variant, str(std))
            o.write(outline)
        o.close()
        sample.add("ready_cnv_bed", cnv_calls_bed)

        os.remove(tmp_calls)
        os.remove(ratio_no_header)
        os.remove(seg_no_header)

    return sample_list

def get_segmented_cnvs(ratio_no_header, tmp_calls):
    '''
    '''
    calls_dict = defaultdict(dict)
    a = pybedtools.BedTool(ratio_no_header)
    b = pybedtools.BedTool(tmp_calls)
    c = a.intersect(b, wo=True, stream=True)
    for line in iter(c):
        line = str(line)
        line = line.rstrip()
        tmp  = line.split('\t')
        variant = "\t".join(tmp[6:12])
        if not variant in calls_dict:
            calls_dict[variant] = list()
        else:
            calls_dict[variant].append(float(tmp[5]))
    return calls_dict
