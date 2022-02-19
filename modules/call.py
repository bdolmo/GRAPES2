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

def export_all_calls(sample_list,analysis_dict):
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

def call_cnvs_2(sample_list, upper_del_threshold, dup_threshold):
    '''
        Calling CNVs debug
    '''
    upper_del_threshold = float(upper_del_threshold)
    dup_threshold = float(dup_threshold)
    for sample in sample_list:
        msg = (" INFO: Calling CNVs on sample {}").format(sample.name)
        logging.info(msg)

        cnv_calls_name = ("{}.calls.bed").format(sample.name)
        cnv_calls_bed  = str(Path(sample.sample_folder)/cnv_calls_name)
        sample.add("cnv_calls_bed", cnv_calls_bed)

        # We will use plain dicts and bedtools stuff, instead of pandas
        o = open(cnv_calls_bed, 'w')
        o.write("chr\tstart\tend\tregions\tlog2_ratio\tcn\tcnvtype\n")
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
                log2_ratio= tmp[4]
                cn = int(tmp[5])
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
        seg.close()
        o.close()
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
