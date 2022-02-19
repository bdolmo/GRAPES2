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


def launch_normalization(sample_list, analysis_dict, ann_dict):
    '''
    '''
    for sample in sample_list:

        # Load a dataframe of raw coverage data
        df = pd.read_csv(analysis_dict['unified_raw_depth'], sep="\t")


        # Calculate median coverage
        median_depth = round(df[sample.name].median(), 2)
        sample.add('median_depth', median_depth)

    # Normalize by GC-content
    norm_factors = ['gc']
    df = normalize(analysis_dict, sample_list,  norm_factors)

    # Export normalized coverage
    normalized_depth_name = ("{}.normalized.depth.bed").format(analysis_dict['output_name'])
    normalized_depth = str(Path(analysis_dict['output_dir']) / normalized_depth_name)
    df.to_csv(normalized_depth, sep="\t", mode='w', index=None)

    analysis_dict['normalized_depth'] = normalized_depth

    return sample_list, analysis_dict


def norm_by_lib(row, sample_name, total_reads):
    '''
    '''
    length = row['end']-row['start']
    norm_by_len = round(10e4*(row[sample_name]/length),3)
    norm_lib = round(10e2*(norm_by_len/total_reads),3)
    return norm_lib

def normalize(analysis_dict, sample_list, fields):
    '''
    '''

    df = pd.read_csv(analysis_dict['unified_raw_depth'], sep="\t")

    fields_str = ','.join(fields)
    msg = (" INFO: Normalizing {}").format(fields_str)
    logging.info(msg)

    for sample in sample_list:
        sample_lib_tag = ("{}_normalized_library").format(sample.name)
        # df[sample_lib_tag] = round(1000000*(df[sample.name]/sample.ontarget_reads), 4)
        df[sample_lib_tag] = df.apply(norm_by_lib, sample_name=sample.name,
            total_reads=sample.ontarget_reads,  axis=1)

        cov_target = sample_lib_tag

        idx = 0
        for field in fields:
            idx+=1
            field_int = ("{}_integer").format(field)

            # Get integer field (gc or map) value
            df[field_int] = df[field].apply(int)

            median_field_cov_autosomes = ("median_{}_cov_autosomes").format(field)
            median_field_cov_chrX = ("median_{}_cov_chrX").format(field)
            median_field_cov_chrY = ("median_{}_cov_chrY").format(field)

            median_cov_autosomes = round(df[(df['chr'] != "chrX") & (df['chr'] != "chrY")][cov_target].median(),3)
            median_cov_chrX = round(df[(df['chr'] != "chrX")][cov_target].median(),3)
            median_cov_chrY = round(df[(df['chr'] == "chrY")][cov_target].median(),3)

            stats_dict = {
                'median_cov_autosomes': median_cov_autosomes,
                'median_cov_chrX': median_cov_chrX,
                'median_cov_chrY': median_cov_chrY
            }

            # Group coverage by field and calculate the median
            df[median_field_cov_autosomes]  = df[(df['chr'] != "chrX") & (df['chr'] != "chrY")]\
                .groupby(field_int)[cov_target].transform(median_by_interval,
                cov_target=cov_target, median=median_cov_autosomes)
            df[median_field_cov_chrX]  = df[(df['chr'] == "chrX")]\
                .groupby(field_int)[cov_target].transform(median_by_interval,
                 cov_target=cov_target, median=median_cov_chrX)
            df[median_field_cov_chrY]  = df[(df['chr'] == "chrY")]\
                .groupby(field_int)[cov_target].transform(median_by_interval,
                cov_target=cov_target, median=median_cov_chrY)

            # Apply rolling median
            normalized_field = ("{}_normalized_{}").format(sample.name, field)
            df[normalized_field] = df.apply(apply_normalization, cov_target=cov_target,
                stats_dict=stats_dict, field=field, axis=1)
            if idx == len(fields):
                normalized_final = ("{}_normalized_final").format(sample.name)
                df[normalized_final] = df[normalized_field]
    return df

def median_by_interval(x, cov_target, median):
    '''
    '''

    median_interval = round(x.median(),3)
    if not median:
        median_interval = median
    return median_interval

def apply_normalization(row, cov_target, stats_dict, field):
    '''
        Normalization
    '''
    median_field_chrX = ("median_{}_cov_chrX").format(field)
    median_field_chrY = ("median_{}_cov_chrY").format(field)
    median_field_autosomes = ("median_{}_cov_autosomes").format(field)

    if row['chr'] == "chrX":
        norm_cov = round((row[cov_target]*stats_dict['median_cov_chrX'])/row[median_field_chrX],2)
    elif row['chr'] == "chrY":
        norm_cov = round((row[cov_target]*stats_dict['median_cov_chrY'])/row[median_field_chrY],2)
    else:
        norm_cov = round((row[cov_target]*stats_dict['median_cov_autosomes'])/row[median_field_autosomes],2)

    return norm_cov
