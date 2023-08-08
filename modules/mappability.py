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
from multiprocessing import Pool



def split_bedfile(filename, num_chunks):
    """
    This function reads a bed file and splits it into smaller bed files.
    """
    with open(filename) as f:
        lines = f.readlines()
    
    chunk_size = len(lines) // num_chunks
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

    bed_chunks = []
    for i, chunk in enumerate(chunks):
        chunk_filename = f"{filename}_chunk_{i}.bed"
        with open(chunk_filename, 'w') as f:
            f.writelines(chunk)
        bed_chunks.append(pybedtools.BedTool(chunk_filename))
    
    return bed_chunks

def intersect_beds(args):
    a, b, output_file = args
    c = a.intersect(b, wo=True, stream=True)
    with open(output_file, "w") as f:
        for line in c:
            f.write(str(line))


def process_results(c):
    map_dict = defaultdict(dict)

    for line in c:
        tmp = line.split("\t")
        mappability = float(tmp[3])
        size = int(tmp[6]) - int(tmp[5])
        bases = int(tmp[-1])
        marginal = ((mappability * bases) / size) * 100
        coordinate = "\t".join(tmp[4:9])
        if not coordinate in map_dict:
            map_dict[coordinate] = marginal
        else:
            map_dict[coordinate] += marginal

    return map_dict

def annotate_mappability_bed(input_bed, mappability_bed, num_processes=3):
    """ """
    output_bed = input_bed.replace(".bed", ".map.bed")
    if not os.path.isfile(output_bed):

        msg = " INFO: Extracting off-target mappability"
        logging.info(msg)

        # Split input bed into chunks
        a = pybedtools.BedTool(mappability_bed)
        b_chunks = split_bedfile(input_bed, num_processes)  # Use the new split function

        temp_files = [f'temp_{i}.bed' for i in range(num_processes)]
        # Use multiprocessing pool to process bed intersections in parallel
        with Pool(processes=num_processes) as pool:  # Assuming num_processes here
            pool.map(intersect_beds, [(a, b_chunk, temp_file) for b_chunk, temp_file in zip(b_chunks, temp_files)])

        # Combine results
        combined_results = {}
        for temp_file in temp_files:
            with open(temp_file, 'r') as f:
                result = [line.strip() for line in f]
                combined_results.update(process_results(result))
            os.remove(temp_file)

        with open(output_bed, "w") as o:
            for region in combined_results:
                o.write(region + "\t" + str(combined_results[region]) + "\n")

    return output_bed


# def annotate_mappability_bed(input_bed, mappability_bed):
#     """ """
#     output_bed = input_bed.replace(".bed", ".map.bed")
#     if not os.path.isfile(output_bed):

#         msg = " INFO: Extracting off-target mappability"
#         logging.info(msg)

#         a = pybedtools.BedTool(mappability_bed)
#         b = pybedtools.BedTool(input_bed)
#         c = a.intersect(b, wo=True, stream=True)
#         map_dict = defaultdict(dict)

#         for line in iter(c):
#             line = str(line)
#             line = line.rstrip()
#             tmp = line.split("\t")
#             mappability = float(tmp[3])
#             size = int(tmp[6]) - int(tmp[5])
#             bases = int(tmp[-1])
#             marginal = ((mappability * bases) / size) * 100
#             coordinate = "\t".join(tmp[4:9])
#             if not coordinate in map_dict:
#                 map_dict[coordinate] = marginal
#             else:
#                 map_dict[coordinate] += marginal

#         o = open(output_bed, "w")
#         for region in map_dict:
#             o.write(region + "\t" + str(map_dict[region]) + "\n")
#         o.close()
#     return output_bed


def annotate_mappability(analysis_dict, ann_dict):
    """ """

    msg = " INFO: Extracting mappability"
    logging.info(msg)

    gc_bed = analysis_dict["gc_bed"]

    map_bed_name = os.path.basename(gc_bed).replace(".gc.bed", ".map.bed")
    map_bed = str(Path(analysis_dict["output_dir"]) / map_bed_name)

    if not os.path.isfile(map_bed):
        a = pybedtools.BedTool(ann_dict["mappability"])
        b = pybedtools.BedTool(gc_bed)
        c = a.intersect(b, wo=True, stream=True)
        map_dict = defaultdict(dict)

        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp = line.split("\t")
            mappability = float(tmp[3])
            size = int(tmp[6]) - int(tmp[5])
            bases = int(tmp[-1])
            marginal = ((mappability * bases) / size) * 100
            coordinate = "\t".join(tmp[4:9])
            if not coordinate in map_dict:
                map_dict[coordinate] = marginal
            else:
                map_dict[coordinate] += marginal

        o = open(map_bed, "w")
        for region in map_dict:
            o.write(region + "\t" + str(map_dict[region]) + "\n")
        o.close()
    else:
        msg = " INFO: Skipping mappability extraction"
        logging.info(msg)

    analysis_dict["map_bed"] = map_bed
    analysis_dict["ready_bed"] = map_bed

    return analysis_dict
