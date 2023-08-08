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


def split_file(file_name, lines_per_chunk):
    """
    Split a file into chunks, each containing a certain number of lines.
    Returns a list of the chunk file names.
    """
    chunk_files = []
    with open(file_name, 'r') as f:
        chunk = []
        i = 0
        for line in f:
            chunk.append(line)
            if len(chunk) >= lines_per_chunk:
                chunk_file_name = f"{file_name}_chunk_{i}.bed"
                with open(chunk_file_name, 'w') as chunk_file:
                    chunk_file.write(''.join(chunk))
                chunk_files.append(chunk_file_name)
                chunk = []
                i += 1
        if chunk:  # If there are leftover lines, write them to a new chunk file
            chunk_file_name = f"{file_name}_chunk_{i}.bed"
            with open(chunk_file_name, 'w') as chunk_file:
                chunk_file.write(''.join(chunk))
            chunk_files.append(chunk_file_name)
    return chunk_files

def annotate_gc_bed(input_bed, output_dir, genome_fasta):
    """
    Compute GC content for a large BED file in chunks.
    Writes the results for each chunk to a separate file in the output directory.
    """
    # Split the input BED file into chunks

    output_bed = input_bed.replace(".bed", ".gc.bed")
    if os.path.isfile(output_bed):
        return output_bed

    lines_per_chunk = 100
    chunk_files = split_file(input_bed, lines_per_chunk)

    chunks_with_gc_list = []
    # Compute GC content for each chunk
    for chunk_file in chunk_files:
        a = pybedtools.BedTool(chunk_file)
        b = a.nucleotide_content(fi=genome_fasta, pattern="CG", C=True)

        # Write the results to a file in the output directory
        output_file = os.path.join(output_dir, os.path.basename(chunk_file).replace(".bed", ".gc.bed"))
        b.saveas(output_file)

        chunks_with_gc_list.append(output_file)

    df_list = []
    for file in chunks_with_gc_list:
        df = pd.read_csv(file, sep="\t", header=None)
        
        df = df.iloc[:, [0,1,2,3,5]]
        df_list.append(df)

    master_df = pd.concat(df_list, ignore_index=True)
    master_df.to_csv(output_bed, sep="\t", header=False, index=False)
    
    tmp_bed = output_bed.replace(".bed", ".tmp.bed")
    o = open(tmp_bed, "w")
    with open(output_bed) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\t")
            tmp = line.split("\t")
            tmp[4] = str(round(float(tmp[4])*100,2))
            out_str = '\t'.join(tmp)
            o.write(out_str+"\n")
    o.close()
    os.remove(output_bed)
    os.rename(tmp_bed, output_bed)

    # Delete the chunk files to clean up
    for chunk_file in chunk_files:
        os.remove(chunk_file)

    for chunk_file in chunks_with_gc_list:
        os.remove(chunk_file)

    return output_bed 



# def annotate_gc_bed(input_bed, genome_fasta):

#     msg = " INFO: Extracting gc content"
#     logging.info(msg)
#     print(msg)


#     output_bed = os.path.basename(input_bed).replace(".bed", ".gc.bed")
#     a = pybedtools.BedTool(input_bed)
#     b = a.nucleotide_content(fi=genome_fasta, pattern="CG", C=True)


#     o = open(output_bed, "w")
#     for line in iter(b):
#         line = str(line)
#         line = line.rstrip()
#         tmp = line.split("\t")
#         tmp[5] = str(round(100 * float(tmp[5]), 2))
#         indices = [0, 1, 2, 3, 5]
#         new_line = "\t".join([tmp[i] for i in indices])
#         o.write(new_line + "\n")
#     o.close()
#     print("done")

#     return output_bed

def annotate_gc(analysis_dict):
    """
    Add gc content
    """
    bed = analysis_dict["bed"]

    gc_bed_name = os.path.basename(bed).replace(".bed", ".gc.bed")
    gc_bed = str(Path(analysis_dict["output_dir"]) / gc_bed_name)

    msg = " INFO: Extracting gc content for roi file {}".format(bed)
    a = pybedtools.BedTool(bed)
    b = a.nucleotide_content(fi=analysis_dict["reference"], pattern="CG", C=True)

    o = open(gc_bed, "w")
    for line in iter(b):
        line = str(line)
        line = line.rstrip()
        tmp = line.split("\t")
        tmp[5] = str(round(100 * float(tmp[5]), 2))
        indices = [0, 1, 2, 3, 5]
        new_line = "\t".join([tmp[i] for i in indices])
        o.write(new_line + "\n")
    o.close()

    analysis_dict["gc_bed"] = gc_bed
    return analysis_dict
