import os
import sys
import logging
import glob
import subprocess
import numpy as np
import pandas as pd


def remove_tmp_files(input_dir):
    """ """
    count_files = glob.glob(input_dir + "/*.bam_counts.bed")
    for f in count_files:
        os.remove(f)

    isize_files = glob.glob(input_dir + "/*.bam_isizes.bed")
    for f in isize_files:
        os.remove(f)

    coverage_files = glob.glob(input_dir + "/*.bam_coverage.bed")
    for f in coverage_files:
        os.remove(f)

    per_base_files = glob.glob(input_dir + "/*.per.base.coverage.bed")
    for f in per_base_files:
        os.remove(f) 

    normalized_base = glob.glob(input_dir + "/*.normalized.per.base.bed")
    for f in normalized_base:
        os.remove(f) 

    tmp_calls = glob.glob(input_dir + "/*.tmp.rawcalls.bed") 
    for f in tmp_calls:
        os.remove(f) 



def sort_bed_file(input_bed):
    """ """

    output_bed = input_bed.replace(".bed", ".tmp.bed")

    # Load the BED file into a DataFrame
    df = pd.read_csv(input_bed, sep='\t', names=['chr', 'start', 'end', 'name'], header=None)

    if df.empty:
        return df

    # Replace 'chrX' and 'chrY' with temporary placeholders
    df['chr'] = df['chr'].replace({'chrM': 'chr0', 'chrX': 'chr23', 'chrY': 'chr24'})

    # Remove 'chr' prefix for sorting but keep it in a separate column
    df['chr_num'] = df['chr'].str.replace('chr', '').astype(int)

    # Sort by chromosomal position
    df = df.sort_values(['chr_num', 'start', 'end'])

    # Drop the temporary column
    df = df.drop(columns=['chr_num'])

    # Replace temporary placeholders with 'chrX' and 'chrY'
    df['chr'] = df['chr'].replace({'chr0':'chrM', 'chr23': 'chrX', 'chr24': 'chrY'})

    # Save the sorted DataFrame back to a BED file
    df.to_csv(output_bed, sep='\t', header=False, index=False)

    os.remove(input_bed)
    os.rename(output_bed, input_bed)


def signal_to_noise(data):
    """ """
    median = np.median(data)
    std = np.std(data)
    if std == 0:
        s2n = 0
    else:
        s2n = abs(median / std)
    return round(s2n, 3)


def remove_bed_header(file, pattern):
    """ """
    no_header_file = file.replace(".bed", ".noheader.bed")
    nh = open(no_header_file, "w")
    with open(file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(pattern):
                continue
            nh.write(line + "\n")
    f.close()
    nh.close()
    return no_header_file


def assign_genotype_based_on_cn(cn):
    """
    Simple genotype assignment based on copy number
    """
    if cn == 0:
        gt = "1/1"
    elif cn == 1:
        gt = "0/1"
    else:
        gt = "./1"
    return gt


def validate_bed(bed):
    """ """
    pass


class MissingInputBamFiles(Exception):
    pass


class NgsUtilsError(Exception):
    pass


def check_executable(program, dump_messages=True):
    """ """
    bashCommand = ("which {}").format(program)
    p1 = subprocess.run(
        bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output = p1.stdout.decode("UTF-8")
    error = p1.stderr.decode("UTF-8")
    executable = False
    if not error:
        if output:
            msg = (" INFO: Found executable for {}").format(os.path.basename(program))
            logging.info(msg)
            pass
        else:
            msg = (" ERROR: Unable to execute {}").format(os.path.basename(program))
            raise NgsUtilsError(msg)
            # print(msg)
    else:
        msg = (" ERROR: Unable to execute {}").format(os.path.basename(program))
        raise NgsUtilsError(msg)


def get_bam_files(input_dir):
    """
    Get bam files from input dir
    """
    bam_list = glob.glob(input_dir + "/*.bam")

    if not bam_list:
        msg = (" ERROR: missing input bam files from {} directory").format(input_dir)
        # logging.error(msg)
        raise MissingInputBamFiles(msg)

    return bam_list
