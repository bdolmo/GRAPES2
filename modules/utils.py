import os
import sys
import logging
import glob
import subprocess
import numpy as np


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
