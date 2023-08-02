import sys
import os
import pysam
import re
import numpy as np
from statistics import median
from Bio import Align
from modules.bed import BedRecord, load_bed_file
from collections import defaultdict
import subprocess

def call_structural_variants(bam, bed, fasta, output_dir, sample_name, ngs_utils_dict, ann_dict):
    """ """

    min_size = 15
    max_size = 1000000
    threads = 1

    command = [ngs_utils_dict["grapes_sv"], 
            '-b', bam,
            '-g', fasta,
            '-n', sample_name,
            '-o', output_dir,
            '--wes', "on",
            '--wgs', "off",
            '-c', str(5),
            '-s', str(11),
            '-r', str(5),
            '-l', str(10),
            '-a', str(25),
            '-m', str(4),
            '-i', str(min_size),
            '-j', str(max_size),
            '--find-small', "on",
            '--find-large', "on",
            '-t', str(threads),
            '-e', ann_dict["blacklist"]]
    
    result = subprocess.run(command, capture_output=True, text=True)
    
    tmp_files = {
        "FR": os.path.join(output_dir, f"{sample_name}.FR.bam"),
        "RF": os.path.join(output_dir, f"{sample_name}.RF.bam"),
        "FF": os.path.join(output_dir, f"{sample_name}.FF.bam"),
        "RR": os.path.join(output_dir, f"{sample_name}.RR.bam"),
        "SR": os.path.join(output_dir, f"{sample_name}.SR.bam"),
        "info": os.path.join(output_dir, f"{sample_name}.discordantInfo.txt"),
        "fastq": os.path.join(output_dir, f"{sample_name}.fastq"),
        "FRcluster": os.path.join(output_dir, f"{sample_name}.FR.clusters.bed"),
        "RFcluster": os.path.join(output_dir, f"{sample_name}.RF.clusters.bed"),
        "FFcluster": os.path.join(output_dir, f"{sample_name}.FF.clusters.bed"),
        "RRcluster": os.path.join(output_dir, f"{sample_name}.RR.clusters.bed"),
    }


    if result.returncode != 0:
        print(f" INFO: Error running GRAPES_SV:\n{result.stderr}")
    else:
        print(f" INFO: GRAPES_SV ran successfully")

        raw_file_name = f"{sample_name}.tmp.rawcalls.bed"
        raw_file = os.path.join(output_dir, raw_file_name)

        bed_out = os.path.join(output_dir, f"{sample_name}.GRAPES2.breakpoints.bed")

        if not os.path.isfile(bed_out):
            o = open(bed_out, "w")
        else:
            o = open(bed_out, "w")

        seen_coord = {}
        with open(raw_file) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                coordinates = '\t'.join(tmp[0:2])
                if not coordinates in seen_coord:
                    seen_coord[coordinates] = True
                else:
                    continue
                o.write(line+"\n")
        o.close()
        f.close()

    for file_type in tmp_files:
        if os.path.isfile(tmp_files[file_type]):
            os.remove(tmp_files[file_type])


if __name__ == "__main__":
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221118-CGC64001/SUDD_147/RB34017/BAM_FOLDER/RB34017.rmdup.bam"
    bam = "/home/bdelolmo/RB33925.rmdup.bam"
    bed = "/home/bdelolmo/BED/gendiag_85.CDS.bed"
    fasta = "/home/bdelolmo/REF_DIR/hg19/ucsc.hg19.fasta"

    # blat = Blat(reference, chr="chr2", start=179431346, end=179437577)

    # call_structural_variants(bam, bed, fasta)
