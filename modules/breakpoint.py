import sys
import os
import pysam
import re
import numpy as np
from statistics import median
from bed import BedRecord
from blat import Blat
from Bio import Align
from assembler import DeBruijnAssembler, OverlapAssembler
# import edlib
# from swalign import swalign


def get_bam_stats(bam: str, bed: str, N=5000):
    regions = []
    with open(bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
            regions.append(rec)
    f.close()

    insert_sizes = []
    count = 0

    bam_file = pysam.AlignmentFile(bam, "rb")
    for region in regions:
        for read in bam_file.fetch(region.chr, region.start, region.end):
            insert_sizes.append(abs(read.template_length))
            if count == N:
                break
            count += 1
    isize_median = median(insert_sizes)
    isize_mad = median([abs(number - isize_median) for number in insert_sizes])
    upper_limit = isize_median + (10 * isize_mad)

    bam_stats = {
        "isize_median": isize_median,
        "isize_mad": isize_mad,
        "isize_threshold": upper_limit,
    }
    return bam_stats


def scan_breakreads(bam: str, bed: str, fasta: str):
    """ """

    bam_stats = get_bam_stats(bam, bed)

    regions = []
    with open(bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
            regions.append(rec)
    f.close()

    bam_file = pysam.AlignmentFile(bam, "rb")
    for region in regions:
        if not region.start == 179424037:
            continue

        breakreads = []
        for read in bam_file.fetch(region.chr, region.start, region.end):
            if not read.cigarstring:
                continue
            if not read.is_proper_pair:
                if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
                    breakreads.append(read.query_sequence)
                if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
                    breakreads.append(read.query_sequence)

        oas = OverlapAssembler(breakreads, 21)
        contigs = oas.compute_overlaps()

        ref = pysam.FastaFile(fasta)
        reference = ref.fetch('chr2', 179431346, 179437577)
        seqs = contigs
        # print(seqs)
        # sys.exit()
        blat = Blat(reference, chr="chr2", start=179431346, end=179437577)
        blat.align(seqs)

    pass


if __name__ == "__main__":
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221118-CGC64001/SUDD_147/RB34017/BAM_FOLDER/RB34017.rmdup.bam"
    bam = "/home/bdelolmo/RB33925.rmdup.bam"
    bed = "/home/bdelolmo/BED/gendiag_85.CDS.bed"
    fasta = "/home/bdelolmo/REF_DIR/hg19/ucsc.hg19.fasta"
    scan_breakreads(bam, bed, fasta)
