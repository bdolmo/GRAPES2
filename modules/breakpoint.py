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

from ssw_aligner import local_pairwise_align_ssw




def test_blat():
    """ """

    ref = "AGGGAAAGCATTACGACTCACGCGGGGCAATCTCAATGACGACGCGCAGAAAAAATTTTTACGACGCTCGGGGAAAAACTGGCATTTCAGCACTAGCGGCGATATACTG"
    seqs = ["TACGACcCACGCGGGGGACGCGCAGAAAAAATTTTTACGACGCTCGGGGAAAAACTGGCAT"]
    alignment = local_pairwise_align_ssw(seqs[0],
                                         ref,
                                         gap_open_penalty=11,
                                         gap_extend_penalty=1,
                                         match_score=2,
                                         mismatch_score=-3)
    print(alignment)
    # seqs = ["GACTCACGCGGGGCAATCCGTCGTAAAAATTTTTTCTGCCGTCGTCATTGACTCGGGGAAAAACTGGCAT"]

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -5
    aligner.open_gap_score = -18
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    # 2,-8,-12,0

    alignments = aligner.align(ref, seqs[0])

    # print(alignments[0])
    # print(alignments[0].aligned)
    # sys.exit()
    blat = Blat(ref)
    blat.align(seqs)
#
# def test_swalign():
#     """ """
#     ref = "TTTAAACCGACGACGTACGACGAGCAGCTAGCACAGGAGGATGAGCTGATCCGCCCCCCAGTCAAGGGGATATTATAACGAGCTAGCGGATATTATAACGAGCTAGCTGAAAACCGAGC"
#     alt = "TTTAAACCGACGACGTACGACGAGCAGCTAGCACAGGAGGATGAGCTGATCCGCCCCCCAGTCAAGGGGATATTATAACGAGCTAGCGGATATTATAACGAGCTAGCTGAAAACCGAGC"
#     match = 2
#     mismatch = -1
#
#     scoring = swalign.NucleotideScoringMatrix(match, mismatch)
#     sw = swalign.LocalAlignment(
#         scoring, gap_penalty=-12, gap_extension_penalty=0, gap_extension_decay=4
#     )
#     aln = sw.align(ref, alt)
#     aln.dump()


def get_bam_stats(bam: str, bed: str, N=5000):
    """ """
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
                # if abs(read.template_length) > bam_stats["isize_threshold"]:
                #     print(read)
                if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
                    # print(read)
                    breakreads.append(read.query_sequence)
                if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
                    # print(read)
                    breakreads.append(read.query_sequence)
        # n = 0
        # for read in breakreads:
        #     print(f">{str(n)}seq")
        #     print(read)
        #     n+=1
        # for read in breakreads[35:45]:
        #     print(read)
        oas = OverlapAssembler(breakreads, 21)
        contigs = oas.compute_overlaps()

        ref = pysam.FastaFile(fasta)
        reference = ref.fetch('chr2', 179431346, 179437577)

        seqs = [contigs[-1]]
        # print(reference)
        blat = Blat(reference)
        blat.align(seqs)
        # print(breakreads[25:35])
        # dbg = DeBruijnAssembler(breakreads[35:45], 25)
        # contig_list = dbg.eulerian_walk()
        # print(contig_list)
    pass


if __name__ == "__main__":
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221118-CGC64001/SUDD_147/RB34017/BAM_FOLDER/RB34017.rmdup.bam"
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221103-CGC64001/SUDD_147/RB33925/BAM_FOLDER/RB33925.rmdup.bam"
    bed = "/home/bdelolmo/Escriptori/BED/gendiag_85.CDS.bed"
    fasta = "/home/bdelolmo/Escriptori/reference/ucsc.hg19.fasta"
    # test_swalign()
    test_blat()
    scan_breakreads(bam, bed, fasta)
