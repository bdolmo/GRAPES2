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

import networkx as nx
from collections import defaultdict
from itertools import combinations



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


def insert_size(reads):
    pos = sorted([read.pos for read in reads])
    return pos[-1] - pos[0]


def create_breakreads_graph(breakreads, mean_insert_size, sd_insert_size):
    G = nx.Graph()
    max_insert_size = mean_insert_size + 10 * sd_insert_size

    for qname, reads in breakreads.items():
        if len(reads) != 2:
            continue

        G.add_node(qname, reads=reads)

    for qname, other_qname in combinations(G.nodes, 2):
        reads = G.nodes[qname]['reads']
        other_reads = G.nodes[other_qname]['reads']

        if is_overlapping(reads, other_reads, max_insert_size):
            G.add_edge(qname, other_qname)

    return G


def is_overlapping(reads_a, reads_b, max_insert_size):
    chromosomes_a = {read.reference_id for read in reads_a}
    chromosomes_b = {read.reference_id for read in reads_b}

    if chromosomes_a != chromosomes_b:
        return False

    pos_a = sorted([read.pos for read in reads_a])
    pos_b = sorted([read.pos for read in reads_b])

    if abs(pos_a[0] - pos_b[0]) > max_insert_size or abs(pos_a[1] - pos_b[1]) > max_insert_size:
        return False


    insert_sizes = [insert_size(pair) for pair in combinations(reads_a + reads_b, 2) if pair[0].reference_id == pair[1].reference_id]

    if not insert_sizes:
        return False

    mean_cluster_insert_size = np.mean(insert_sizes)
    sd_insert_size = np.std(insert_sizes)
    within_limit = [size for size in insert_sizes if size <= mean_cluster_insert_size + 3 * sd_insert_size]

    return len(within_limit) == len(insert_sizes)


def cluster_breakreads(G):
    cliques = nx.algorithms.community.k_clique_communities(G, 2)
    clusters = [G.subgraph(clique) for clique in cliques]
    return clusters


def scan_breakreads(bam: str, bed: str, fasta: str):
    """ """

    regions = []
    with open(bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
            regions.append(rec)
    f.close()

    breakreads = defaultdict(list)
    bam_file = pysam.AlignmentFile(bam, "rb")
    for read in bam_file:
        if not read.cigarstring:
            continue
        if not read.is_proper_pair:
            if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
                # breakreads.append(read.query_sequence)
                breakreads[read.qname].append(read)

            if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
                # breakreads.append(read.query_sequence)
                breakreads[read.qname].append(read)

    return breakreads


    # sys.exit()

    # bam_file = pysam.AlignmentFile(bam, "rb")
    # for region in regions:
    #     if not region.start == 179424037:
    #         continue

    #     breakreads = []
    #     for read in bam_file.fetch(region.chr, region.start, region.end):
    #         if not read.cigarstring:
    #             continue
    #         if not read.is_proper_pair:
    #             if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
    #                 breakreads.append(read.query_sequence)
    #             if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
    #                 breakreads.append(read.query_sequence)

    #     oas = OverlapAssembler(breakreads, 21)
    #     contigs = oas.compute_overlaps()
    #     ref = pysam.FastaFile(fasta)
    #     reference = ref.fetch('chr2', 179431346, 179437577)
    #     seqs = contigs
    #     # print(seqs)
    #     # sys.exit()

    #     # for seq in seqs:
    #     #     print(str(seq))


    #     blat = Blat(reference, chr="chr2", start=179431346, end=179437577)
    #     vars = blat.align(seqs)
    #     print(vars)

    # pass


def call_structural_variants(bam: str, bed: str, fasta: str):
    """ """

    bam_stats = get_bam_stats(bam, bed)

    breakreads = scan_breakreads(bam, bed, fasta)
    breakreads_graph = create_breakreads_graph(breakreads, bam_stats['isize_median'], bam_stats['isize_mad'])
    clusters = cluster_breakreads(breakreads_graph)

    print('Total clusters:', len(clusters))
    seqs = []
    for i, cluster in enumerate(clusters, start=1):
        print(f'Cluster {i}:')
        for node, data in cluster.nodes(data=True):
            reads = data['reads']
            print(f'  {node}:')
            for read in reads:
                print(f'    {read.reference_name}:{read.pos}-{read.pos + read.qlen}')
                seqs.append(read.seq)
  

    oas = OverlapAssembler(seqs, 21)
    contigs = oas.compute_overlaps()
    for s in contigs:
        print(s)
    ref = pysam.FastaFile(fasta)
    reference = ref.fetch('chr2', 179431346, 179437577)

    blat = Blat(k=21, ref=reference, chr="chr2", start=179431346, end=179437577)

    # blat = Blat(reference, chr="chr2", start=179431346, end=179437577)
    for i in range(0, 1):
        for seq in contigs:
            print("aligning", seq)
            vars = blat.align(seqs=[seq])
            print(vars)


if __name__ == "__main__":
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221118-CGC64001/SUDD_147/RB34017/BAM_FOLDER/RB34017.rmdup.bam"
    bam = "/home/bdelolmo/RB33925.rmdup.bam"
    bed = "/home/bdelolmo/BED/gendiag_85.CDS.bed"
    fasta = "/home/bdelolmo/REF_DIR/hg19/ucsc.hg19.fasta"
    call_structural_variants(bam, bed, fasta)
