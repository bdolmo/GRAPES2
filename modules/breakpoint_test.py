import sys
import os
import pysam
import re
import numpy as np
from statistics import median
from Bio import Align
from modules.bed import BedRecord, load_bed_file
from modules.blat import Blat
from modules.assembler import DeBruijnAssembler, OverlapAssembler
from intervaltree import Interval, IntervalTree
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


# def create_breakreads_graph(breakreads, mean_insert_size, sd_insert_size):

#     msg = " INFO: Creating breakreads graph"
#     print(msg)

#     G = nx.Graph()
#     max_insert_size = mean_insert_size + 10 * sd_insert_size

#     for qname, reads in breakreads.items():
#         G.add_node(qname, reads=reads)
    
#     # sys.exit()

#     for qname, other_qname in combinations(G.nodes, 2):
#         reads = G.nodes[qname]['reads']
#         other_reads = G.nodes[other_qname]['reads']
#         if len(reads) > 1:
#             if is_overlapping_discordants(reads, other_reads, max_insert_size):
#                 G.add_edge(qname, other_qname)
#         else:
#             if is_overlapping_softclipped(reads, other_reads):
#                 G.add_edge(qname, other_qname)
#     return G



def create_breakreads_graph(breakreads, mean_insert_size, sd_insert_size):
    msg = " INFO: Creating breakreads graph"
    print(msg)

    G = nx.Graph()
    max_insert_size = mean_insert_size + 10 * sd_insert_size

    # Create an IntervalTree
    itree = IntervalTree()
    itrees = defaultdict(IntervalTree)
    # Add reads to the graph and the interval tree
    for qname, reads in breakreads.items():
        G.add_node(qname, reads=reads)

        # Compute the interval for this read and add it to the tree
        # Here, we consider an interval to span from the start of the first read to the end of the last read
        # start = min(read.pos for read in reads)
        # end = max(read.pos + read.alen for read in reads)
        for read in reads:
            start = read.pos
            end = read.pos + read.alen
            itrees[read.reference_id].add(Interval(start, end, qname))

    # Query the interval tree instead of looping over all pairs of nodes
    for qname in G.nodes:
        reads = G.nodes[qname]['reads']

        for read in reads:
            start = read.pos-50
            end = read.pos + read.alen
            tree = itrees[read.reference_id]

            # Use the interval tree to get all overlapping nodes
            for interval in tree[start:end]:
                other_qname = interval.data
                other_reads = G.nodes[other_qname]['reads']

                if len(reads) > 1:
                    if is_overlapping_discordants(reads, other_reads, max_insert_size):
                        G.add_edge(qname, other_qname)
                else:
                    if is_overlapping_softclipped(reads, other_reads):
                        G.add_edge(qname, other_qname)
    return G

def is_overlapping_softclipped(reads_a, reads_b):
    # adjust for soft-clipped portions

    chromosomes_a = {read.reference_id for read in reads_a}
    chromosomes_b = {read.reference_id for read in reads_b}

    if chromosomes_a != chromosomes_b:
        return False

    pos_a = [get_breakpoint(read) for read in reads_a]
    pos_b = [get_breakpoint(read) for read in reads_b]

    # loop over each position in reads_a and reads_b
    for pos_a in pos_a:
        for pos_b in pos_b:
            if abs(pos_a - pos_b) <= 200: # within 10bps
                return True
    return False


def is_overlapping_discordants(reads_a, reads_b, max_insert_size):
    chromosomes_a = {read.reference_id for read in reads_a}
    chromosomes_b = {read.reference_id for read in reads_b}

    if chromosomes_a != chromosomes_b:
        return False

    pos_a = sorted([read.pos for read in reads_a])
    pos_b = sorted([read.pos for read in reads_b])

    if len(pos_a) < 2 or len(pos_b) < 2:
        return False

    if abs(pos_a[0] - pos_b[0]) > max_insert_size or abs(pos_a[1] - pos_b[1]) > max_insert_size:
        return False

    insert_sizes = [insert_size(pair) for pair in combinations(reads_a + reads_b, 2) if pair[0].reference_id == pair[1].reference_id]

    if not insert_sizes:
        return False

    mean_cluster_insert_size = np.mean(insert_sizes)
    sd_insert_size = np.std(insert_sizes)
    within_limit = [size for size in insert_sizes if size <= mean_cluster_insert_size + 3 * sd_insert_size]

    return len(within_limit) == len(insert_sizes)


def get_breakpoint(read):
    # if softclipping at the start, breakpoint is the start of the aligned sequence
    if re.match(r'^[0-9]+S', read.cigarstring):
        return read.pos
    # if softclipping at the end, breakpoint is the end of the aligned sequence
    elif re.match(r'[0-9]+M([0-9]+I)?[0-9]+S$', read.cigarstring):
        # print(read.cigarstring, read.pos, read.pos + get_aligned_length(read.cigarstring))

        return read.pos + get_aligned_length(read.cigarstring)
    # if no softclipping, return the start of the aligned sequence
    else:
        return read.pos

def get_aligned_length(cigar):
    # extract all 'M' (match) operations from the CIGAR string and sum their lengths
    match_lengths = [int(x[:-1]) for x in re.findall(r'[0-9]+M', cigar)]
    return sum(match_lengths)


# def cluster_breakreads(G):
#     cliques = nx.algorithms.community.k_clique_communities(G, 20)
#     clusters = [G.subgraph(clique) for clique in cliques]
#     return clusters



# def cluster_breakreads_naive(bam: str, breakreads: dict, window=5):
#     clusters = defaultdict(list)
#     bam_file = pysam.AlignmentFile(bam, 'rb')
#     # For each breakread group
#     for qname, reads in breakreads.items():
#         for read in reads:
#             chromosome_name = bam_file.get_reference_name(read.reference_id)
#             if not chromosome_name in clusters:
#                 clusters[chromosome_name] = defaultdict(list)
#             # Consider both the start and end position of the read for clustering
#             # for pos in [read.reference_start, read.reference_end]:
            
#             pos = get_breakpoint(read)
#             if pos is None:
#                 continue
#             # Find a cluster that this position belongs to
#             found = False
#             for key in clusters[chromosome_name]:
#                 if abs(pos - key) <= window:
#                     clusters[chromosome_name][key].append(read)
#                     found = True
#                     break
#             # If no appropriate cluster is found, start a new one
#             if not found:
#                 clusters[chromosome_name][pos].append(read)

#     # Collect the keys of the clusters to delete
#     to_delete = [(chromosome, cluster) for chromosome in clusters
#                  for cluster in clusters[chromosome] if len(clusters[chromosome][cluster]) < 5]

#     # Delete the collected clusters
#     for chromosome, cluster in to_delete:
#         del clusters[chromosome][cluster]
    
#     for chromosome in clusters:
#         for pos in clusters[chromosome]:

#             print("#####", chromosome, pos, len(clusters[chromosome][pos]))
#             for read in clusters[chromosome][pos]:
#                 print(read)

#     bam_file.close()
#     return clusters


def overlaps(read, blacklisted_regions):
    for region in blacklisted_regions:
        if read.reference_start < region.end and read.reference_end > region.start:
            return True
    return False


def get_total_softclipped_bases(cigar_string):
    # Regex to match 'S' operations in the CIGAR string
    pattern = r'(\d+)S'
    
    # Find all 'S' operations in the CIGAR string
    matches = re.findall(pattern, cigar_string)
    
    # Convert lengths to integers and sum them to get the total soft-clipped bases
    total_softclipped_bases = sum(int(length) for length in matches)
    
    return total_softclipped_bases

def scan_breakreads(bam: str, bed: str, blacklist_bed: str,  fasta: str):
    """ """
    msg = " INFO: Scanning breakreads"
    print(msg)

    regions = []
    analysis_chromosomes = []
    with open(bed) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
            regions.append(rec)
            analysis_chromosomes.append(tmp[0])
    f.close()

    blacklisted_regions = load_bed_file(blacklist_bed)

    breakreads = defaultdict(lambda: defaultdict(list))  # dictionary of dictionary
    bam_file = pysam.AlignmentFile(bam, "rb")
    for read in bam_file:

        ref_name = bam_file.get_reference_name(read.reference_id)
        if ref_name not in analysis_chromosomes:
            continue
        if blacklisted_regions[ref_name][read.reference_start:read.reference_end]:
            continue
        if read.mapping_quality < 20:
            continue
        if not read.cigarstring:
            continue
        if get_total_softclipped_bases(read.cigarstring) < 10:
            continue

        if read.is_reverse:
            read.query_sequence = read.get_forward_sequence()
        if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
            breakreads[ref_name][read.reference_start].append(read)
        if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
            ops = re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)
            position = int(ops[0][0])+read.reference_start
            breakreads[ref_name][position].append(read)
    return breakreads


def cluster_breakreads(breakreads, distance=500):
    clustered_breakreads = defaultdict(lambda: defaultdict(list))

    for ref_name in breakreads:
        # Get a sorted list of positions
        positions = sorted(breakreads[ref_name].keys())
        cluster_start = positions[0]

        for position in positions[1:]:
            # If the difference between the current position and the start of the cluster is larger than the threshold
            if position - cluster_start > distance:
                # Start a new cluster
                cluster_start = position
            # Add the reads to the current cluster
            clustered_breakreads[ref_name][cluster_start].extend(breakreads[ref_name][position])

    return clustered_breakreads



# def scan_breakreads(bam: str, bed: str, blacklist_bed: str,  fasta: str):
#     """ """
#     msg = " INFO: Scanning breakreads"
#     print(msg)

#     regions = []
#     analysis_chromosomes = []
#     with open(bed) as f:
#         for line in f:
#             line = line.rstrip("\n")
#             tmp = line.split("\t")
#             rec = BedRecord(tmp[0], int(tmp[1]), int(tmp[2]))
#             regions.append(rec)
#             analysis_chromosomes.append(tmp[0])
#     f.close()

#     # canonical_chromosomes = set(str("chr" + str(i)) for i in range(1, 23)) | {'chrX', 'chrY'}
#     blacklisted_regions = load_bed_file(blacklist_bed)

#     breakreads = defaultdict(list)
#     bam_file = pysam.AlignmentFile(bam, "rb")
#     for read in bam_file:

#         ref_name = bam_file.get_reference_name(read.reference_id)
#         if ref_name not in analysis_chromosomes:
#             continue
#         if blacklisted_regions[ref_name][read.reference_start:read.reference_end]:
#             continue
#         if read.mapping_quality < 20:
#             continue
#         if not read.cigarstring:
#             continue
#         if read.is_reverse:
#             read.query_sequence = read.get_forward_sequence()
#         if re.search(r"^[0-9]+S([0-9]+I)?[0-9]+M$", read.cigarstring):
#             breakreads[read.qname].append(read)
#         if re.search(r"^[0-9]+M([0-9]+I)?[0-9]+S$", read.cigarstring):
#             breakreads[read.qname].append(read)            
#     return breakreads


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


def call_structural_variants(bam: str, bed: str, blacklist_bed: str, fasta: str):
    """ """

    bam_stats = get_bam_stats(bam, bed)

    breakreads = scan_breakreads(bam, bed, blacklist_bed, fasta)
    clustered_breakreads = cluster_breakreads(breakreads)


    # clusters = cluster_breakreads_naive(bam, breakreads, 500)
    # print(' INFO: Total clusters:', len(clusters))
    # breakreads_graph = create_breakreads_graph(breakreads, bam_stats['isize_median'], bam_stats['isize_mad'])
    # clusters = cluster_breakreads(breakreads_graph)
    variants = []
    for chromosome in clustered_breakreads:

        for cluster_position in clustered_breakreads[chromosome]:
            seqs = []
            pos = int(cluster_position)
            for read in clustered_breakreads[chromosome][cluster_position]:
                seqs.append(read.seq)

            oas = OverlapAssembler(seqs, 21)
            contigs = oas.compute_overlaps()

            # print(f" INFO: {chromosome}:{cluster_position} ### CONTIGS ###")
            # for idx,c in enumerate(contigs):
            #     print(idx, c)

            ref = pysam.FastaFile(fasta)
            reference = ref.fetch(chromosome, pos-250, pos+250)
            blat = Blat(k=21, ref=reference, chr=chromosome, start=pos-250, end=pos+250)

            for seq in contigs:
                vars_list = blat.align(seqs=[seq])
                for variant in vars_list:
                    if not variant in variants:
                        variants.append(variant)
    return variants

def export_sv_calls(sample, output_dir, variants_list):
    """ """
    bed_name = f'{sample.name}.GRAPES2.bed'
    bed_out = os.path.join(output_dir, bed_name)
    if not os.path.isfile(bed_out):
        o = open(bed_out, "w")
    else:
        o = open(bed_out, "a")
    for variant in variants_list:
        if int(variant["size"]) < 50:
            continue
        print(variant)
        call = f'{variant["chr"]}\t{variant["pos"]}\t{variant["end"]}\t{variant["size"]};{variant["type"]};{variant["contig"]}'
        o.write(call+"\n")
    o.close()


if __name__ == "__main__":
    bam = "/home/bdelolmo/Escriptori/NGS_APP/data/RUN20221118-CGC64001/SUDD_147/RB34017/BAM_FOLDER/RB34017.rmdup.bam"
    bam = "/home/bdelolmo/RB33925.rmdup.bam"
    bed = "/home/bdelolmo/BED/gendiag_85.CDS.bed"
    fasta = "/home/bdelolmo/REF_DIR/hg19/ucsc.hg19.fasta"

    # blat = Blat(reference, chr="chr2", start=179431346, end=179437577)

    call_structural_variants(bam, bed, fasta)
