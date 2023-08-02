import os
import sys
import random
# import edlib
from Bio import Align
from Bio.Seq import Seq
from modules.align import Aligner
from modules.affine.align import semiglobal_alignment, local_alignment


aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -5
aligner.open_gap_score = -20
aligner.extend_gap_score = -1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0

class Contig():

    def __init__(self, seq):
        self._seq = seq
        self._layout = []

class OverlapAssembler():
    """ """
    def __init__(self, reads, k):
        self.reads = list(set(reads))
        self.k = k

    def get_max_overlaps(self, reads):

        for read1 in reads:
            max_overlap = 0
            max_read = None
            contig = None
            for read2 in reads:
                if read1 == read2:
                    continue
               
                # if not read1 and not read2:
                #   continue
                # alignment = semiglobal_alignment(read1, read2) # New alignment function
                # if alignment["r_pos"] == 0 and alignment["q_end"] == len(read2):
                #     subseq1 = read2
                #     subseq2 = read1[alignment["r_end"]:]
                #     overlap = alignment["r_end"]
                # elif alignment["r_end"] == len(read1) and alignment["q_pos"] == 0:
                #     subseq2 = read2
                #     subseq1 = read1[0:alignment["r_pos"]]
                #     overlap = alignment["q_end"]
                # else:
                #     continue
                # if overlap > max_overlap and overlap > 20:
                #     max_overlap = overlap
                #     max_read = read2
                #     contig = subseq1+subseq2

                alignments = aligner.align(read1, read2)
                optimal_aln = alignments[0]
                seg1 = optimal_aln.aligned[0]
                seg2 = optimal_aln.aligned[1]
                if len(seg1) == 1 and len(seg2) == 1:
                    # print(optimal_aln, seg1,seg1[0], seg2, seg2[0])
                    if seg1[0][0] == 0 and seg2[0][1] == len(read2):
                        subseq1 = read2
                        subseq2 = read1[seg1[0][1]:]
                        overlap = seg1[0][1]
                    elif seg1[0][1] == len(read1) and seg2[0][0] == 0:
                        subseq2 = read2
                        subseq1 = read1[0:seg1[0][0]]
                        overlap = seg2[0][1]
                    else:
                        continue
                    if overlap > max_overlap and overlap > 20:
                        max_overlap = overlap
                        max_read = read2
                        contig = subseq1+subseq2

            return read1, max_read, max_overlap, contig

    def compute_overlaps(self):
        reads = self.reads

        bag = []
        for read in reads:
            c_read = Contig(read)
            bag.append(c_read)

        contigs = []
        while reads:
            read_a, read_b, olen, contig = self.get_max_overlaps(reads)
            is_rev = False
            if not read_b:
                seq = Seq(read_a)
                reads.remove(read_a)
                read_a = seq.reverse_complement()
                reads.insert(0, read_a)
                read_a, read_b, olen, contig = self.get_max_overlaps(reads)

            reads.remove(read_a)
            if not read_b:
                contigs.append(read_a)
            else:
                reads.remove(read_b)
                reads.append(contig)
            if not reads:
                break
 
        return contigs

class Node():
  '''
  '''
  def __init__(self, km1mer):
    self.km1mer = km1mer
    self.nin = 0
    self.nout = 0

class DeBruijnAssembler():
  '''
  '''
  def __init__(self, reads, k):
    self.reads = reads
    self.k = k
    self.g = self.build_graph()

  def build_graph(self) -> dict:
    '''
      Build a de Bruijn Graph
    '''

    G = dict()
    nodes = dict()

    # unique_reads = set(self.reads)
    kmer_distribution = {}
    for read in self.reads:
      i = 0
      for i in range(0, len(read)-self.k):
        km1_l = read[i:i+self.k]
        km1_r = read[i+1:i+self.k+1]
        node_l = None
        node_r = None
        if km1_l in nodes:
          node_l = nodes[km1_l]
          node_l.nout+=1
        else:
          node_l = Node(km1_l)
          node_l.nout+=1
          nodes[km1_l] = node_l
          G[node_l] = list()
        if km1_r in nodes:
          node_r = nodes[km1_r]
          node_r.nin+=1
        else:
          node_r = Node(km1_r)
          node_r.nin+=1
          nodes[km1_r] = node_r
          G[node_r] = list()
        if node_r not in G[node_l]:
          G[node_l].append(node_r)

    return G

  def eulerian_walk(self):
    '''
    '''
    contig_list = list()

    if len(self.g.keys()) == 0:
      return ""


    start = list(self.g.keys())[0]
    for node in self.g:
      for subnode in self.g[node]:
        if node.nin < start.nin:
          start = node
    current = start
    contig = current.km1mer[:-1]

    tour = []
    def _visit(current):
      while len(self.g[current]) > 0:
        dst = self.g[current].pop()
        _visit(dst)
      tour.append(current)

    _visit(current)
    tour = tour[::-1]
    for n in tour:
      if n.nout == 0:
        contig+=n.km1mer[-1]
        contig_list.append(contig)
        contig = ""
      else:
        contig+=n.km1mer[-1]

    return contig_list
