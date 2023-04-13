import os
import sys
# import edlib
# from Bio import Align
from align import Aligner
from Bio.Seq import Seq
import re


class Seed:
    """ """

    def __init__(self, q_pos, q_end, r_pos, r_end, seq, strand):
        self.q_pos = q_pos
        self.q_end = q_end
        self.r_pos = r_pos
        self.r_end = r_end
        self.seq = seq
        self.strand = strand
        self.cigar = ""

        self.seeds = []

        self.visited = False
        if self.strand == "+":
            self.diag = r_pos - q_pos
        else:
            self._minusqpos = q_pos-1
            self._minusrpos = q_end
            self.diag = self._minusrpos-self._minusqpos

    def __repr__(self):
        back = (
            f"q_pos:{self.q_pos},q_end:{self.q_end} {(self.seq)} "
            f"r_pos:{self.r_pos},r_end:{self.r_end} strand:{self.strand} cigar:{self.cigar} diag:{self.diag}"
        )
        return back

class Blat:
    """Blat-like implementation"""

    def __init__(self, ref: str, chr: str, start:int, end: int,  k=11):
        self._ref = ref
        # self._query_seq = query_seq
        self._k = k
        self._chr = chr
        self._start = start
        self._end = end
        self._seeds = []

        seq = Seq(ref)
        self._ref_rev = seq.reverse_complement()
        self._ref_index = self.build_index(self._ref)
        self._ref_rev_index = self.build_index(self._ref_rev)
        self._len_ref = len(ref)

    def extend_seeds(self, seq):
        """ """

        if len(self.seeds) == 1:
            init_seed = self.seeds[0]
            # print("single seed")
            if init_seed.q_pos > 0:
                query_prefix = seq[0:init_seed.q_end+1]
                aligner = Aligner(self._ref, query_prefix)
                aln = aligner.affine_semiglobal_alignment()

                # Update seed 
                init_seed.cigar = aln["cigar"]
                init_seed.q_pos = aln["q_pos"]
                init_seed.r_pos = aln["r_pos"]
                init_seed.seq = query_prefix
                self.seeds[0] = init_seed

            if init_seed.q_end-init_seed.q_pos < len(seq)-1:
                query_prefix = seq[init_seed.q_pos:]
                aligner = Aligner(self._ref, query_prefix)
                aln = aligner.affine_semiglobal_alignment()

                # Update seed 
                init_seed.cigar = aln["cigar"]
                init_seed.q_pos = aln["q_pos"]
                init_seed.q_end = aln["q_end"]
                init_seed.r_pos = aln["r_pos"]
                init_seed.r_end = aln["r_end"]
                init_seed.seq = query_prefix
                self.seeds[0] = init_seed
            return self.seeds

        extended_seeds = []

        stack = self.seeds
        idx = 0

        while len(stack) > 1:
            curr_seed = stack[0]
            next_seed = stack[1]
     
            if next_seed.q_pos > curr_seed.q_end+1:

                if idx == 0:
                    query_between_seeds = seq[0:next_seed.q_end]
                else:
                    query_between_seeds = seq[curr_seed.q_pos:next_seed.q_end]

                ref = self._ref
                if next_seed.strand != curr_seed.strand:
                    ref = self._ref_rev

                aligner = Aligner(ref, query_between_seeds)
                if next_seed.strand != curr_seed.strand:
                    aln = aligner.affine_local_alignment()
                else:
                   aln = aligner.affine_semiglobal_alignment()
    

                new_seed = next_seed
                new_seed.cigar = aln["cigar"]
                new_seed.q_pos = aln["q_pos"]
                new_seed.q_end = aln["q_end"]
                new_seed.r_pos = aln["r_pos"]
                new_seed.r_end = aln["r_end"]
                new_seed.seq = seq[new_seed.q_pos:new_seed.q_end]

                if new_seed.strand != curr_seed.strand:
                    extended_seeds.append(curr_seed)
                    extended_seeds.append(new_seed)
                    stack.remove(curr_seed)
                    stack.remove(next_seed)
                else:
                    extended_seeds.append(new_seed)
                    extended_seeds.append(curr_seed)

                    stack.remove(curr_seed)
                    stack.remove(next_seed)
            else:
                new_seed = next_seed

                extended_seeds.append(curr_seed)
                extended_seeds.append(new_seed)
                stack.remove(curr_seed)
                stack.remove(next_seed)
            idx+=1
 
        if len(stack) == 1:
            extended_seeds.append(stack[0])

        # Extend the initial seed
        if extended_seeds[0].q_pos > 0:

            query_suffix = seq[0:extended_seeds[0].q_end]
          
            aligner = Aligner(self._ref, query_suffix)

            if extended_seeds[1].strand != extended_seeds[0].strand:
                aln = aligner.affine_local_alignment()
            else:
                aln = aligner.affine_semiglobal_alignment()

            new_seed = extended_seeds[0]
            new_seed.cigar = aln["cigar"]
            new_seed.q_pos = aln["q_pos"]
            new_seed.r_pos = aln["r_pos"]
            new_seed.seq = query_suffix
            extended_seeds[0] = new_seed

        # Now extend the final seed
        if extended_seeds[-1].q_end < len(seq) - 1:
            query_suffix = seq[extended_seeds[-1].q_pos:]

            aligner = Aligner(self._ref, query_suffix)

            if extended_seeds[-1].strand != extended_seeds[-2].strand:
                aln = aligner.affine_local_alignment()
            else:
                aln = aligner.affine_semiglobal_alignment()

            new_seed = extended_seeds[-1]
            new_seed.cigar = aln["cigar"]
            new_seed.q_pos = aln["q_pos"]
            new_seed.r_pos = aln["r_pos"]
            new_seed.seq = seq[extended_seeds[-1].q_pos:extended_seeds[-1].q_end]
            extended_seeds[-1] = new_seed

        # print(extended_seeds)

        # extended_seeds = self.merge_overlapping_seeds(extended_seeds)
        return extended_seeds
            
    def merge_overlapping_seeds(self, seeds):
        """ """
        merged_seeds = []

        if len(seeds) <= 1:
            return seeds

        sorted_seeds = sorted(seeds, key=lambda x: x.q_pos, reverse=False)

        current_seed = sorted_seeds[0]

        for next_seed in sorted_seeds[1:]:

            # print(current_seed, next_seed)
            if current_seed.q_end >= next_seed.q_pos and current_seed.q_end >= next_seed.q_pos:
                # Merge overlapping seeds
                current_seed.seq = current_seed.seq[0:] + next_seed.seq[current_seed.q_end:]
                current_seed.q_end = max(current_seed.q_end, next_seed.q_end)
                current_seed.r_end = max(current_seed.r_end, next_seed.r_end)
                current_seed.cigar = self.merge_and_simplify_cigars(current_seed.cigar, next_seed.cigar)
            else:
                # Add the current seed to the merged list
                merged_seeds.append(current_seed)
                current_seed = next_seed

        # Add the last seed to the merged list
        merged_seeds.append(current_seed)

        return merged_seeds

    def merge_and_simplify_cigars(self, cigar1, cigar2):
        """ """
        cigar_operations1 = re.findall(r'\d+[MIDNSHP=X]', cigar1)
        cigar_operations2 = re.findall(r'\d+[MIDNSHP=X]', cigar2)

        combined_cigar_operations = cigar_operations1 + cigar_operations2
        simplified_operations = []

        prev_op = combined_cigar_operations[0]
        count = int(prev_op[:-1])

        for curr_op in combined_cigar_operations[1:]:
            if curr_op[-1] == prev_op[-1]:
                count += int(curr_op[:-1])
            else:
                simplified_operations.append(f"{count}{prev_op[-1]}")
                prev_op = curr_op
                count = int(curr_op[:-1])

        simplified_operations.append(f"{count}{prev_op[-1]}")

        return ''.join(simplified_operations)

    def build_index(self, ref):
        """ """
        ref_index = {}
        for i in range(0, len(ref) - self._k + 1):
            kmer = ref[i : i + self._k]
            if not kmer in ref_index:
                ref_index[kmer] = []
            ref_index[kmer].append(i)
        return ref_index

    def complement(self, ntd):
        complement_dict = {
            "A" : "T",
            "a" : "t",
            "T" : "A",
            "t" : "a",
            "C" : "G",
            "c" : "g",
            "G" : "C",
            "g" : "c"
        }

        return complement_dict[ntd]

    @staticmethod
    def cigar_to_variant(seed):
        """Given a list of cigar operations extract a variant"""
        cigar_op_list = re.findall(r'\d+[MIDNSHP=X]', seed.cigar)
        var_list = []
        cumulative_pos = 0
        for op in cigar_op_list:
            tmp = re.split(r"[XMID]", op)
            if op.endswith("M"):
                cumulative_pos+=int(tmp[0])
                continue
            if op.endswith("X"):
                cumulative_pos+=int(tmp[0])
                continue
            if op.endswith("D"):
                vartype = "deletion"
                end = cumulative_pos+int(tmp[0])+seed.r_pos
                seq = ""
            if op.endswith("I"):
                vartype = "insertion"
                end = cumulative_pos+seed.r_pos
                seq = ""
            var_dict = {
                "type": vartype,
                "size": tmp[0],
                "pos": cumulative_pos+seed.r_pos,
                "end": end,
                "seq": seq
            }
            cumulative_pos+=int(tmp[0])
            var_list.append(var_dict)
        return var_list

    def align(self, seqs):
        """ """
        var_list = []
        for seq in seqs:
            self.seed(seq)
            seeds = self.extend_seeds(seq)

            if len(seeds) == 1:
                # print(seeds[0].cigar)
                var_list = self.cigar_to_variant(seeds[0])
                    
            for i in range(0, len(seeds)-1):
                curr_seed = seeds[i]
                next_seed = seeds[i+1]
                # print(curr_seed, next_seed)

                vars1 = self.cigar_to_variant(curr_seed)
                vars2 = self.cigar_to_variant(next_seed)
                var_dict = {}

                if vars1 and vars1 not in var_list:
                    var_list.append(vars1)
                if vars2 and vars2 not in var_list:
                    var_list.append(vars2)

                if next_seed.r_pos > curr_seed.r_end+1:
                    deletion_size = next_seed.r_pos - curr_seed.r_end+1
                    var_dict =  {
                        "type": "deletion",
                        "size" : deletion_size,
                        "chr": self._chr,
                        "pos": curr_seed.r_end+1+self._start,
                        "end": next_seed.r_pos+self._start,
                        "alt": self._ref[curr_seed.r_end+1:next_seed.r_pos]
                    }
                else:
                    if next_seed.strand != curr_seed.strand:
                        inv_seed = curr_seed
                        if next_seed.strand == "-":
                            inv_seed = next_seed

                        inversion_size = inv_seed.r_end - inv_seed.r_pos
                        var_dict =  {
                            "type": "inversion",
                            "size" : inversion_size,
                            "chr": self._chr,
                            "pos": inv_seed.r_pos+self._start,
                            "end": inv_seed.r_end+self._start,
                            "alt": self._ref[inv_seed.r_pos:inv_seed.r_end]
                        }

                if var_dict and var_dict not in var_list:
                    var_list.append(var_dict)

        for var in var_list:
            print(var)

    def seed(self, seq):
        """ """
        seeds = []
        for i in range(0, len(seq) - self._k + 1, 1):
            kmer = seq[i : i + self._k]
            if kmer in self._ref_index:
                if len(self._ref_index[kmer]) > 1:
                    continue
                for hit in self._ref_index[kmer]:
                    seed = Seed(
                        q_pos=i,
                        q_end=i + self._k-1,
                        r_pos=hit,
                        r_end=hit + self._k-1,
                        strand="+",
                        seq=kmer,
                    )
                    seeds.append(seed)
            else:
                if kmer in self._ref_rev_index:
                    if len(self._ref_rev_index[kmer]) > 1:
                        continue

                    kseq = Seq(kmer)
                    fwd_kmer = kseq.reverse_complement()

                    for hit in self._ref_rev_index[kmer]:
                        seed = Seed(
                            q_pos=i,
                            q_end=i + self._k-1,
                            r_pos=self._len_ref-hit,
                            r_end=self._len_ref-hit - self._k-1,
                            strand="-",
                            seq=fwd_kmer,
                        )
                        seeds.append(seed)
        chains = []
        if not seeds:
            return chains

        # merge consecutive seeds
        idx = 0
        for i in range(1, len(seeds)):
            if seeds[i].q_pos <= seeds[idx].q_end+1 and seeds[i].diag == seeds[idx].diag:
                seeds[idx].q_end = seeds[i].q_end
                seeds[idx].r_end = seeds[i].r_end
                seeds[idx].seq = seq[seeds[idx].q_pos:seeds[idx].q_end+1]
                seeds[idx].cigar = f"{str(len(seeds[idx].seq))}M"
            else:
                idx += 1
                seeds[idx] = seeds[i]

        for i in range(0, idx+1):
            if seeds[i].strand == "-":
                rpos = seeds[i].r_pos
                rend = seeds[i].r_end
                seeds[i].r_pos = rend
                seeds[i].r_end = rpos
                seeds[i].cigar = f"{str(len(seeds[i].seq))}M"
            chains.append(seeds[i])
        schains = sorted(chains, key=lambda x: x.q_pos, reverse=False)
        self.seeds = schains

        # Fix positions if applicable
        for i in range(0, len(self.seeds)-1):
            curr_seed = self.seeds[i]
            next_seed = self.seeds[i+1]
            next_seed.r_pos+=1
            if curr_seed.strand == next_seed.strand:
                if curr_seed.q_end > next_seed.q_pos:
                    offset = curr_seed.q_end-next_seed.q_pos
                    next_seed.q_pos = curr_seed.q_end+1
                    next_seed.r_pos = next_seed.r_pos+offset
                    next_seed.seq = next_seed.seq[offset+1:]
                    next_seed.cigar = f"{str(len(next_seed.seq))}M"

        return self.seeds
  

if __name__ == "__main__":

    ref_seq   = Seq("ttttaggacggtacgagacagcgtttttttagtacgcccctagtacCCCCCAGAGAAGATAGTATGAAGGGGGGAGATATGAGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagctacccattaggggaaaagttttaccgt")

    # 1. This ons is inverted in the middle!
    query_seq = Seq("acgagacagcgtttttttagtacgcccctagtacCCCCCAxxxxxxagGAGAAGATAGTATGAAGGGGGGAGATATGAGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagctacccattaggggaaaagttttaccgt")

    # query_seq = Seq("tacgcccctagtacCCCCCAGAGAAGATAGGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagct")
    print("ref_seq:", ref_seq)
    print("query_seq:", query_seq)
    
    blat = Blat(k=11, ref=ref_seq)
    blat.align(seqs=[query_seq])
