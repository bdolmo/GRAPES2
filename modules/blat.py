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

    def __init__(self, ref: str, query_seq:str, k=11):
        self._ref = ref
        self._query_seq = query_seq
        self._k = k

        seq = Seq(ref)
        self._ref_rev = seq.reverse_complement()
        self._ref_index = self.build_index(self._ref)
        self._ref_rev_index = self.build_index(self._ref_rev)
        self._len_ref = len(ref)

    def extend_seeds(self):
        """ """

        if len(self.seeds) == 1:
            # print("single seed")
            init_seed = self.seeds[0]

            if init_seed.q_pos > 0:
                query_prefix = self._query_seq[0:init_seed.q_end+1]
                aligner = Aligner(self._ref, query_prefix)
                aln = aligner.affine_semiglobal_alignment()

                # print(aln["pretty_aln"])
                # Update seed 
                init_seed.cigar = aln["cigar"]
                init_seed.q_pos = aln["q_pos"]
                init_seed.r_pos = aln["r_pos"]
                init_seed.seq = query_prefix
                self.seeds[0] = init_seed

            if init_seed.q_end-init_seed.q_pos < len(self._query_seq)-1:
                query_prefix = self._query_seq[init_seed.q_pos:]
                aligner = Aligner(self._ref, query_prefix)
                aln = aligner.affine_semiglobal_alignment()
                # print(aln['pretty_aln'])

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
                    
            if next_seed.q_pos > curr_seed.q_end and next_seed.strand == curr_seed.strand:
                if idx == 0:
                    query_between_seeds = self._query_seq[0:next_seed.q_end]
                else:
                    query_between_seeds = self._query_seq[curr_seed.q_pos:next_seed.q_end]

                aligner = Aligner(self._ref, query_between_seeds)
                aln = aligner.affine_semiglobal_alignment()

                stack.remove(curr_seed)
                stack.remove(next_seed)

                next_seed.cigar = aln["cigar"]
                next_seed.q_pos = aln["q_pos"]
                next_seed.r_pos = aln["r_pos"]
                next_seed.seq = query_between_seeds
                extended_seeds.append(next_seed)
            else:
                extended_seeds.append(curr_seed)
                extended_seeds.append(next_seed)
                stack.remove(curr_seed)
                stack.remove(next_seed)
            idx+=1
 
        # Now extend initial seed
        if extended_seeds[0].q_pos > 0:
            query_suffix = self._query_seq[:extended_seeds[0].q_pos]

            aligner = Aligner(self._ref, query_suffix)
            aln = aligner.affine_semiglobal_alignment()

            next_seed.cigar = aln["cigar"]
            next_seed.q_pos = aln["q_pos"]
            next_seed.r_pos = aln["r_pos"]
            next_seed.seq = query_suffix
            extended_seeds[0] = next_seed

        # Now extend final seed
        if extended_seeds[-1].q_end < len(self._query_seq) - 1:
            query_suffix = self._query_seq[extended_seeds[-1].q_pos:]

            aligner = Aligner(self._ref, query_suffix)
            aln = aligner.affine_semiglobal_alignment()
            print(aln['pretty_aln'])
            next_seed.cigar = aln["cigar"]
            next_seed.q_pos = aln["q_pos"]
            next_seed.r_pos = aln["r_pos"]
            next_seed.seq = query_suffix
            extended_seeds[-1] = next_seed

        # self.seeds = self.merge_overlapping_seeds(self.seeds)
        # for seed in self.seeds:
        #     print(seed)
        print(extended_seeds)
            
    def merge_overlapping_seeds(self, seeds):
        merged_seeds = []
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

    def align(self, seqs):
        """ """
        for seq in seqs:
            print(f"query_seq: {seq} {len(seq)} ")
            chains = self.seed(seq)
            chains = self.extend(seq, chains)

            cigar = []

            if len(chains) == 1:
                print(chains)
                pass

            len_chains = len(chains)-1
            for i in range(0, len(chains) - 1):
                chain = chains[i]
                next_chain = chains[i + 1]
                ops = []
                cigar_dict = {}

                print(chain, next_chain, sep="\t")
                for n in range(0, len(chain.seq)):
                    chain_ntd = chain.seq[n]

                    ref_ntd = self._ref[chain.r_pos+n]
                    chain_ntd = chain.seq[n]

                    if chain.strand == "-":
                        pass
                        # ops.append("V")
                    else:
                        ref_ntd = self._ref[chain.r_pos+n]
                        if chain_ntd == ref_ntd:
                            ops.append("M")
                        else:
                            ops.append("X")

                if next_chain.diag != 0:

                    if next_chain.r_pos > chain.r_end+1:
                        deletion = self._ref[chain.r_end : next_chain.r_pos-1]

                        for n in range(0, len(deletion)):
                            ops.append("D")

                    if next_chain.r_pos == chain.r_end:

                        if next_chain.strand != chain.strand:
                            for n in range(0, len(deletion)):
                                ops.append("V")

                        insertion = seq[chain.q_end : next_chain.q_pos]
                        for n in range(0, len(insertion)):
                            ops.append("I")
                        print("INSERTION", insertion)
                if i == len(chains)-2:
                    for n in range(0, len(next_chain.seq)):
                        chain_ntd = next_chain.seq[n]
                        if next_chain.strand == "-":
                            ops.append("V")
                        else:
                            ref_ntd = self._ref[next_chain.r_pos+n]
                            if chain_ntd == ref_ntd:
                                ops.append("M")
                            else:
                                ops.append("X")

                cigar_dict = {
                    "M": 0,
                    "X": 0,
                    "D": 0,
                    "I": 0,
                    "V": 0
                }
                print(ops)
                for n in range(0, len(ops)-1):
                    curr_op = ops[n]
                    next_op = ops[n+1]
                    if next_op == curr_op:
                        cigar_dict[curr_op]+=1
                    if next_op != curr_op or n+1 == len(ops)-1:
                        print("fgggg")
                        cigar_dict[curr_op]+=1
                        cigar.append(str(cigar_dict[curr_op])+curr_op)
                        cigar_dict = {
                            "M": 0,
                            "X": 0,
                            "D": 0,
                            "I": 0,
                            "V": 0
                        }
            print(''.join(cigar))

    def seed(self, seq):
            """ """
            seeds = []
            for i in range(0, len(seq) - self._k + 1, self._k):
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
            # print(self.seeds)
            # sys.exit()

if __name__ == "__main__":

    # ref_seq   = Seq("tacgcccctagtacCCCCCAGAGAAGATAGTATGAAGGGGGGAGATATGAGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagct")
    # query_seq = Seq("tacgcccctagtacCCCCCAGAGAAGATAGTATGAAGGGGGGAGATATGAGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagct")

    # 1. This ons is inverted in the middle!
    # query_seq = Seq("tacgcccctagtacCCCCCAGAGAAGATATTTCTCATATCTCCCCCCTTCATACATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagct")

    # 2. This one has a simple deletion in the middle!
    # query_seq = Seq("tacgcccctagtacCCCCCAGAGAAGATAGTATcattagggggatgagaccccaggtagagagtacccccgggttatgxxagaaaaggcccc")

    ref_seq   = Seq("tacgcccctagtacCCCCCAGAGAAGATAGTATGAAGGGGGGAGATATGAGAAAATAGAATCCCACCTAGCCCGAGAcgatgtacgtgataaaaagct")
    query_seq = Seq("tacgcccctagtacCCCCCAGAGAAGATAGTATxxxxxxxxxxxxxxxxxxACCTAGCCCGAGAcgatgtacgtgataaaaagct")

    print("ref_seq:", ref_seq)
    print("query_seq:", query_seq)
    
    blat = Blat(k=11, ref=ref_seq, query_seq=query_seq)
    blat.seed(query_seq)
    blat.extend_seeds()
