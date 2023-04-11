import os
import sys
# import edlib
from Bio import pairwise2, Align
from Bio.Seq import Seq

class Seed:
    def __init__(self, q_pos, q_end, r_pos, r_end, strand, seq):
        self.q_pos = q_pos
        self.q_end = q_end
        self.r_pos = r_pos
        self.r_end = r_end
        self.strand = strand
        self.seq = seq

    @property
    def diag(self):
        """ """
        if self.strand == "+":
            self._diag = self.r_pos - self.q_pos
        else:
            self._minusqpos = self.q_pos-1
            self._minusrpos = self.q_end
            self._diag = self._minusrpos-self._minusqpos

        return self._diag

    def __repr__(self):
        return (
            f"Seed(q_pos={self.q_pos}, q_end={self.q_end}, r_pos={self.r_pos}"
            f", r_end={self.r_end}, strand={self.strand}, seq={self.seq})"
        )


class BLAT:
    def __init__(self, k, ref_seq):
        self._k = k
        self._ref_seq = ref_seq
        self._len_ref = len(ref_seq)
        self._max_hits = 1
        self._ref_index, self._ref_rev_index = self._build_ref_indices()

    def _build_ref_indices(self):
        """ """
        ref_index = {}
        ref_rev_index = {}
        ref_seq_rc = self._ref_seq.reverse_complement()

        for i in range(self._len_ref - self._k + 1):
            kmer = self._ref_seq[i : i + self._k]
            rev_kmer = ref_seq_rc[i : i + self._k]

            if kmer in ref_index:
                ref_index[kmer].append(i)
            else:
                ref_index[kmer] = [i]

            if rev_kmer in ref_rev_index:
                ref_rev_index[rev_kmer].append(i)
            else:
                ref_rev_index[rev_kmer] = [i]

        return ref_index, ref_rev_index

    def seed(self, query_seq):
        """ """
        seeds = []
        for i in range(0, len(query_seq) - self._k + 1):
            kmer = query_seq[i : i + self._k]
            if kmer in self._ref_index:
                if len(self._ref_index[kmer]) > self._max_hits:
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
                    if len(self._ref_rev_index[kmer]) > self._max_hits:
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

        idx = 0
        for i in range(1, len(seeds)):
            seed1 = seeds[idx]
            seed2 = seeds[i]
            if seed2.q_pos <= seed1.q_end+1 and seed2.diag == seed1.diag:
                seed1.q_end = seed2.q_end
                seed1.r_end = seed2.r_end
                seed1.seq = query_seq[seed1.q_pos:seed1.q_end+1]
            else:
                idx += 1
                seeds[idx] = seed2

        for i in range(0, idx+1):
            seed = seeds[i]
            if seed.strand == "-":
                rpos = seed.r_pos
                rend = seed.r_end
                seed.r_pos = rend
                seed.r_end = rpos
            chains.append(seed)

        return chains

    def extend_gapped(self, ref_seq, query_seq, seeds, match=1, mismatch=-1, 
        gap_open=-1, gap_extend=-1, x_dropoff=10):
        """ """
        extended_alignments = []

        for seed in seeds:
            q_left = seed.q_pos
            q_right = len(query_seq) - seed.q_end - 1
            r_left = seed.r_pos
            r_right = len(ref_seq) - seed.r_end - 1

            if seed.strand == "-":
                query_seq = query_seq.reverse_complement()

            # print(q_left, q_right)
            # print("seq", seed.seq, ref_seq[:r_left][::-1], query_seq[:q_left][::-1])


            # Extend left
            left_score, left_path = self.extend_gapped_x_dropoff(
                ref_seq[:r_left][::-1],
                query_seq[:q_left][::-1],
                x_dropoff,
                match,
                mismatch,
                gap_open,
                gap_extend,
            )

            # Extend right
            right_score, right_path = self.extend_gapped_x_dropoff(
                ref_seq[seed.r_end + 1:],
                query_seq[seed.q_end + 1:],
                x_dropoff,
                match,
                mismatch,
                gap_open,
                gap_extend,
            )

            # Combine left and right extensions with the seed
            alignment = {
                "q_start": seed.q_pos - len(left_path),
                "q_end": seed.q_end + len(right_path),
                "r_start": seed.r_pos - len(left_path),
                "r_end": seed.r_end + len(right_path),
                "strand": seed.strand,
                "score": left_score + right_score + 1,
                "alignment": left_path[::-1] + seed.seq + right_path,
            }

            extended_alignments.append(alignment)

            if seed.strand == "-":
                query_seq = query_seq.reverse_complement()

        return extended_alignments

    def extend_gapped_x_dropoff(self, ref_tail, query_tail, x_dropoff, match, mismatch, gap_open, gap_extend):
        """ """
        H = [[0] * (len(query_tail) + 1) for _ in range(len(ref_tail) + 1)]
        E = [[0] * (len(query_tail) + 1) for _ in range(len(ref_tail) + 1)]
        F = [[0] * (len(query_tail) + 1) for _ in range(len(ref_tail) + 1)]

        max_score = 0
        max_pos = (0, 0)

        for i in range(1, len(ref_tail) + 1):
            for j in range(1, len(query_tail) + 1):
                match_score = H[i - 1][j - 1] + (match if ref_tail[i - 1] == query_tail[j - 1] else mismatch)
                E[i][j] = max(E[i][j - 1] + gap_extend, H[i][j - 1] + gap_open + gap_extend)
                F[i][j] = max(F[i - 1][j] + gap_extend, H[i - 1][j] + gap_open + gap_extend)
                H[i][j] = max(0, match_score, E[i][j], F[i][j])

                if H[i][j] > max_score:
                    max_score = H[i][j]
                    max_pos = (i, j)

        if max_score == 0:
            return 0, ""

        i, j = max_pos
        path = []

        while H[i][j] != 0:
            if max_score - H[i][j] > x_dropoff:
                break

            if H[i][j] == H[i - 1][j - 1] + (match if ref_tail[i - 1] == query_tail[j - 1] else mismatch):
                path.append(ref_tail[i - 1])
                i -= 1
                j -= 1
            elif H[i][j] == E[i][j]:
                path.append('-')
                j -= 1
            else:
                path.append(ref_tail[i - 1].lower())
                i -= 1

        return max_score, "".join(path)


if __name__ == "__main__":
    ref_seq = Seq("GGAGATATGAGAAAATAGAATCCCACCTAGCCCATTAGCAGCGAGTGTGTACCGAATAGCCGAC")
    query_seq = Seq("GGAGGAGATATGAGAAAATAGAATCCCACCTAGCCCATTAGCAGCGAGTGTGTACCGAATAGCCGACGATATGAGAAAATAGAGAGTGTGTACCGAATAGCCGAC")

    print("ref_seq:", ref_seq)
    print("query_seq:", query_seq)
    
    blat = BLAT(k=11, ref_seq=ref_seq)
    seed_chains = blat.seed(query_seq)

    # print("Seed chains:")
    # for chain in seed_chains:
    #     print(chain)

    extended_alignments = blat.extend_gapped(ref_seq, query_seq, seed_chains)

    print("\nExtended alignments:")
    for alignment in extended_alignments:
        print(alignment)

