import os
import re
import numpy as np


class Aligner:

    def __init__(self, seq1, seq2, match=2, mismatch=-5, gap_open=-10, gap_extend=-1):
        """ """
        self._seq1 = seq1.upper()
        self._seq2 = seq2.upper()
        self._match = match
        self._mismatch = mismatch
        self._gap_open = gap_open
        self._gap_extend = gap_extend

    @staticmethod
    def matprint(mat, fmt="g"):
        col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
        for x in mat:
            for i, y in enumerate(x):
                print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
            print("")

    @staticmethod
    def score_match(a, b, match_score, mismatch_score) -> int:
        return match_score if a == b else mismatch_score

    @staticmethod
    def compact_cigar_string(cigar_string) -> str:
        """Build a CIGAR string give a list of operations in linear space"""

        compact_cigar = ""
        current_char = cigar_string[0]
        count = 1

        for char in cigar_string[1:]:
            if char == current_char:
                count += 1
            else:
                compact_cigar += f"{count}{current_char}"
                current_char = char
                count = 1

        compact_cigar += f"{count}{current_char}"
        return compact_cigar

    def affine_semiglobal_alignment(self) -> dict:
        """ Semi-global alignment (ends-free alignment) with
            affine gap penalties
        """

        m = len(self._seq1)
        n = len(self._seq2)

        M = np.zeros((m+1, n+1))
        X = np.zeros((m+1, n+1))
        Y = np.zeros((m+1, n+1))

        # Initialize the first column of M, X, and Y
        for i in range(1, m+1):
            M[i, 0] = 0 
            X[i, 0] = self._gap_open + self._gap_extend-i
            Y[i, 0] = float('-inf')
         
        # Initialize the first row of M, X, and Y
        for j in range(1, n+1):
            M[0, j] = 0 
            X[0, j] = float('-inf')
            Y[0, j] = self._gap_open + self._gap_extend-j
        M[0][0] = 0
        X[0][0] = float('-inf')
        Y[0][0] = float('-inf')

        # filling matrices
        max_score = 0
        max_i = ""
        max_j = ""
        max_diag_score = 0
        xdrop_threshold = 50
        is_terminated = False

        for i in range(1, len(self._seq1)+1):
            for j in range(1, len(self._seq2)+1):
                M[i, j] = self.score_match(self._seq1[i-1], 
                            self._seq2[j-1], 
                            self._match, 
                            self._mismatch) + max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1])
                X[i, j] = max(M[i-1, j] + self._gap_open, 
                            X[i-1, j] + self._gap_extend, 
                            Y[i, j-1] + self._gap_open)
                Y[i, j] = max(M[i, j-1] + self._gap_open, 
                            Y[i, j-1] + self._gap_extend, 
                            X[i, j-1] + self._gap_open)

                # Update the maximum score in the current diagonal
                if M[i, j] > max_diag_score:
                    max_diag_score = M[i, j]

                # If the current score is below the X-drop threshold, set the score to -inf
                xdrop_threshold = abs(i-j)+self._match+40
                # print(xdrop_threshold)
                if max_diag_score - M[i, j] > xdrop_threshold:
                    M[i, j] = float('-inf')
                    X[i, j] = float('-inf')
                    Y[i, j] = float('-inf')
                    is_terminated = True

                max_value = max(M[i][j], X[i][j], Y[i][j])

                if max_value >= max_score:
                    max_score = max_value
                    max_i = i
                    max_j = j

        seq1_align = ""
        seq2_align = ""
        spacer = ""
        if is_terminated:
            i = max_i
            j = max_j
        cigar_str = ""

        # Traceback. Returning the single best optimal alignment
        while i > 0 or j > 0:
            # Check which matrix the current value came from
            current_matrix = None
            current_value = M[i][j]

            max_value = max(M[i, j], X[i, j], Y[i, j])

            if max_value == M[i, j]:
                current_matrix = "M"
            if max_value == X[i, j]:
                current_matrix = "X"
            if max_value == Y[i, j]:
                current_matrix = "Y"
           
            if max_value == float('-inf'):
                i -= 1
                j -= 1
                continue

            # Update alignment and indices based on the current matrix
            if current_matrix == "M":
                if i > 0 and j > 0:
                    seq1_align+=self._seq1[i-1]
                    seq2_align+=self._seq2[j-1]
                    if self._seq1[i-1] != self._seq2[j-1]:
                        spacer += "*"
                        cigar_str+="X"
                    else:
                        spacer += "|"
                        cigar_str+="M"
                    i -= 1
                    j -= 1
                else:
                    break
            elif current_matrix == "X":
                if i > 0:
                    seq2_align+="-"
                    seq1_align+=self._seq1[i-1]
                    spacer += " "
                    i -= 1
                    cigar_str+="D"
                else:
                    break
            elif current_matrix == "Y":
                if j > 0:
                    seq1_align+="-"
                    seq2_align+=self._seq2[j-1]
                    spacer += " "
                    j -= 1
                    cigar_str+="I"
                else:
                    break

        seq1_align = seq1_align[::-1]
        seq2_align = seq2_align[::-1]
        spacer = spacer[::-1]
        cigar_str = cigar_str[::-1]

        cigar = self.compact_cigar_string(cigar_str)
        cigar_operations = re.findall(r'\d+[MIDNSHP=X]', cigar)

        # Remove deletions (D) from the beginning and end of the CIGAR string
        if cigar_operations[0][-1] == 'D':
            cigar_operations.pop(0)
        if cigar_operations[-1][-1] == 'D':
            cigar_operations.pop()

        # Reconstruct the modified CIGAR string
        cigar = ''.join(cigar_operations)

        # offsets
        ref_end = max_i
        query_end = max_j
        ref_pos = i
        query_pos = j

        # alignment depiction
        pretty_aln = f"{seq1_align}\n{spacer}\n{seq2_align}"

        alignment = {
            "q_pos": query_pos,
            "q_end": query_end,
            "q_span": query_end-query_pos,
            "r_pos": ref_pos,
            "r_end": ref_end,
            "score": max_score,
            "cigar": cigar,
            "pretty_aln": pretty_aln
        }
        print(alignment)

        return alignment


    def affine_local_alignment(self) -> dict:
        """ Semi-global alignment (ends-free alignment) with
            affine gap penalties
        """

        m = len(self._seq1)
        n = len(self._seq2)

        M = np.zeros((m+1, n+1))
        X = np.zeros((m+1, n+1))
        Y = np.zeros((m+1, n+1))

        # Initialize the first column of M, X, and Y
        for i in range(1, m+1):
            M[i, 0] = 0 
            X[i, 0] = self._gap_open + self._gap_extend-i
            Y[i, 0] = float('-inf')
         
        # Initialize the first row of M, X, and Y
        for j in range(1, n+1):
            M[0, j] = 0 
            X[0, j] = float('-inf')
            Y[0, j] = self._gap_open + self._gap_extend-j
        M[0][0] = 0
        X[0][0] = float('-inf')
        Y[0][0] = float('-inf')

        # filling matrices
        max_score = 0
        max_i = ""
        max_j = ""
        max_diag_score = 0
        xdrop_threshold = 30
        is_terminated = False

        for i in range(1, len(self._seq1)+1):
            for j in range(1, len(self._seq2)+1):
                M[i, j] = max(
                    0,
                    self.score_match(self._seq1[i - 1], self._seq2[j - 1], self._match, self._mismatch)
                    + max(M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1]),
                )
                X[i, j] = max(0, M[i-1, j] + self._gap_open, 
                            X[i-1, j] + self._gap_extend, 
                            Y[i, j-1] + self._gap_open)
                Y[i, j] = max(0, M[i, j-1] + self._gap_open, 
                            Y[i, j-1] + self._gap_extend, 
                            X[i, j-1] + self._gap_open)

                # Update the maximum score in the current diagonal
                if M[i, j] > max_diag_score:
                    max_diag_score = M[i, j]

                # If the current score is below the X-drop threshold, set the score to -inf

                max_value = max(M[i][j], X[i][j], Y[i][j])

                if max_value >= max_score:
                    max_score = max_value
                    max_i = i
                    max_j = j

        seq1_align = ""
        seq2_align = ""
        spacer = ""
        i = max_i
        j = max_j
        cigar_str = ""


        # Traceback. Returning the single best optimal alignment
        while (i > 0 or j > 0) and (M[i][j] > 0 or X[i][j] > 0 or Y[i][j] > 0):
            # Check which matrix the current value came from
            current_matrix = None
            current_value = M[i][j]

            max_value = max(M[i, j], X[i, j], Y[i, j])

            if max_value == M[i, j]:
                current_matrix = "M"
            if max_value == X[i, j]:
                current_matrix = "X"
            if max_value == Y[i, j]:
                current_matrix = "Y"
           
            if max_value == float('-inf'):
                i -= 1
                j -= 1
                continue

            # Update alignment and indices based on the current matrix
            if current_matrix == "M":
                if i > 0 and j > 0:
                    seq1_align+=self._seq1[i-1]
                    seq2_align+=self._seq2[j-1]
                    if self._seq1[i-1] != self._seq2[j-1]:
                        spacer += "*"
                        cigar_str+="X"
                    else:
                        spacer += "|"
                        cigar_str+="M"
                    i -= 1
                    j -= 1
                else:
                    break
            elif current_matrix == "X":
                if i > 0:
                    seq2_align+="-"
                    seq1_align+=self._seq1[i-1]
                    spacer += " "
                    i -= 1
                    cigar_str+="D"
                else:
                    break
            elif current_matrix == "Y":
                if j > 0:
                    seq1_align+="-"
                    seq2_align+=self._seq2[j-1]
                    spacer += " "
                    j -= 1
                    cigar_str+="I"
                else:
                    break

        seq1_align = seq1_align[::-1]
        seq2_align = seq2_align[::-1]
        spacer = spacer[::-1]
        cigar_str = cigar_str[::-1]

        cigar = self.compact_cigar_string(cigar_str)
        
        cigar_operations = re.findall(r'\d+[MIDNSHP=X]', cigar)

        # Remove deletions (D) from the beginning and end of the CIGAR string
        if cigar_operations[0][-1] == 'D':
            cigar_operations.pop(0)
        if cigar_operations[-1][-1] == 'D':
            cigar_operations.pop()

        # Reconstruct the modified CIGAR string
        cigar = ''.join(cigar_operations)
        print(cigar)

        # offsets
        ref_end = max_i
        query_end = max_j
        ref_pos = i
        query_pos = j

        # alignment depiction
        pretty_aln = f"{seq1_align}\n{spacer}\n{seq2_align}"

        alignment = {
            "q_pos": query_pos,
            "q_end": query_end,
            "q_span": query_end-query_pos,
            "r_pos": ref_pos,
            "r_end": ref_end,
            "score": max_score,
            "cigar": cigar,
            "pretty_aln": pretty_aln
        }
        return alignment


if __name__ == "__main__":

    seq1 = "afaefageanufjmferh"
    seq2 = "afaefajmferh"

    aln = Aligner(seq1, seq2)
    alignment = aln.affine_semiglobal_alignment()
    print(alignment['pretty_aln'])
    print(alignment['score'])