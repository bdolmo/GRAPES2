import numpy as np
import re
cimport numpy as np
cimport cython
from libc.math cimport INFINITY as INF

def score_match(a, b, match_score, mismatch_score) -> int:
    return match_score if a == b else mismatch_score

def STRINGIFY_CIGAR(cigar_string) -> str:
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


def compact_cigar_string(sequence):
    cigar = []
    count = 1
    prev_char = sequence[0]

    for char in sequence[1:]:
        if char == prev_char:
            count += 1
        else:
            cigar.append(f"{count}{prev_char}")
            count = 1
            prev_char = char

    cigar.append(f"{count}{prev_char}")
    return "".join(cigar)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def semiglobal_alignment(seq1, seq2, gap_open=-10, gap_extend=-1, match=2, mismatch=-1, band_width=30) -> dict:
    """ Semi-global alignment (ends-free alignment) with
        affine gap penalties
    """
    cdef:
        int m = len(seq1)
        int n = len(seq2)
        int i, j, max_i, max_j, j_lower_bound, j_upper_bound
        double[:, :] M = np.zeros((m + 1, n + 1), dtype=np.float64)
        double[:, :] X = np.zeros((m + 1, n + 1), dtype=np.float64)
        double[:, :] Y = np.zeros((m + 1, n + 1), dtype=np.float64)
        double max_score, max_diag_score, xdrop_threshold
        bint is_terminated
        str seq1_align, seq2_align, spacer, cigar_str, cigar
        int seq1_len, seq2_len
        object cigar_operations
        double neg_inf = -INF

    for i in range(m + 1):
        if i == 0:
            X[i, 0] = 0
            Y[0, i] = 0
            M[i, 0] = 0
            M[0, i] = 0
            X[0, i] = neg_inf
            Y[i, 0] = neg_inf
        else:
            X[i, 0] = gap_open + gap_extend * i
            Y[0, i] = gap_open + gap_extend * i
            M[i, 0] = 0
            M[0, i] = 0
            X[0, i] = neg_inf
            Y[i, 0] = neg_inf

    X[0, 0] = neg_inf
    Y[0, 0] = neg_inf

    seq1_len = len(seq1)
    seq2_len = len(seq2)

    for i in range(1, seq1_len + 1):

        # for j in range(1, seq2_len + 1):

        j_lower_bound = max(1, i - band_width)
        j_upper_bound = min(seq2_len + 1, i + band_width+1)
        # for j in range(j_lower_bound, j_upper_bound):
            
        for j in range(1, seq2_len + 1):
            M[i, j] = score_match(seq1[i-1], 
                        seq2[j-1], 
                        match, 
                        mismatch) + max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1])
            X[i, j] = max(M[i-1, j] + gap_open, 
                        X[i-1, j] + gap_extend, 
                        Y[i, j-1] + gap_open)
            Y[i, j] = max(M[i, j-1] + gap_open, 
                        Y[i, j-1] + gap_extend, 
                        X[i, j-1] + gap_open)

            # Update the maximum score in the current diagonal
            if M[i, j] > max_diag_score:
                max_diag_score = M[i, j]

            # If the current score is below the X-drop threshold, set the score to -inf
            xdrop_threshold = abs(i-j)+match+30
            # print(xdrop_threshold)
            if max_diag_score - M[i, j] > xdrop_threshold:
                M[i, j] = float('-inf')
                X[i, j] = float('-inf')
                Y[i, j] = float('-inf')
                is_terminated = True

            # else:
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
                seq1_align+=seq1[i-1]
                seq2_align+=seq2[j-1]
                if seq1[i-1] != seq2[j-1]:
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
                seq1_align+=seq1[i-1]
                spacer += " "
                i -= 1
                cigar_str+="D"
            else:
                break
        elif current_matrix == "Y":
            if j > 0:
                seq1_align+="-"
                seq2_align+=seq2[j-1]
                spacer += " "
                j -= 1
                cigar_str+="I"
            else:
                break

    seq1_align = seq1_align[::-1]
    seq2_align = seq2_align[::-1]
    spacer = spacer[::-1]
    cigar_str = cigar_str[::-1]

    if not cigar_str:
        alignment = {
            "q_pos": 0,
            "q_end": 0,
            "q_span": 0,
            "r_pos": 0,
            "r_end": 0,
            "score": 0,
            "cigar": "0M",
            "pretty_aln": "."
        }
        return alignment


    cigar = compact_cigar_string(cigar_str)
    cigar_operations = re.findall(r'\d+[MIDNSHP=X]', cigar)
    # print(cigar)

    # Remove deletions (D) from the beginning and end of the CIGAR string

    if isinstance(cigar_operations, list) and len(cigar_operations) > 0 and isinstance(cigar_operations[0], list) and len(cigar_operations[0]) > 0:
        # The list is multidimensional, and the expression is valid.
        if cigar_operations[0][-1] == 'D':
            cigar_operations.pop(0)


    if isinstance(cigar_operations, list) and len(cigar_operations) > 0 and isinstance(cigar_operations[0], list) and len(cigar_operations[0]) > 0:
        # The list is multidimensional, and the expression is valid.
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

    return alignment

def local_alignment(seq1, seq2, gap_open=-10, gap_extend=-1, match=2, mismatch=-4, band_width=50) -> dict:
    """
        Local alignment with affine gap penalties
    """
    cdef:
        int m = len(seq1)
        int n = len(seq2)
        int i, j, max_i, max_j, j_lower_bound, j_upper_bound
        double[:, :] M = np.zeros((m + 1, n + 1), dtype=np.float64)
        double[:, :] X = np.zeros((m + 1, n + 1), dtype=np.float64)
        double[:, :] Y = np.zeros((m + 1, n + 1), dtype=np.float64)
        double max_score, max_diag_score, xdrop_threshold
        bint is_terminated
        str seq1_align, seq2_align, spacer, cigar_str, cigar
        int seq1_len, seq2_len
        object cigar_operations
        double neg_inf = -INF

        # Initialize the first column of M, X, and Y

    for i in range(1, m+1):
        M[i, 0] = 0 
        X[i, 0] = gap_open + gap_extend-i
        Y[i, 0] = float('-inf')
        
    # Initialize the first row of M, X, and Y
    for j in range(1, n+1):
        M[0, j] = 0 
        X[0, j] = float('-inf')
        Y[0, j] = gap_open + gap_extend-j
    M[0][0] = 0
    X[0][0] = float('-inf')
    Y[0][0] = float('-inf')


    for i in range(1, len(seq1)+1):

        j_lower_bound = max(1, i - band_width)
        j_upper_bound = min(n + 1, i + band_width+1)
        # for j in range(j_lower_bound, j_upper_bound):
        for j in range(1, len(seq2)+1):
            # if abs(i - j) <= band_width:
            M[i, j] = max(
                0,
                score_match(seq1[i - 1], seq2[j - 1], match, mismatch)
                + max(M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1]),
            )
            X[i, j] = max(0, M[i-1, j] + gap_open, 
                        X[i-1, j] + gap_extend, 
                        Y[i, j-1] + gap_open)
            Y[i, j] = max(0, M[i, j-1] + gap_open, 
                        Y[i, j-1] + gap_extend, 
                        X[i, j-1] + gap_open)

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
    cigar_str = ""


    i = max_i 
    j = max_j
    # print(i, j)


    if i == 0 and j == 0:
        alignment = {
            "q_pos": 0,
            "q_end": 0,
            "q_span": 0,
            "r_pos": 0,
            "r_end": 0,
            "score": 0,
            "cigar": "0M",
            "pretty_aln": "."
        }
        return alignment


    # Traceback. Returning the single best optimal alignment
    #while (i > 0 or j > 0) and (M[i][j] > 0 or X[i][j] > 0 or Y[i][j] > 0):
    while (i > 0 and i < m + 1 and j > 0 and j < n + 1) and (M[i][j] > 0 or X[i][j] > 0 or Y[i][j] > 0):

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
                seq1_align+=seq1[i-1]
                seq2_align+=seq2[j-1]
                if seq1[i-1] != seq2[j-1]:
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
                seq1_align+=seq1[i-1]
                spacer += " "
                i -= 1
                cigar_str+="D"
            else:
                break
        elif current_matrix == "Y":
            if j > 0:
                seq1_align+="-"
                seq2_align+=seq2[j-1]
                spacer += " "
                j -= 1
                cigar_str+="I"
            else:
                break

    seq1_align = seq1_align[::-1]
    seq2_align = seq2_align[::-1]
    spacer = spacer[::-1]
    cigar_str = cigar_str[::-1]

    # print(seq1_align)
    # print(spacer)
    # print(seq2_align)
    if not cigar_str:
        alignment = {
            "q_pos": 0,
            "q_end": 0,
            "q_span": 0,
            "r_pos": 0,
            "r_end": 0,
            "score": 0,
            "cigar": "0M",
            "pretty_aln": "."
        }
        return alignment
   

    cigar = compact_cigar_string(cigar_str)
    cigar_operations = re.findall(r'\d+[MIDNSHP=X]', cigar)

    if isinstance(cigar_operations, list) and len(cigar_operations) > 0 and isinstance(cigar_operations[0], list) and len(cigar_operations[0]) > 0:
        # The list is multidimensional, and the expression is valid.
        if cigar_operations[0][-1] == 'D':
            cigar_operations.pop(0)


    if isinstance(cigar_operations, list) and len(cigar_operations) > 0 and isinstance(cigar_operations[0], list) and len(cigar_operations[0]) > 0:
        # The list is multidimensional, and the expression is valid.
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
    return alignment