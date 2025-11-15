from typing import Callable, Tuple
import copy

DEBUG = False

def score_fun(a: str, 
              b: str,
              match_score: int = 5, 
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score

def needleman_wunsch_affine(seq1: str, 
                            seq2: str, 
                            score_fun: Callable = score_fun, 
                            gap_open: int = -10, 
                            gap_extend: int = -1) -> Tuple[str, str, int]:
    '''
    Inputs:
    seq1 - first sequence
    seq2 - second sequence
    score_fun - function that takes two characters and returns score
    gap_open - gap open penalty
    gap_extend - gap extend penalty
    Outputs:
    aln1 - first aligned sequence
    aln2 - second aligned sequence
    score - score of the alignment
    '''

    n = len(seq1)+1
    m = len(seq2)+1
    infinity = float('-inf')

    # 1. Initialize matrices

    score_match = [[0 for _ in range(m)] for _ in range(n)]
    score_insertion = copy.deepcopy(score_match)
    score_deletion = copy.deepcopy(score_match)
    result = copy.deepcopy(score_match)

    score_match[0][0] = 0
    result[0][0] = 0
    for i in range(1, n):
        score_match[i][0] = infinity
        score_insertion[i][0] = infinity
        score_deletion[i][0] = gap_open + (i-1)*gap_extend
        result[i][0] = score_deletion[i][0]
    for j in range(1, m):
        score_match[0][j] = infinity
        score_insertion[0][j] = gap_open + (j-1)*gap_extend
        score_deletion[0][j] = infinity
        result[0][j] = score_insertion[0][j]

    # 2. Fill matrices
    for i in range(1, n):
        for j in range(1, m):
            score_match[i][j] = max(
                score_match[i-1][j-1] + score_fun(seq1[i-1], seq2[j-1]),
                score_insertion[i-1][j-1] + score_fun(seq1[i-1], seq2[j-1]),
                score_deletion[i-1][j-1] + score_fun(seq1[i-1], seq2[j-1]),
            )
            score_insertion[i][j] = max(
                score_insertion[i][j-1] + gap_extend,
                score_match[i][j-1] + gap_open,
                score_deletion[i][j-1] + gap_extend,
            )
            score_deletion[i][j] = max(
                score_deletion[i-1][j] + gap_extend,
                score_match[i-1][j] + gap_open,
                score_insertion[i-1][j] + gap_extend,
            )
            result[i][j] = max(
                score_deletion[i][j],
                score_insertion[i][j],
                score_match[i][j],
            )


    # 3. Traceback
    aln1 = list()
    aln2 = list()
    i = n-1
    j = m-1
    while i > 0 or j > 0:
        current_score = result[i][j]
        if i > 0 and j > 0 and score_match[i][j] == current_score:
            aln1.append(seq1[i-1])
            aln2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif j > 0 and score_insertion[i][j] == current_score:
            aln1.append('-')
            aln2.append(seq2[j-1])
            j -= 1
        else:
            aln1.append(seq1[i-1])
            aln2.append('-')
            i -= 1
            
    return ''.join(aln1[::-1]), ''.join(aln2[::-1]), result[n-1][m-1]

def print_array(matrix: list):
    for row in matrix:
        for element in row:
            print(f"{element:6}", end="")
        print()

def main():
    aln1, aln2, score = needleman_wunsch_affine("ACGT", "TAGT", gap_open=-10, gap_extend=-1) 
    print(f'str 1: {aln1}')
    print(f'str 2: {aln2}')
    print(f'score: {score}')
    


if __name__ == "__main__":
    main()
