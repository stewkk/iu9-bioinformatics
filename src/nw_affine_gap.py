from typing import Callable, Tuple

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

    #infinity = 2 * gap_open + (n + m - 2) * gap_extend + 1
    infinity = float('-inf')

    # 1. Initialize matrices


    # 2. Fill matrices
    # We assume that consecutive gaps on different sequences are not allowed
    for i in range(1, n):
        for j in range(1, m):
            pass


    # 3. Traceback
    aln1 = ''
    aln2 = ''
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        pass

    
    return aln1[::-1], aln2[::-1], score

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