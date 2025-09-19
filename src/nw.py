from typing import Callable, Tuple
import argparse
import sys
from abc import ABC, abstractmethod


INF = sys.maxsize//3

class MatrixStorage(ABC):
    @abstractmethod
    def get(self, i: int, j: int) -> int:
        pass

    @abstractmethod
    def set(self, i: int, j: int, value: int) -> None:
        pass

    @abstractmethod
    def is_valid(self, i: int, j: int) -> bool:
        pass

    @abstractmethod
    def get_k(self) -> int:
        pass


class FullMatrix(MatrixStorage):
    def __init__(self, n: int, k: int) -> None:
        self.matrix = [[0 for _ in range(n)] for _ in range(n)]

    def get(self, i: int, j: int) -> int:
        if not self.is_valid(i, j):
            return -INF
        return self.matrix[i][j]

    def set(self, i: int, j: int, value: int) -> None:
        if not self.is_valid(i, j):
            return
        self.matrix[i][j] = value

    def is_valid(self, i: int, j: int) -> bool:
        return i >= 0 and i < len(self.matrix) and j >= 0 and j < len(self.matrix[0])

    def get_k(self) -> int:
        return len(self.matrix)


class BandedMatrix(MatrixStorage):
    def __init__(self, n: int, k: int) -> None:
        self.k = k
        self.matrix = [[0 for _ in range(2*k+1)] for _ in range(n)]

    def translate_coords(self, i: int, j: int) -> Tuple[int, int]:
        return j, i-j

    def get(self, i: int, j: int) -> int:
        if not self.is_valid(i, j):
            return -INF
        if abs(i-j) > self.k:
            return -INF
        i, j = self.translate_coords(i, j)
        return self.matrix[i][j]

    def set(self, i: int, j: int, value: int) -> None:
        if not self.is_valid(i, j):
            return
        if abs(i-j) > self.k:
            return
        i, j = self.translate_coords(i, j)
        self.matrix[i][j] = value

    def is_valid(self, i: int, j: int) -> bool:
        return i >= 0 and i < len(self.matrix) and j >= 0 and j < len(self.matrix)

    def get_k(self) -> int:
        return self.k


def score_fun(a: str, 
              b: str,
              match_score: int = 5, 
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score

def needleman_wunsch(seq1: str,
                     seq2: str,
                     matrix: MatrixStorage,
                     score_fun: Callable[[str, str], int] = score_fun,
                     gap_penalty: int = -10) -> int:

    """Given two sequences, aligns them using the Needleman-Wunsch algorithm.

    This function takes two sequences and optionally a scoring function and a
    gap penalty value as arguments. 
    The function returns the optimal alignment score.

    Args:
        seq1: The first sequence, e.g. 'ACCGT'
        seq2: The second sequence, e.g. 'ACGT'
        matrix: Subclass of MatrixStorage for storing score
        score_fun: The scoring function, e.g. score_fun('A', 'A') returns 5
        gap_penalty: The gap penalty value, e.g. -10

    Returns:
        score: The optimal alignment score, e.g. 10
    """
    # Initialize the score matrix.
    ncol = len(seq2)+1
    nrow = len(seq1)+1

    # Fill base elements.
    for i in range(nrow):
        matrix.set(i, 0, i*gap_penalty)
    for j in range(ncol):
        matrix.set(0, j, j*gap_penalty)

    k = matrix.get_k()

    # Fill score and backtracking matrices.
    for i in range(1, nrow):
        for j in range(max(1, i-k), min(ncol, i+k+1)):
            matrix.set(i, j, max(
                matrix.get(i-1, j)+gap_penalty,
                matrix.get(i, j-1)+gap_penalty,
                matrix.get(i-1, j-1)+score_fun(seq1[i-1], seq2[j-1]),
            ))

    return matrix.get(nrow-1, ncol-1)

def print_results(score: int, file = None):
    """Prints the results of the Needleman-Wunsch algorithm.

    This function takes two aligned sequences and the optimal alignment score
    as arguments. It prints the sequences and the score to the standard output
    or to a file.

    Args:
        score: The optimal alignment score, e.g. 10
        file: The file to print to. If None, prints to the standard output.

    Returns:
        None
    """
    if file is None:
        file = sys.stdout

    print("Score: %s" % score, file=file)

def main():
    parser = argparse.ArgumentParser(description='Needleman-Wunsch algorithm')
    parser.add_argument('seq1', help='first sequence')
    parser.add_argument('seq2', help='second sequence')
    parser.add_argument('--match', type=int, help='match score')
    parser.add_argument('--mismatch', type=int, help='mismatch score')
    parser.add_argument('--gap', type=int, default=-10, help='gap penalty')
    parser.add_argument('--debug', action='store_true', help='debug mode')
    args = parser.parse_args()

    global DEBUG
    DEBUG = args.debug

    if args.match and args.mismatch:
        score = needleman_wunsch(args.seq1,
                                             args.seq2,
                                             FullMatrix(len(args.seq1)+1, len(args.seq2)+1),
                                             score_fun=lambda x, y: args.match if x == y else args.mismatch,
                                             gap_penalty=args.gap)
    else:
        assert not args.match and not args.mismatch, "match and mismatch must be specified together"
        score = needleman_wunsch(args.seq1,
                                             args.seq2,
                                             FullMatrix(len(args.seq1)+1, len(args.seq2)+1),
                                             score_fun=score_fun,
                                             gap_penalty=args.gap)
    print_results(score)

    return score

if __name__ == '__main__':
    main()
