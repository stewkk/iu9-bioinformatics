import src.nw as align
import pytest


@pytest.mark.parametrize('sequences,expected_score', [
    # identical sequences
    (('ACGT', 'ACGT'), 20),
    # completely different sequences
    (('AAAA', 'TTTT'), -16),
    # single mismatch
    (('ACGT', 'ACAT'), 11),
    # one match
    (('ACGT', 'TTGG'), -7),
    # all mismatches
    (('ACG', 'TGC'), -12),
    # empty sequences
    (('', ''), 0),
    # single character match
    (('A', 'A'), 5),
    # single character mismatch
    (('A', 'T'), -4),
])
def test_nw_full_matrix(sequences, expected_score):
    """FullMatrix (match=5, mismatch=-4, gap=-10)"""
    seq1, seq2 = sequences
    score = align.needleman_wunsch(
        seq1, seq2,
        align.FullMatrix(len(seq1)+1, len(seq2)+1),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    assert score == expected_score


@pytest.mark.parametrize('sequences', [
    ('ACGT', 'ACGT'),
    ('AAAA', 'TTTT'),
    ('ACGT', 'ACAT'),
    ('ACGT', 'TTGG'),
    ('ACG', 'TGC'),
    ('A', 'A'),
    ('A', 'T'),
    ('ATCG', 'ATCG'),
    ('GATTACA', 'GATTACA'),
    ('ABCDEFG', 'ABCDEFG'),
])
def test_full_vs_banded(sequences):
    seq1, seq2 = sequences
    n = max(len(seq1), len(seq2))
    k = n
    
    score_full = align.needleman_wunsch(
        seq1, seq2,
        align.FullMatrix(len(seq1)+1, len(seq2)+1),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    
    score_banded = align.needleman_wunsch(
        seq1, seq2,
        align.BandedMatrix(len(seq1)+1, k),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    
    assert score_full == score_banded


@pytest.mark.parametrize('sequences,k', [
    (('ACGT', 'ACGT'), 1),
    (('ACGT', 'ACGT'), 2),
    (('ACGT', 'ACAT'), 1),
    (('ACGT', 'ACAT'), 2),
    (('AAAA', 'TTTT'), 1),
    (('AAAA', 'TTTT'), 2),
    (('ATCGATCG', 'ATCGATCG'), 2),
    (('ATCGATCG', 'ATCGATCG'), 3),
])
def test_k_vs_k_plus_1(sequences, k):
    seq1, seq2 = sequences
    
    score_k = align.needleman_wunsch(
        seq1, seq2,
        align.BandedMatrix(len(seq1)+1, k),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    
    score_k_plus_1 = align.needleman_wunsch(
        seq1, seq2,
        align.BandedMatrix(len(seq1)+1, k+1),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    
    assert score_k_plus_1 >= score_k


@pytest.mark.parametrize('sequences,k,expected_score,should_pass', [
    (('ACGT', 'ACGT'), 2, 20, True),
    (('ACGT', 'ACAT'), 2, 11, True),
    (('A', 'A'), 1, 5, True),
    (('AAAA', 'TTTT'), 1, -16, True),
    
    (('ACGT', 'TTGG'), 1, 2, False),
    (('', ''), 0, 0, True),
])
def test_banded(sequences, k, expected_score, should_pass):
    seq1, seq2 = sequences
    
    score = align.needleman_wunsch(
        seq1, seq2,
        align.BandedMatrix(len(seq1)+1, k),
        score_fun=lambda x, y: 5 if x == y else -4,
        gap_penalty=-10
    )
    
    if should_pass:
        assert score == expected_score
    else:
        assert score < expected_score
