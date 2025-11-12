import src.nw as align

import pytest


@pytest.mark.parametrize('sequences,expected_score,expected_aligned_sequences',
                         [
                             # identical sequences
                             (('ACGT', 'ACGT'), 20, ('ACGT', 'ACGT')),
                             # empty and non-empty sequences
                             (('ACGT', ''), -40, ('ACGT', '----')),
                             # two empty sequences
                             (('', ''), 0, ('', '')),
                             # first letter mismatch
                             (('AACGT', 'ACGT'), 10, ('AACGT', '-ACGT')),
                             # middle letter mismatch
                             (('AATCG', 'AACG'), 10, ('AATCG', 'AA-CG')),
                             # last letter mismatch
                             (('AATCG', 'AATC'), 10, ('AATCG', 'AATC-')),
                             # adjacent letters swap
                             (('AATCG', 'AACTG'), 7, ('AATCG', 'AACTG')),
                             # cutted prefix
                             (('AATCG', 'CG'), -20, ('AATCG', '---CG')),
                         ])
def test_nw(sequences, expected_score, expected_aligned_sequences):
    # match=5, mismatch=-4, gap=-10
    seq1, seq2 = sequences
    expected_seq1, expected_seq2 = expected_aligned_sequences
    score, aligned_seq1, aligned_seq2 = align.needleman_wunsch(seq1,
                                                 seq2, 
                                                 score_fun=lambda x, y: 5 if x == y else -4, 
                                                 gap_penalty=-10)
    assert score == expected_score
    assert aligned_seq1 == expected_seq1
    assert aligned_seq2 == expected_seq2


def test_nw_custom_score_fun():
    # match=1, mismatch=-1, gap=-5
    seq1 = 'ACGT'
    seq2 = 'ACG'
    score, aligned_seq1, aligned_seq2 = align.needleman_wunsch(seq1, seq2,
                                                               score_fun=lambda x, y: align.score_fun(x, y, match_score=1, mismatch_score=-1),
                                                               gap_penalty=-5)

    assert score == -2
    assert aligned_seq1 == 'ACGT'
    assert aligned_seq2 == 'ACG-'
