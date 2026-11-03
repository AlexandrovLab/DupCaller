"""call --mate-merge helpers: merged-family detection, per-read 3' trim, collapsing the two mates of a read pair to one
observation (SNVs and indels), and the fragment-end trim of a merged family's window."""
from types import SimpleNamespace

import numpy as np

from DupCaller_sub.funcs.call import _merged_trim
from DupCaller_sub.funcs.indels import INDEL_ALT, INDEL_REF, INDEL_UNINFORMATIVE
from DupCaller_sub.funcs.prob import (
    _blank_read_3prime,
    _collapse_mates,
    _collapse_mates_indel,
    is_merged_family,
)

N = 4  # no base


def _read(name, start, end, forward, tlen):
    return SimpleNamespace(query_name=name, reference_start=start, reference_end=end, is_forward=forward,
                           is_reverse=not forward, template_length=tlen)


def test_is_merged_family():
    fwd = _read("p1", 100, 150, True, 80)
    rev = _read("p1", 130, 180, False, -80)
    assert is_merged_family([fwd, rev])
    assert not is_merged_family([fwd, _read("p2", 100, 150, True, 80)])


def test_blank_read_3prime_trims_each_read_at_its_own_3_end():
    reads = [_read("p1", 100, 110, True, 20), _read("p1", 104, 114, False, -20)]
    seq = np.zeros([2, 14], dtype=int)
    qual = np.full([2, 14], 30)
    _blank_read_3prime(seq, qual, reads, 100, trim3=3)
    assert (seq[0, 7:10] == N).all() and (seq[0, :7] == 0).all()  # forward read: last 3 bases of 100-110
    assert (seq[1, 4:7] == N).all() and (seq[1, 7:] == 0).all()  # reverse read: first 3 bases of 104-114
    assert (qual[seq == N] == 0).all()


def test_blank_read_3prime_zero_trim_is_a_no_op():
    reads = [_read("p1", 100, 110, True, 20)]
    seq, qual = np.zeros([1, 10], dtype=int), np.full([1, 10], 30)
    _blank_read_3prime(seq, qual, reads, 100, trim3=0)
    assert (seq == 0).all() and (qual == 30).all()


def test_collapse_mates_one_row_per_pair():
    reads = [_read("p1", 0, 5, True, 5), _read("p1", 0, 5, False, -5), _read("p2", 0, 5, True, 5)]
    #                 agree  disagree  disagree,  only     only
    #                         (both ok) mate1 low  mate1    mate2
    seq = np.array([[0,     1,        2,         3,       N],
                    [0,     2,        3,         N,       1],
                    [1,     1,        1,         1,       1]])
    qual = np.array([[30,   30,       10,        30,      0],
                     [40,   30,       30,        0,       30],
                     [30,   30,       30,        30,      30]])
    out_seq, out_qual, out_reads = _collapse_mates(seq, qual, reads, min_bq=20)
    assert [r.query_name for r in out_reads] == ["p1", "p2"]
    assert out_seq[0].tolist() == [0, N, 3, 3, 1]
    assert out_qual[0].tolist() == [40, 0, 30, 30, 30]
    assert out_seq[1].tolist() == seq[2].tolist() and out_qual[1].tolist() == qual[2].tolist()


def test_collapse_mates_unpaired_family_unchanged():
    reads = [_read("p1", 0, 3, True, 3), _read("p2", 0, 3, True, 3)]
    seq, qual = np.zeros([2, 3], dtype=int), np.full([2, 3], 30)
    out = _collapse_mates(seq, qual, reads, min_bq=20)
    assert out[0] is seq and out[1] is qual and out[2] is reads


def test_collapse_mates_indel():
    reads = [_read("p1", 0, 5, True, 5), _read("p1", 0, 5, False, -5)]
    U = INDEL_UNINFORMATIVE
    #                   agree      disagree   one informative  neither
    states = np.array([[INDEL_ALT, INDEL_ALT, U,               U],
                       [INDEL_ALT, INDEL_REF, INDEL_REF,       U]])
    quals = np.array([[20, 30, 0, 0],
                      [35, 30, 25, 0]])
    out_s, out_q = _collapse_mates_indel(states, quals, reads, INDEL_ALT, INDEL_REF, U)
    assert out_s.shape == (1, 4)
    assert out_s[0].tolist() == [INDEL_ALT, U, INDEL_REF, U]
    assert out_q[0].tolist() == [35, 0, 25, 0]


def test_merged_trim_masks_read_through_tails_and_fragment_ends():
    # Read-through pair: fragment 100-140, reads run 10 bp past the other mate's 5' end.
    fwd = _read("p1", 100, 150, True, 40)
    rev = _read("p1", 90, 140, False, -40)
    left, right = _merged_trim([fwd, rev], 90, 150, trim5=8)
    assert (left, right) == (10 + 8, 10 + 8)


def test_merged_trim_overlapping_pair_without_read_through():
    fwd = _read("p1", 100, 150, True, 80)
    rev = _read("p1", 130, 180, False, -80)
    assert _merged_trim([fwd, rev], 100, 180, trim5=8) == (8, 8)
