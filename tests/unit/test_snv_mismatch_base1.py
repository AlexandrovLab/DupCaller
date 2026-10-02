"""genotypeDSSnv: base1 is the mismatch (non-reference) base, base2 the
reference; positions with more than one non-reference base are masked."""

import itertools
from types import SimpleNamespace

import numpy as np
import pytest

from DupCaller_sub.funcs.misc import (
    build_trinuc64_order,
    fallback_error_file,
    load_error_matrices,
)
from DupCaller_sub.funcs.prob import genotypeDSSnv, strand_pvalue

BASE2NUM = {"A": 0, "T": 1, "C": 2, "G": 3}
TRINUC2NUM, _ = build_trinuc64_order()


@pytest.fixture(scope="module")
def params():
    pre = fallback_error_file(".amp.tn.srd.txt")[: -len(".amp.tn.srd.txt")]
    p = {
        "amperr_file": pre + ".amp.tn.srd.txt",
        "dmgerr_file": pre + ".dmg.tn.txt",
        "amperri_file": pre + ".amp.id.txt",
        "dmgerri_file": pre + ".dmg.id.txt",
        "pseudocount": 0.5,
        "minBq": 10,
    }
    load_error_matrices(p)
    return p


def _read(seq, top):
    # top strand = F1R2 (read1 forward), bottom = F2R1 (read2 forward)
    return SimpleNamespace(
        is_read1=top,
        is_read2=not top,
        is_forward=True,
        is_reverse=False,
        query_alignment_qualities=[37] * len(seq),
        query_alignment_sequence=seq,
        cigartuples=[(0, len(seq))],
        reference_start=0,
        reference_length=len(seq),
    )


def _genotype(params, ref_trinuc, top_bases, bot_bases):
    """One family over the 3bp ref_trinuc; only the middle base varies."""
    reads = [_read(ref_trinuc[0] + b + ref_trinuc[2], True) for b in top_bases]
    reads += [_read(ref_trinuc[0] + b + ref_trinuc[2], False) for b in bot_bases]
    ref = np.array([BASE2NUM[b] for b in ref_trinuc])
    trinuc = np.array([64, TRINUC2NUM[ref_trinuc], 64])
    out = genotypeDSSnv(
        reads, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    _, LR, _, mut_antimask, base1, antimask = out[:6]
    return LR, mut_antimask[1], base1[1], antimask[1]


def test_ref_c_with_only_a_and_t_is_masked(params):
    LR, is_mut, b1, anti = _genotype(params, "ACA", "AAAA", "TTTT")
    assert not anti and not is_mut and LR.size == 0


def test_ref_a_with_only_c_and_g_is_masked(params):
    _, is_mut, _, anti = _genotype(params, "CAC", "CCCG", "CCCG")
    assert not anti and not is_mut


def test_ref_plus_two_alts_is_masked(params):
    _, is_mut, _, anti = _genotype(params, "ACA", "CCCA", "CCCT")
    assert not anti and not is_mut


def test_ref_plus_one_alt_is_a_candidate(params):
    LR, is_mut, b1, anti = _genotype(params, "ACA", "TTTT", "TTTT")
    assert anti and is_mut and b1 == BASE2NUM["T"] and LR[0] > 0


def test_minority_alt_is_base1_but_rejected_to_ref(params):
    # 1 alt read among refs: LR < 0, reported as reference coverage.
    LR, is_mut, b1, anti = _genotype(params, "ACA", "CCCT", "CCCC")
    assert anti and not is_mut and b1 == BASE2NUM["C"] and LR.size == 0


def test_negative_lr_site_returned_for_mu(params):
    # The same LR < 0 site is still handed back (position, mismatch base,
    # LR) for the per-channel mu solve.
    reads = [_read("A" + b + "A", True) for b in "CCCT"]
    reads += [_read("ACA", False) for _ in range(4)]
    ref = np.array([BASE2NUM[b] for b in "ACA"])
    trinuc = np.array([64, TRINUC2NUM["ACA"], 64])
    out = genotypeDSSnv(
        reads, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    neg_pos, neg_alt, neg_LR = out[-3:]
    assert neg_pos.tolist() == [1]
    assert neg_alt.tolist() == [BASE2NUM["T"]]
    assert neg_LR.size == 1 and neg_LR[0] < 0


def test_positive_candidate_not_in_negative_outputs(params):
    reads = [_read("ATA", True) for _ in range(4)]
    reads += [_read("ATA", False) for _ in range(4)]
    ref = np.array([BASE2NUM[b] for b in "ACA"])
    trinuc = np.array([64, TRINUC2NUM["ACA"], 64])
    out = genotypeDSSnv(
        reads, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    assert out[-3].size == 0 and out[-1].size == 0


def test_all_reference_is_reference(params):
    _, is_mut, b1, anti = _genotype(params, "ACA", "CCC", "CCC")
    assert anti and not is_mut and b1 == BASE2NUM["C"]


@pytest.mark.parametrize("alt_a, alt_b", list(itertools.permutations("ATG", 2)))
def test_two_alt_masking_is_base_order_invariant(params, alt_a, alt_b):
    _, is_mut, _, anti = _genotype(params, "ACA", alt_a * 3, alt_b * 3)
    assert not anti and not is_mut


def test_all_reference_family_returns_empty_negative_outputs(params):
    # No mismatch anywhere: the early-return path must return the same
    # number of values as the full path.
    reads = [_read("ACA", True) for _ in range(3)]
    reads += [_read("ACA", False) for _ in range(3)]
    ref = np.array([BASE2NUM[b] for b in "ACA"])
    trinuc = np.array([64, TRINUC2NUM["ACA"], 64])
    out = genotypeDSSnv(
        reads, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    assert len(out) == 13 and out[-3].size == 0 and out[-1].size == 0


def test_strand_test_uses_actual_bq_below_minbq(params):
    # C>T call; top strand T(BQ 37) plus two A reads at BQ 8 (below the
    # fixture's minBq 10). Genotyping ignores the A reads, but the strand
    # test counts them at BQ 8, so the top-strand p drops below 1.
    top = [_read("ATA", True), _read("AAA", True), _read("AAA", True)]
    for r in top[1:]:
        r.query_alignment_qualities = [37, 8, 37]
    bot = [_read("ATA", False) for _ in range(3)]
    ref = np.array([BASE2NUM[b] for b in "ACA"])
    trinuc = np.array([64, TRINUC2NUM["ACA"], 64])
    out = genotypeDSSnv(
        top + bot, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    LR, is_mut, ev_top, ev_bot = out[1], out[3][1], out[8], out[9]
    assert is_mut and LR[0] > 0
    quals, k, _ = ev_top[0]
    assert sorted(quals.tolist()) == [8, 8, 37] and k == 2
    q = 10**-0.8
    # P(both BQ-8 reads non-alt and the BQ-37 read alt) dominates the tail.
    assert 0.9 * q * q < strand_pvalue(*ev_top[0]) < 1.1 * (q * q + 3 * q * 2e-4)
    assert strand_pvalue(*ev_bot[0]) == 1.0
