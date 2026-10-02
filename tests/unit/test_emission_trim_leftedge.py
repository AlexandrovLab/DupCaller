"""N1: read emission in calculateSSPosterior is conditioned on the read
showing one of the two modelled alleles -- P(other allele) =
(e/3)/(1 - 2e/3), P(true allele) = 1 minus that -- with the Pseq == 0
sentinel exactly uninformative.
N6: trim_mask with right == 0 trims nothing (a bare mask[-0:] masked the
whole family).
N5: an indel anchored at the window's first base is uncallable in both
learning and calling (left-alignment can't see past the window start)."""

import numpy as np
import pysam
import pytest

from DupCaller_sub.funcs.indels import indel_has_context
from DupCaller_sub.funcs.learn import profileTriNucMismatches
from DupCaller_sub.funcs.misc import (
    fallback_error_file,
    load_error_matrices,
    trim_mask,
)
from DupCaller_sub.funcs.prob import calculateSSPosterior, genotypeDSIndel


# ------------------------------------------------------------------ N1


def _p_other(q):
    e = 10 ** (-q / 10)
    return (e / 3) / (1 - 2 * e / 3)


@pytest.mark.parametrize("q", [11, 18, 25, 37])
def test_emission_is_two_allele_conditional(q):
    pseq = np.full([1, 2], -q / 10 * np.log(10))
    shows_tested = np.array([[True, False]])
    zero = np.zeros(2)
    h1, h2 = calculateSSPosterior(zero, zero, shows_tested, pseq)
    p_other = _p_other(q)
    # No amp error: under H1 (strand carries the tested allele) a read
    # showing it has probability 1 - p_other, one showing the other allele
    # p_other; H2 mirrors that.
    assert h1 == pytest.approx(np.log([1 - p_other, p_other]))
    assert h2 == pytest.approx(np.log([p_other, 1 - p_other]))


def test_emission_q18_value():
    # e = 10^-1.8 = 0.015849; e/3 = 0.0052831; 1 - 2e/3 = 0.98943
    assert _p_other(18) == pytest.approx(0.0053394, rel=1e-4)


def test_emission_with_amp_error_mixes_both_outcomes():
    q, P, P_rev = 25, 1e-3, 2e-3
    pseq = np.full([1, 1], -q / 10 * np.log(10))
    h1, h2 = calculateSSPosterior(
        np.array([P]), np.array([P_rev]), np.array([[True]]), pseq
    )
    p_other = _p_other(q)
    assert h1[0] == pytest.approx(np.log((1 - P) * (1 - p_other) + P * p_other))
    assert h2[0] == pytest.approx(np.log((1 - P_rev) * p_other + P_rev * (1 - p_other)))


def test_emission_sentinel_is_uninformative_and_input_untouched():
    pseq = np.zeros([2, 1])
    before = pseq.copy()
    zero = np.zeros(1)
    h1, h2 = calculateSSPosterior(zero, zero, np.array([[True], [False]]), pseq)
    assert h1[0] == pytest.approx(2 * np.log(0.5))
    assert h2[0] == pytest.approx(2 * np.log(0.5))
    np.testing.assert_array_equal(pseq, before)


# ------------------------------------------------------------------ N6


def test_trim_mask_right_zero_trims_nothing():
    assert not trim_mask(10, 0, 0).any()
    np.testing.assert_array_equal(np.nonzero(trim_mask(10, 2, 0))[0], [0, 1])


def test_trim_mask_left_and_right():
    np.testing.assert_array_equal(np.nonzero(trim_mask(10, 0, 1))[0], [9])
    np.testing.assert_array_equal(np.nonzero(trim_mask(10, 2, 3))[0], [0, 1, 7, 8, 9])


# ------------------------------------------------------------------ N5

HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
)
WIN = 100


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


def _read(cigar, seq, read1, name="r"):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    rec.query_sequence = seq
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    rec.reference_id = 0
    rec.reference_start = WIN
    rec.mapping_quality = 60
    rec.cigar = cigar
    rec.flag = 1 | (64 if read1 else 128)
    return rec


def _deletion_family(ref, del_at):
    lead = del_at - WIN
    cigar = [(0, lead), (2, 1), (0, len(ref) - lead - 1)]
    seq = ref[:lead] + ref[lead + 1 :]
    return [_read(cigar, seq, True, f"t{i}") for i in range(3)] + [
        _read(cigar, seq, False, f"b{i}") for i in range(3)
    ]


def _genotype(params, ref, family):
    n = len(ref)
    ref_int = np.array(["ATCG".index(b) for b in ref])
    hp = np.zeros([3, n])
    hp[0] = 1
    return genotypeDSIndel(
        family, WIN, WIN + n, ref_int, np.ones(n, bool), hp, np.zeros([3, n]), params
    )


def test_has_context_excludes_window_first_base():
    assert not indel_has_context(0, 10)
    assert indel_has_context(1, 10)
    assert not indel_has_context(9, 10)


def test_indel_anchored_at_window_start_dropped(params):
    ref = "ACGTTGCAAGCTTACG"
    # Deleting ref[1] ('C') is anchored at the window's first base.
    assert len(_genotype(params, ref, _deletion_family(ref, WIN + 1))[2]) == 0
    # Deleting ref[2] ('G') is anchored one base in: scored.
    assert _genotype(params, ref, _deletion_family(ref, WIN + 2))[2] == [
        f"{WIN + 1}:-1"
    ]


def test_truncated_left_alignment_in_leading_run_dropped(params):
    # The family starts inside an A run: a deletion of any A left-aligns
    # into the window's first base, where the run's true start is unknown.
    ref = "AAAAGCTTACGGTCAG"
    assert len(_genotype(params, ref, _deletion_family(ref, WIN + 3))[2]) == 0


def test_learn_skips_indel_anchored_at_window_start():
    ref = "ACGTTGCAAGCTTACG"
    n = len(ref)
    ref_int = np.array(["ATCG".index(b) for b in ref])
    cut = np.ones(n, dtype=bool)
    cut[1:] = ref_int[1:] != ref_int[:-1]
    run_id = np.cumsum(cut) - 1
    hp_raw = np.vstack((np.bincount(run_id)[run_id], cut)).astype(float)

    def hp_events(del_at):
        # One deletion read per strand among clean reads: an amp event.
        lead = del_at - WIN
        dcig = [(0, lead), (2, 1), (0, n - lead - 1)]
        dseq = ref[:lead] + ref[lead + 1 :]
        fam = []
        for read1 in (True, False):
            fam += [_read([(0, n)], ref, read1, f"{read1}{i}") for i in range(3)]
            fam.append(_read(dcig, dseq, read1, f"{read1}d"))
        res = profileTriNucMismatches(
            seqs=fam,
            reference_start=WIN,
            reference_int=ref_int,
            trinuc_int=np.zeros(n, dtype=int),
            hp_raw=hp_raw,
            str_raw=np.zeros([3, n]),
            antimask=np.ones(n, dtype=bool),
            params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
        )
        return res[1][:, [0, 3, 6, 9]].sum()

    assert hp_events(WIN + 1) == 0
    assert hp_events(WIN + 2) > 0
