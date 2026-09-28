"""Indel coordinate consistency between learning and calling.

- left_align_indel must never shift an anchor left of reference_start
  (it used to return reference_start - 1 at the window's left edge).
- profileTriNucMismatches (learn) and genotypeDSIndel (call) must apply
  the same locus mask (learn used to test the base BEFORE the anchor) and
  classify the same STR context (call used to look at the base AFTER a
  deletion, learn at the first deleted base).
"""
import numpy as np
import pysam
import pytest

from DupCaller_sub.funcs import prob
from DupCaller_sub.funcs.indels import (
    indel_context_index,
    indel_mask_span,
    indel_passes_mask,
    left_align_indel,
)
from DupCaller_sub.funcs.learn import profileTriNucMismatches

CHROM = "chr1"
WIN_START = 100
_B2N = {"A": 0, "T": 1, "C": 2, "G": 3}


def _ref_int(seq):
    return np.array([_B2N[b] for b in seq])


# ---------------------------------------------------------- left_align_indel


@pytest.mark.parametrize(
    "indel, expected",
    [
        ("100:-1", "100:-1"),
        ("101:-1", "100:-1"),
        ("100:1:A", "100:1:A"),
        ("101:1:A", "100:1:A"),
    ],
)
def test_left_align_stops_at_window_start_homopolymer(indel, expected):
    ref = _ref_int("A" * 20)
    assert left_align_indel(indel, ref, WIN_START) == expected


@pytest.mark.parametrize(
    "indel, expected",
    [
        ("100:-2", "100:-2"),
        ("103:-2", "100:-2"),
        ("100:2:TA", "100:2:TA"),
        ("103:2:AT", "100:2:TA"),
    ],
)
def test_left_align_stops_at_window_start_dinucleotide(indel, expected):
    ref = _ref_int("ATATATATATATGC")
    assert left_align_indel(indel, ref, WIN_START) == expected


def test_left_align_interior_unchanged():
    # Room to the left: normal canonicalization still reaches the run start.
    ref = _ref_int("GC" + "A" * 8 + "GC")
    assert left_align_indel("108:-1", ref, WIN_START) == "101:-1"
    assert left_align_indel("108:1:A", ref, WIN_START) == "101:1:A"


# ------------------------------------------------------------ mask helpers


def test_mask_span_insertion_and_deletions():
    assert indel_mask_span(5, 1) == (5, 6)
    assert indel_mask_span(5, 3) == (5, 6)
    assert indel_mask_span(5, -1) == (5, 7)
    assert indel_mask_span(5, -5) == (5, 11)


def test_passes_mask_first_position_and_edges():
    antimask = np.ones(10, dtype=bool)
    assert indel_passes_mask(antimask, 0, 1)
    assert indel_passes_mask(antimask, 0, -1)
    assert indel_passes_mask(antimask, 0, -3)
    # Outside the window fails instead of wrapping/truncating.
    assert not indel_passes_mask(antimask, -1, 1)
    assert not indel_passes_mask(antimask, 8, -2)


def test_passes_mask_uses_anchor_not_preceding_base():
    antimask = np.ones(10, dtype=bool)
    antimask[3] = False
    assert indel_passes_mask(antimask, 4, 1)  # preceding base masked only
    assert not indel_passes_mask(antimask, 3, 1)  # anchor masked
    assert not indel_passes_mask(antimask, 2, -1)  # deleted base masked


def test_context_index_is_first_base_after_anchor():
    assert indel_context_index(7) == 8


# ------------------------------------------- learn vs call on a real family

# Window 100..129: G C [A T A T] G C G G ... ; (AT)x2 annotated as an STR
# at 102..105. A 4bp deletion of ATAT is anchored at 101 ("101:-4"): its
# first deleted base (102) is inside the STR, the base after it (106) is
# not -- the case where the old call-side anchor disagreed with learn.
REF = "GC" + "ATAT" + "GCGGCCGCGGCCGCGGCCGCGGCC"
STR_START, STR_END = 102, 106
DEL_POS, DEL_LEN = 101, 4
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
)


def _hp_raw(seq):
    ref = _ref_int(seq)
    cut = np.ones(len(ref), dtype=bool)
    cut[1:] = ref[1:] != ref[:-1]
    run_id = np.cumsum(cut) - 1
    return np.vstack((np.bincount(run_id)[run_id], cut)).astype(float)


def _str_raw():
    s = np.zeros([3, len(REF)])
    s[0, STR_START - WIN_START : STR_END - WIN_START] = 2
    s[1, STR_START - WIN_START : STR_END - WIN_START] = 2
    s[2, STR_START - WIN_START] = 1
    return s


def _read(name, is_read1, with_del, ref=REF, del_pos=DEL_POS, del_len=DEL_LEN):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    if with_del:
        lead = del_pos + 1 - WIN_START
        rec.query_sequence = ref[:lead] + ref[lead + del_len :]
        rec.cigar = [(0, lead), (2, del_len), (0, len(ref) - lead - del_len)]
    else:
        rec.query_sequence = ref
        rec.cigar = [(0, len(ref))]
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(rec.query_sequence))
    rec.reference_id = 0
    rec.reference_start = WIN_START
    rec.mapping_quality = 60
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False
    return rec


def _family(**kw):
    # F1R2 all carry the deletion, F2R1 all reference: a single-strand
    # (damage-type) event, so learn credits the dmg matrices.
    return [_read(f"t{i}", True, True, **kw) for i in range(3)] + [
        _read(f"b{i}", False, False, **kw) for i in range(3)
    ]


def _learn_event_rows(antimask):
    result = profileTriNucMismatches(
        seqs=_family(),
        reference_start=WIN_START,
        reference_int=_ref_int(REF),
        trinuc_int=np.zeros(len(REF), dtype=int),
        hp_raw=_hp_raw(REF),
        str_raw=_str_raw(),
        antimask=antimask,
        params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
    )
    str_dmg_count = result[5]
    col = -DEL_LEN + 5
    return set(np.nonzero(str_dmg_count[:, col])[0].tolist())


def _call_contexts(monkeypatch, family, ref, antimask, str_raw):
    seen = []

    def fake_probs(hp, str_bin, *args):
        seen.append((int(hp), int(str_bin)))
        return (1e-3, 1 - 1e-3, 1e-3, 1 - 1e-3, 1e-3, 1 - 1e-3)

    monkeypatch.setattr(prob, "indelErrorProbs", fake_probs)
    prob.genotypeDSIndel(
        family,
        WIN_START,
        WIN_START + len(ref),
        _ref_int(ref),
        antimask,
        _hp_raw(ref),
        str_raw,
        {
            "ampmat_hp": None,
            "dmgmat_hp": None,
            "ampmat_str": None,
            "dmgmat_str": None,
            "minBq": 10,
        },
    )
    return seen


def _call_event_rows(antimask, monkeypatch):
    seen = _call_contexts(monkeypatch, _family(), REF, antimask, _str_raw())
    return {str_bin for _, str_bin in seen}


def _antimask(masked=()):
    a = np.ones(len(REF), dtype=bool)
    for p in masked:
        a[p - WIN_START] = False
    return a


@pytest.mark.parametrize(
    "masked, admitted",
    [
        ((), True),
        ((DEL_POS - 1,), True),  # only the preceding base is masked
        ((DEL_POS,), False),  # anchor masked
        ((DEL_POS + 2,), False),  # a deleted base masked
    ],
)
def test_learn_and_call_agree_on_mask_and_str_bin(masked, admitted, monkeypatch):
    learn_rows = _learn_event_rows(_antimask(masked))
    call_rows = _call_event_rows(_antimask(masked), monkeypatch)
    assert learn_rows == call_rows
    # (AT)x2 -> total length 4 -> STR bin 1.
    assert learn_rows == ({1} if admitted else set())


def test_learn_and_call_agree_on_hp_length(monkeypatch):
    # G A A A A A [C] C G T: deleting the C at 106 is anchored at the last A
    # (105). Its own run is CC (length 2); the old call-side max over
    # pos..pos+2 picked up the neighbouring A-run (length 5) instead.
    ref = "G" + "A" * 5 + "CC" + "GT" * 11
    kw = {"ref": ref, "del_pos": 105, "del_len": 1}
    family = _family(**kw)
    antimask = np.ones(len(ref), dtype=bool)
    str_raw = np.zeros([3, len(ref)])

    result = profileTriNucMismatches(
        seqs=family,
        reference_start=WIN_START,
        reference_int=_ref_int(ref),
        trinuc_int=np.zeros(len(ref), dtype=int),
        hp_raw=_hp_raw(ref),
        str_raw=str_raw,
        antimask=antimask,
        params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
    )
    hp_dmg_count = result[4]
    c = _B2N["C"]
    learn_hp = {int(r) + 1 for r in np.nonzero(hp_dmg_count[:, c * 3 + 0] > 0)[0]}

    call_hp = {
        hp for hp, _ in _call_contexts(monkeypatch, family, ref, antimask, str_raw)
    }
    assert learn_hp == call_hp == {2}
