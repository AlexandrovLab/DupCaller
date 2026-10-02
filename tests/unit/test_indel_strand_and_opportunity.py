"""Indel amp-error opportunity gating (learn) and the indel strand
independence test (call).

- profileTriNucMismatches credits indel amp opportunity only for reads
  whose base at the context position clears minBq -- the same reads that
  can be ALT/REF events (getIndelArr's BQ gate) -- and moves an ALT read
  out of the idLen=0 column only if it was credited there.
- indelReversionRate looks up the amp error that undoes the indel in the
  mutant molecule's context; indel_strand_pvalue is the binomial upper
  tail of REF reads under that rate.
"""
import numpy as np
import pysam
import pytest
from scipy.stats import binom

from DupCaller_sub.funcs.learn import profileTriNucMismatches
from DupCaller_sub.funcs.prob import (
    indel_strand_evidence,
    indel_strand_pvalue,
    indelReversionRate,
)

CHROM = "chr1"
WIN_START = 100
_B2N = {"A": 0, "T": 1, "C": 2, "G": 3}
# G C [A A A] G C G G ...: the only A/T homopolymer of length 3 starts at
# 102, so the HP3 A/T cells hold nothing but this run.
REF = "GC" + "AAA" + "GCGGCCGCGGCCGCGGCCGCGGCC"
RUN_START = 102
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
)


def _ref_int(seq):
    return np.array([_B2N[b] for b in seq])


def _hp_raw(seq):
    ref = _ref_int(seq)
    cut = np.ones(len(ref), dtype=bool)
    cut[1:] = ref[1:] != ref[:-1]
    run_id = np.cumsum(cut) - 1
    return np.vstack((np.bincount(run_id)[run_id], cut)).astype(float)


def _read(name, is_read1, del_at=None, low_bq_at=None):
    """Full-window read; del_at deletes that one reference base (CIGAR
    placement as given, not left-aligned); low_bq_at gives that
    reference position's base BQ 2."""
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    if del_at is None:
        seq = REF
        rec.cigar = [(0, len(REF))]
        ref_positions = list(range(WIN_START, WIN_START + len(REF)))
    else:
        lead = del_at - WIN_START
        seq = REF[:lead] + REF[lead + 1 :]
        rec.cigar = [(0, lead), (2, 1), (0, len(REF) - lead - 1)]
        ref_positions = [
            p for p in range(WIN_START, WIN_START + len(REF)) if p != del_at
        ]
    rec.query_sequence = seq
    quals = [40] * len(seq)
    if low_bq_at is not None:
        quals[ref_positions.index(low_bq_at)] = 2
    rec.query_qualities = pysam.qualitystring_to_array(
        "".join(chr(q + 33) for q in quals)
    )
    rec.reference_id = 0
    rec.reference_start = WIN_START
    rec.mapping_quality = 60
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False
    return rec


def _hp3_a_cells(alt_del_at):
    # F1R2: 3 REF reads (one with BQ 2 at the run start) + 1 read deleting
    # one A of the run. F2R1: 4 clean REF reads.
    family = [
        _read("t0", True),
        _read("t1", True),
        _read("t2", True, low_bq_at=RUN_START),
        _read("t3", True, del_at=alt_del_at),
    ] + [_read(f"b{i}", False) for i in range(4)]
    result = profileTriNucMismatches(
        seqs=family,
        reference_start=WIN_START,
        reference_int=_ref_int(REF),
        trinuc_int=np.zeros(len(REF), dtype=int),
        hp_raw=_hp_raw(REF),
        str_raw=np.zeros([3, len(REF)]),
        antimask=np.ones(len(REF), dtype=bool),
        params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
    )
    hp_alt_count = result[1]
    a = _B2N["A"]
    row = 3 - 1
    return hp_alt_count[row, a * 3 + 0], hp_alt_count[row, a * 3 + 1]


@pytest.mark.parametrize(
    "alt_del_at",
    [
        RUN_START,  # CIGAR deletes the run's first base (the context base)
        RUN_START + 2,  # CIGAR deletes the last A; left-aligns to the same event
    ],
)
def test_hp_opportunity_counts_only_bq_qualifying_reads(alt_del_at):
    n_del, n_ref = _hp3_a_cells(alt_del_at)
    # One deletion event either way.
    assert n_del == 1
    # idLen=0 column: F1R2's 2 high-BQ REF reads + F2R1's 4. The BQ-2 read
    # is not an opportunity (it can't be an event either), and the ALT read
    # is an event, not a reference observation, whichever run base its
    # CIGAR removed.
    assert n_ref == 6


# ------------------------------------------------------ indelReversionRate


def _mats():
    hp = np.zeros([10, 12])
    for row in range(10):
        for base in range(4):
            hp[row, base * 3 + 0] = (row + 1) * 1e-4 + base * 1e-6  # -1
            hp[row, base * 3 + 1] = 1.0
            hp[row, base * 3 + 2] = (row + 1) * 1e-5 + base * 1e-7  # +1
    st = np.zeros([5, 11])
    for row in range(5):
        for col in range(11):
            st[row, col] = (row + 1) * 1e-3 + col * 1e-5
    return hp, st


def test_reversion_rate_deletion_uses_shorter_run_insertion_rate():
    hp, st = _mats()
    t = _B2N["T"]
    # Deleting one T from TTTTTTT (7) leaves 6; reverting is +1 in HP6.
    assert indelReversionRate(7, 0, -1, t, t, hp, st) == hp[5, t * 3 + 2]


def test_reversion_rate_deletion_from_isolated_base():
    hp, st = _mats()
    g = _B2N["G"]
    # Deleting an isolated G leaves no run: reverting is a non-repeat +1.
    assert indelReversionRate(1, 0, -1, g, g, hp, st) == st[0, 6]


def test_reversion_rate_run_extending_insertion_uses_longer_run():
    hp, st = _mats()
    t = _B2N["T"]
    # Inserting a T into TTTTTT (6) makes 7; reverting is -1 in HP7.
    assert indelReversionRate(6, 0, 1, t, t, hp, st) == hp[6, t * 3 + 0]
    # Capped at HP10+.
    assert indelReversionRate(10, 0, 1, t, t, hp, st) == hp[9, t * 3 + 0]


def test_reversion_rate_mismatched_insertion_is_isolated_base_deletion():
    hp, st = _mats()
    a, t = _B2N["A"], _B2N["T"]
    # T inserted into an A run (HP=10 at the context base): the new T is an
    # isolated base, so reverting is -1 of an HP1 T, not anything A-run.
    assert indelReversionRate(10, 0, 1, a, t, hp, st) == hp[0, t * 3 + 0]


def test_reversion_rate_str_uses_opposite_length():
    hp, st = _mats()
    assert indelReversionRate(1, 2, -2, 0, 0, hp, st) == st[2, 2 + 5]
    assert indelReversionRate(1, 1, 3, 0, 0, hp, st) == st[1, -3 + 5]


def test_reversion_rate_without_matrices_is_nan():
    assert np.isnan(indelReversionRate(3, 0, -1, 0, 0, None, None))


# ------------------------------------------------------ indel_strand_pvalue


def test_indel_strand_evidence_counts_informative_reads():
    assert indel_strand_evidence(5, 2, 1e-3) == (7, 2, 1e-3)


def test_indel_strand_pvalue_no_ref_reads_is_one():
    assert indel_strand_pvalue(6, 0, 1e-3) == 1.0


@pytest.mark.parametrize("n, k, q", [(5, 1, 1e-4), (5, 2, 3.7e-4), (9, 1, 2.9e-3)])
def test_indel_strand_pvalue_is_binomial_upper_tail(n, k, q):
    assert indel_strand_pvalue(n, k, q) == pytest.approx(
        binom.sf(k - 1, n, q), rel=1e-6
    )


def test_indel_strand_pvalue_nan_rate():
    assert np.isnan(indel_strand_pvalue(5, 1, float("nan")))
