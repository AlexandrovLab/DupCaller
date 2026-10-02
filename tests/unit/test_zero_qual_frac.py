"""genotypeDSSnv masks a position when at least --maxZeroQualFrac of the
family's reads have no usable base there (deleted, N, off the read, or
BQ <= minBq), for both coverage (antimask) and candidates."""

from types import SimpleNamespace

import numpy as np
import pytest

from DupCaller_sub.funcs.misc import (
    build_trinuc64_order,
    fallback_error_file,
    load_error_matrices,
)
from DupCaller_sub.funcs.prob import genotypeDSSnv

BASE2NUM = {"A": 0, "T": 1, "C": 2, "G": 3}
TRINUC2NUM, _ = build_trinuc64_order()
REF = "ACG"


@pytest.fixture(scope="module")
def base_params():
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


def _read(top, middle=None, bq=37):
    """Read over REF; middle=None deletes the middle base (CIGAR 1M1D1M),
    otherwise that base is shown there."""
    if middle is None:
        seq = REF[0] + REF[2]
        cigar = [(0, 1), (2, 1), (0, 1)]
    else:
        seq = REF[0] + middle + REF[2]
        cigar = [(0, 3)]
    quals = [37] * len(seq)
    if middle is not None:
        quals[1] = bq
    return SimpleNamespace(
        is_read1=top,
        is_read2=not top,
        is_forward=True,
        is_reverse=False,
        query_alignment_qualities=quals,
        query_alignment_sequence=seq,
        cigartuples=cigar,
        reference_start=0,
        reference_length=3,
    )


def _genotype(params, reads):
    ref = np.array([BASE2NUM[b] for b in REF])
    trinuc = np.array([64, TRINUC2NUM[REF], 64])
    out = genotypeDSSnv(
        reads, 0, ref, trinuc, np.ones(3, bool), np.ones(3, bool), params, None
    )
    _, _, _, mut_antimask, _, antimask = out[:6]
    return mut_antimask[1], antimask[1]


def _deletion_family():
    # Per strand: 2 reads deleting the middle base, 1 showing C>T.
    return [_read(True), _read(True), _read(True, "T")] + [
        _read(False),
        _read(False),
        _read(False, "T"),
    ]


def test_default_threshold_is_0_9(base_params):
    # 4 of 6 reads deleted (0.67): kept at the default 0.9.
    candidate, covered = _genotype(base_params, _deletion_family())
    assert candidate and covered
    # 9 of 10 reads deleted: masked.
    reads = [_read(True) for _ in range(5)] + [_read(False) for _ in range(4)]
    reads.append(_read(False, "T"))
    candidate, covered = _genotype(base_params, reads)
    assert not candidate and not covered


def test_mostly_deleted_position_is_masked_at_half(base_params):
    params = dict(base_params, maxZeroQualFrac=0.5)
    candidate, covered = _genotype(params, _deletion_family())
    assert not candidate and not covered


def test_mostly_deleted_position_kept_when_filter_relaxed(base_params):
    params = dict(base_params, maxZeroQualFrac=1.0)
    candidate, covered = _genotype(params, _deletion_family())
    assert candidate and covered


def test_one_strand_all_low_bq_is_masked_at_half(base_params):
    # 2+2 family whose bottom-strand bases are all BQ <= minBq: half the
    # reads have no usable base.
    params = dict(base_params, maxZeroQualFrac=0.5)
    reads = [
        _read(True, "C"),
        _read(True, "C"),
        _read(False, "C", 5),
        _read(False, "C", 5),
    ]
    _, covered = _genotype(params, reads)
    assert not covered


def test_clean_position_is_not_masked(base_params):
    reads = [_read(True, "C") for _ in range(3)] + [_read(False, "C") for _ in range(3)]
    _, covered = _genotype(base_params, reads)
    assert covered
