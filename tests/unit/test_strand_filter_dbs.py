"""Strand independence filter bookkeeping in Caller.py:

- _apply_strand_filter compares the unrounded p-values (MSP prints 4
  significant digits) and skips records with a non-finite p;
- _fail_dbs_of_failed_sbs propagates a failure only within one read
  family (TAG1/TAG2/SP), and a failed DBS takes its other SBS with it;
- --p_threshold must be in [0, 1).
"""

import argparse
import importlib.util
import os

import pytest

from DupCaller_sub.Caller import (
    _apply_strand_filter,
    _dbs_family_keys,
    _fail_dbs_of_failed_sbs,
    _sbs_family_key,
)


def _sbs(
    pos, tag1="AAA", tag2="CCC", sp=100, tl=150, filt="PASS", strand=((0.5,), (0.5,))
):
    return {
        "chrom": "chr1",
        "pos": pos,
        "filter": filt,
        "infos": {
            "TAG1": tag1,
            "TAG2": tag2,
            "SP": sp,
            "TL": tl,
            "MSP": ".",
            "_strand": strand,
        },
    }


def _dbs(pos, tag1="AAA", tag2="CCC", sp=100, tl=150):
    return {
        "chrom": "chr1",
        "pos": pos,
        "filter": "PASS",
        "infos": {"TAG1": tag1, "TAG2": tag2, "SP": sp, "TL": tl},
    }


def _identity(p):
    return p


def test_threshold_compares_unrounded_p():
    # 0.0500029 prints as "0.05" but is above a 0.05 threshold.
    rec = _sbs(10, strand=((0.0500029,), (1.0,)))
    failed = _apply_strand_filter([rec], _identity, 0.05)
    assert failed == []
    assert rec["filter"] == "PASS"
    assert rec["infos"]["MSP"] == "0.05,1"


def test_threshold_fails_at_or_below():
    rec = _sbs(10, strand=((1.0,), (0.05,)))
    assert _apply_strand_filter([rec], _identity, 0.05) == [rec]
    assert rec["filter"] == "strand_independence"


def test_zero_threshold_reports_msp_without_failing():
    rec = _sbs(10, strand=((1e-9,), (1e-9,)))
    assert _apply_strand_filter([rec], _identity, 0) == []
    assert rec["filter"] == "PASS"
    assert rec["infos"]["MSP"] == "1e-09,1e-09"


def test_nonfinite_p_is_not_tested():
    rec = _sbs(10, strand=((float("nan"),), (1e-9,)))
    assert _apply_strand_filter([rec], _identity, 0.05) == []
    assert rec["filter"] == "PASS"
    assert rec["infos"]["MSP"] == "."


def test_record_without_strand_inputs_is_skipped():
    rec = _sbs(10, strand=None)
    assert _apply_strand_filter([rec], _identity, 0.05) == []
    assert rec["infos"]["MSP"] == "."


def test_dbs_keys_are_its_two_sbs():
    dbs = _dbs(10)
    assert _dbs_family_keys(dbs) == [
        _sbs_family_key(_sbs(10)),
        _sbs_family_key(_sbs(11)),
    ]


def test_dbs_fails_with_own_family_sbs_and_fails_partner():
    s10, s11 = _sbs(10), _sbs(11)
    s10["filter"] = "strand_independence"
    dbs = _dbs(10)
    n_dbs, n_partner = _fail_dbs_of_failed_sbs([dbs], [s10, s11], [s10])
    assert (n_dbs, n_partner) == (1, 1)
    assert dbs["filter"] == "strand_independence"
    assert s11["filter"] == "strand_independence"


def test_other_family_failure_at_same_position_is_ignored():
    # A different molecule (other barcodes) fails at position 10; this
    # family's DBS and SBS stay PASS.
    other = _sbs(10, tag1="GGG", tag2="TTT")
    other["filter"] = "strand_independence"
    s10, s11 = _sbs(10), _sbs(11)
    dbs = _dbs(10)
    assert _fail_dbs_of_failed_sbs([dbs], [other, s10, s11], [other]) == (0, 0)
    assert dbs["filter"] == "PASS"
    assert s10["filter"] == s11["filter"] == "PASS"


def test_same_barcodes_and_start_other_template_length_is_another_family():
    # A fragment no longer than the read: R1-fwd (+TL) and R2-rev (-TL)
    # form two families with the same barcodes and start. Only the family
    # whose SBS failed loses its DBS and partner.
    a10, a11 = _sbs(10, tl=120), _sbs(11, tl=120)
    b10, b11 = _sbs(10, tl=-120), _sbs(11, tl=-120)
    a10["filter"] = "strand_independence"
    da, db = _dbs(10, tl=120), _dbs(10, tl=-120)
    n = _fail_dbs_of_failed_sbs([da, db], [a10, a11, b10, b11], [a10])
    assert n == (1, 1)
    assert (
        da["filter"] == "strand_independence" and a11["filter"] == "strand_independence"
    )
    assert db["filter"] == "PASS" and b10["filter"] == b11["filter"] == "PASS"


def test_three_adjacent_sbs_fail_both_dbs():
    s10, s11, s12 = _sbs(10), _sbs(11), _sbs(12)
    s10["filter"] = "strand_independence"
    d10, d11 = _dbs(10), _dbs(11)
    n_dbs, n_partner = _fail_dbs_of_failed_sbs([d11, d10], [s10, s11, s12], [s10])
    assert (n_dbs, n_partner) == (2, 2)
    assert d10["filter"] == d11["filter"] == "strand_independence"
    assert s11["filter"] == s12["filter"] == "strand_independence"


def _p_threshold_value():
    path = os.path.join(os.path.dirname(__file__), "..", "..", "src", "DupCaller.py")
    spec = importlib.util.spec_from_file_location("dupcaller_cli", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.p_threshold_value


@pytest.mark.parametrize("value", ["0", "0.05", "0.999"])
def test_p_threshold_accepts_unit_interval(value):
    assert _p_threshold_value()(value) == float(value)


@pytest.mark.parametrize("value", ["1", "1.5", "-0.01", "nan"])
def test_p_threshold_rejects_outside_unit_interval(value):
    with pytest.raises(argparse.ArgumentTypeError):
        _p_threshold_value()(value)
