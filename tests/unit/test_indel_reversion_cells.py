"""_indel_error_cells reads the reversion (alt->ref) rate in the mutant
molecule's own context:

- a 1bp deletion from a homopolymer of 11+ bases leaves a 10+ run, so
  its reversion is the HP10+ "+1" rate, not HP9's (the run length is
  capped only for the row lookup, after subtracting the deleted base);
- an indel of 2bp or more reverts in the STR bin of the mutant tract
  (reference tract length + indel length), which can differ from the
  reference tract's bin at the 10/25/40bp edges.
"""

import numpy as np
import pytest

from DupCaller_sub.funcs.prob import (
    _indel_error_cells,
    indelErrorProbs,
    indelReversionRate,
    str_length_bin,
)

A = 0


@pytest.mark.parametrize(
    "run_len, reversion_row",
    [
        (2, 0),  # mutant run of 1 -> HP1
        (10, 8),  # mutant run of 9 -> HP9
        (11, 9),  # mutant run of 10 -> HP10+
        (15, 9),  # mutant run of 14 -> HP10+
        (127, 9),
    ],
)
def test_deletion_reversion_row_uses_uncapped_run(run_len, reversion_row):
    reversion, forward = _indel_error_cells(run_len, 0, -1, A, A)
    assert reversion == ("hp", reversion_row, A * 3 + 2, 1 * 3 + 2)
    assert forward == ("hp", min(run_len, 10) - 1, A * 3 + 0, 1 * 3 + 0)


@pytest.mark.parametrize("run_len", [9, 10, 11, 15])
def test_insertion_reversion_row_is_capped_longer_run(run_len):
    reversion, forward = _indel_error_cells(run_len, 0, 1, A, A)
    assert reversion == ("hp", min(run_len + 1, 10) - 1, A * 3 + 0, 1 * 3 + 0)
    assert forward == ("hp", min(run_len, 10) - 1, A * 3 + 2, 1 * 3 + 2)


def _mats():
    amp_hp = np.arange(120, dtype=float).reshape(10, 12) + 1000
    dmg_hp = np.arange(120, dtype=float).reshape(10, 12) + 2000
    amp_str = np.arange(55, dtype=float).reshape(5, 11) + 3000
    dmg_str = np.arange(55, dtype=float).reshape(5, 11) + 4000
    return amp_hp, dmg_hp, amp_str, dmg_str


def test_long_run_deletion_rates_match_hp10_row():
    amp_hp, dmg_hp, amp_str, dmg_str = _mats()
    Pamp, Pamp_rev, *_ = indelErrorProbs(
        15, 0, -1, A, A, amp_hp, dmg_hp, amp_str, dmg_str
    )
    assert Pamp == amp_hp[9, A * 3 + 2]
    assert Pamp_rev == amp_hp[9, A * 3 + 0]
    assert indelReversionRate(15, 0, -1, A, A, amp_hp, amp_str) == Pamp


@pytest.mark.parametrize(
    "unit, total, expected",
    [
        (1, 20, 0),  # homopolymer, not an STR
        (2, 3, 0),  # fewer than two units
        (2, 4, 1),
        (2, 9, 1),
        (2, 10, 2),
        (2, 24, 2),
        (2, 25, 3),
        (2, 26, 3),
        (2, 39, 3),
        (2, 40, 4),
    ],
)
def test_str_length_bin(unit, total, expected):
    assert str_length_bin(unit, total) == expected


def test_str_insertion_reverts_in_mutant_bin():
    # +2 into a 24bp tract (bin 2) makes a 26bp tract (bin 3).
    reversion, forward = _indel_error_cells(1, 2, 2, 0, 0, strs_mut=3)
    assert reversion == ("str", 3, 3, 3)
    assert forward == ("str", 2, 7, 7)


def test_str_deletion_reverts_in_mutant_bin():
    # -2 from a 26bp tract (bin 3) leaves a 24bp tract (bin 2).
    amp_hp, dmg_hp, amp_str, dmg_str = _mats()
    Pamp, Pamp_rev, *_ = indelErrorProbs(
        1, 3, -2, 0, 0, amp_hp, dmg_hp, amp_str, dmg_str, strs_mut=2
    )
    assert Pamp == amp_str[2, 2 + 5]
    assert Pamp_rev == amp_str[3, -2 + 5]


def test_str_mutant_bin_defaults_to_reference_bin():
    # Per-context callers (power tables, thresholds) don't know the tract
    # length; both cells then use the reference bin.
    reversion, forward = _indel_error_cells(1, 2, 3, 0, 0)
    assert reversion == ("str", 2, 2, 2)
    assert forward == ("str", 2, 8, 8)
