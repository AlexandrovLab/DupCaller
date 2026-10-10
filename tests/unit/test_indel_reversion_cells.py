"""_indel_error_cells reads the reversion (alt->ref) rate in the mutant
molecule's own context:

- a 1bp deletion from a homopolymer of 11+ bases leaves a 10+ run, so
  its reversion is the HP10+ "+1" rate, not HP9's (the run length is
  capped only for the row lookup, after subtracting the deleted base);
- an indel of 2bp or more reverts in the STR bin of the mutant tract
  (reference tract length + indel length), which can differ from the
  reference tract's bin at the 10/25/40bp edges;
- str.txt columns are a (deletion, reference, insertion) triple per
  repeat-unit group: the tract's unit in STR bins 1-4 (any number of
  units), the indel length in bin 0.
"""

import numpy as np
import pytest

from DupCaller_sub.funcs.indels import STR_NCOL, str_col
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
    amp_str = np.arange(5 * STR_NCOL, dtype=float).reshape(5, STR_NCOL) + 3000
    dmg_str = np.arange(5 * STR_NCOL, dtype=float).reshape(5, STR_NCOL) + 4000
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
    reversion, forward = _indel_error_cells(1, 2, 2, 0, 0, strs_mut=3, unit_len=2)
    u2_del, u2_ins = str_col(2, -1), str_col(2, 1)
    assert reversion == ("str", 3, u2_del, u2_del)
    assert forward == ("str", 2, u2_ins, u2_ins)


def test_str_deletion_reverts_in_mutant_bin():
    # -2 from a 26bp tract (bin 3) leaves a 24bp tract (bin 2).
    amp_hp, dmg_hp, amp_str, dmg_str = _mats()
    Pamp, Pamp_rev, *_ = indelErrorProbs(
        1, 3, -2, 0, 0, amp_hp, dmg_hp, amp_str, dmg_str, strs_mut=2, unit_len=2
    )
    assert Pamp == amp_str[2, str_col(2, 1)]
    assert Pamp_rev == amp_str[3, str_col(2, -1)]


def test_str_mutant_bin_defaults_to_reference_bin():
    # Per-context callers (power tables, thresholds) don't know the tract
    # length; both cells then use the reference bin.
    reversion, forward = _indel_error_cells(1, 2, 3, 0, 0, unit_len=3)
    u3_del, u3_ins = str_col(3, -1), str_col(3, 1)
    assert reversion == ("str", 2, u3_del, u3_del)
    assert forward == ("str", 2, u3_ins, u3_ins)


@pytest.mark.parametrize(
    "unit, id_len, group",
    [
        (2, -2, 2),
        (2, -4, 2),  # two units of a dinucleotide: the dinucleotide's column
        (2, 6, 2),
        (3, -3, 3),
        (4, 4, 4),
        (5, -5, 5),
        (6, 6, 5),  # 5+ pools longer units
        (8, -5, 5),  # capped length, unit 8
    ],
)
def test_str_cells_use_tract_unit(unit, id_len, group):
    reversion, forward = _indel_error_cells(1, 2, id_len, 0, 0, unit_len=unit)
    sign = 1 if id_len > 0 else -1
    assert forward[2] == str_col(group, sign)
    assert reversion[2] == str_col(group, -sign)


@pytest.mark.parametrize("id_len, group", [(2, 2), (-3, 3), (4, 4), (5, 5), (-5, 5)])
def test_str0_cells_use_indel_length(id_len, group):
    # Not a slip (bin 0): the unit group is the indel length.
    reversion, forward = _indel_error_cells(1, 0, id_len, 0, 0, unit_len=2)
    sign = 1 if id_len > 0 else -1
    assert forward == ("str", 0, str_col(group, sign), str_col(group, sign))
    assert reversion == ("str", 0, str_col(group, -sign), str_col(group, -sign))


def test_mismatched_insertion_uses_str0_u1():
    reversion, forward = _indel_error_cells(3, 0, 1, A, 2)
    assert forward == ("str", 0, str_col(1, 1), str_col(1, 1))
    reversion, _ = _indel_error_cells(3, 0, 1, A, -1)
    assert reversion == ("str", 0, str_col(1, -1), str_col(1, -1))


def test_str_bin_needs_unit():
    with pytest.raises(ValueError):
        _indel_error_cells(1, 2, 2, 0, 0)
