"""indelErrorProbs fills the posterior slots like genotypeDSSnv does:
Pamp/Pdmg/Pdmg_bot = alt->ref (reversion in the mutant molecule's
context), Pamp_rev/Pdmg_rev/Pdmg_rev_bot = ref->alt (the indel itself in
the reference context)."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from DupCaller_sub.funcs.prob import (  # noqa: E402
    indelErrorProbs,
    indelMaxLR,
    indelReversionRate,
)

A, T, C, G = 0, 1, 2, 3
RC = {A: T, T: A, C: G, G: C}


@pytest.fixture
def mats():
    # Distinct values in every cell so each lookup is identifiable.
    amp_hp = np.arange(120, dtype=float).reshape(10, 12) + 1000
    dmg_hp = np.arange(120, dtype=float).reshape(10, 12) + 2000
    amp_str = np.arange(55, dtype=float).reshape(5, 11) + 3000
    dmg_str = np.arange(55, dtype=float).reshape(5, 11) + 4000
    return amp_hp, dmg_hp, amp_str, dmg_str


def _probs(mats, hps, strs, id_len, ref, ins):
    amp_hp, dmg_hp, amp_str, dmg_str = mats
    return indelErrorProbs(
        hps, strs, id_len, ref, ins, amp_hp, dmg_hp, amp_str, dmg_str
    )


def test_hp_deletion(mats):
    amp_hp, dmg_hp, _, _ = mats
    Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot = _probs(
        mats, 7, 0, -1, T, T
    )
    # reversion: +1 in a run of 6; forward: -1 in a run of 7
    assert Pamp == amp_hp[5, T * 3 + 2]
    assert Pamp_rev == amp_hp[6, T * 3 + 0]
    assert Pdmg == dmg_hp[5, T * 3 + 2]
    assert Pdmg_rev == dmg_hp[6, T * 3 + 0]
    assert Pdmg_bot == dmg_hp[5, RC[T] * 3 + 2]
    assert Pdmg_rev_bot == dmg_hp[6, RC[T] * 3 + 0]


def test_isolated_base_deletion_reverts_as_mismatched_insertion(mats):
    amp_hp, dmg_hp, amp_str, dmg_str = mats
    Pamp, Pamp_rev, Pdmg, Pdmg_rev, _, _ = _probs(mats, 1, 0, -1, G, G)
    assert Pamp == amp_str[0, 6]
    assert Pdmg == dmg_str[0, 6]
    assert Pamp_rev == amp_hp[0, G * 3 + 0]
    assert Pdmg_rev == dmg_hp[0, G * 3 + 0]


def test_run_extending_insertion(mats):
    amp_hp, dmg_hp, _, _ = mats
    Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot = _probs(mats, 3, 0, 1, C, C)
    # reversion: -1 in a run of 4; forward: +1 in a run of 3
    assert Pamp == amp_hp[3, C * 3 + 0]
    assert Pamp_rev == amp_hp[2, C * 3 + 2]
    assert Pdmg_bot == dmg_hp[3, RC[C] * 3 + 0]
    assert Pdmg_rev_bot == dmg_hp[2, RC[C] * 3 + 2]


def test_insertion_into_capped_run(mats):
    amp_hp, _, _, _ = mats
    Pamp, Pamp_rev, *_ = _probs(mats, 10, 0, 1, A, A)
    assert Pamp == amp_hp[9, A * 3 + 0]
    assert Pamp_rev == amp_hp[9, A * 3 + 2]


def test_mismatched_insertion(mats):
    amp_hp, dmg_hp, amp_str, dmg_str = mats
    Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot = _probs(mats, 5, 0, 1, T, G)
    # reversion: -1 of the inserted G as a run of 1; forward: str row 0 +1
    assert Pamp == amp_hp[0, G * 3 + 0]
    assert Pdmg == dmg_hp[0, G * 3 + 0]
    assert Pdmg_bot == dmg_hp[0, RC[G] * 3 + 0]
    assert Pamp_rev == amp_str[0, 6]
    assert Pdmg_rev == dmg_str[0, 6]
    assert Pdmg_rev_bot == dmg_str[0, 6]


def test_str_indel(mats):
    _, _, amp_str, dmg_str = mats
    Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot = _probs(
        mats, 1, 2, -3, A, A
    )
    assert Pamp == amp_str[2, 3 + 5]
    assert Pamp_rev == amp_str[2, -3 + 5]
    assert Pdmg == Pdmg_bot == dmg_str[2, 3 + 5]
    assert Pdmg_rev == Pdmg_rev_bot == dmg_str[2, -3 + 5]


@pytest.mark.parametrize(
    "hps,strs,id_len,ref,ins",
    [
        (7, 0, -1, T, T),
        (1, 0, -1, G, G),
        (3, 0, 1, C, C),
        (5, 0, 1, T, G),
        (1, 2, -3, A, A),
    ],
)
def test_reversion_rate_matches_pamp(mats, hps, strs, id_len, ref, ins):
    amp_hp, _, amp_str, _ = mats
    Pamp = _probs(mats, hps, strs, id_len, ref, ins)[0]
    assert indelReversionRate(hps, strs, id_len, ref, ins, amp_hp, amp_str) == Pamp


def test_forward_damage_rate_lowers_deletion_ceiling():
    # A higher ref->alt deletion damage rate at the candidate's own run
    # must make a deletion call harder (lower max LR), as for SNVs.
    amp_hp = np.full((10, 12), 1e-5)
    amp_str = np.full((5, 11), 1e-5)
    dmg_str = np.full((5, 11), 1e-6)
    low = np.full((10, 12), 1e-6)
    high = low.copy()
    high[3, T * 3 + 0] = 1e-3
    high[3, RC[T] * 3 + 0] = 1e-3

    def max_lr(dmg_hp):
        _, _, Pd, Pdr, Pdb, Pdrb = indelErrorProbs(
            4, 0, -1, T, T, amp_hp, dmg_hp, amp_str, dmg_str
        )
        return indelMaxLR(Pd, Pdr, Pdb, Pdrb)

    assert max_lr(high) < max_lr(low)
