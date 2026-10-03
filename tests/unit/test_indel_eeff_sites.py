"""Indel Eeff (depth_by_hpstr) site set: one site per homopolymer run (its
first base, where every left-aligned event in it is called), runs cut by
the window edge excluded; STR channels count the starts of tracts whose
unit divides the indel length; STR0 insertions every position, STR0
deletions every non-slip deletion site; deletions and insertions are
counted separately; the anchor and context base (and every deleted base)
must pass the indel masks; reads on both strands.
Also the mu solve's root search and the -mr table version check."""

import numpy as np
import pandas as pd
import pytest

from DupCaller_sub.Caller import INDEL_EEFF_SITE_SET, _load_mutation_rate_override
from DupCaller_sub.funcs.indels import str_unit_multiple
from DupCaller_sub.funcs.misc import (
    DEPTH_HPSTR_NCAT,
    DEPTH_HPSTR_STR0_ANY,
    _channel_eeff_at_threshold,
    depth_hpstr_str0_bucket,
    depth_hpstr_str_bucket,
    indel_eeff_site_masks,
    indel_mu_channel,
    init_refine_worker,
    MU0_NO_ROOT_RATE,
    refine_channel_task,
    str_tract_valid,
)


def _run_boundary_valid(cut):
    # Same rule as call.py's _run_boundary_valid (homopolymer runs).
    valid = np.ones(cut.shape[0], dtype=bool)
    pos = np.nonzero(cut)[0]
    if pos.size:
        valid[pos[-1] :] = False
        if pos[0] == 0:
            first_end = pos[1] if pos.size > 1 else cut.shape[0]
        else:
            first_end = pos[0]
        valid[:first_end] = False
    return valid


def _masks(ref, antimask, tracts=()):
    """tracts: (start, unit, repeat count) STR annotations."""
    n = len(ref)
    ref_int = np.array(["ATCG".index(b) for b in ref])
    hp_cut = np.ones(n, dtype=bool)
    hp_cut[1:] = ref_int[1:] != ref_int[:-1]
    unit = np.zeros(n, dtype=int)
    total_len = np.zeros(n, dtype=int)
    str_cut = np.zeros(n, dtype=bool)
    for start, u, count in tracts:
        unit[start : start + u * count] = u
        total_len[start : start + u * count] = u * count
        str_cut[start] = True
    return indel_eeff_site_masks(
        antimask,
        ref_int,
        hp_cut,
        _run_boundary_valid(hp_cut) & antimask,
        unit,
        total_len,
        str_cut,
        str_tract_valid(str_cut, total_len) & antimask,
    )


def test_hp_run_counts_once_at_its_first_base():
    #      0123456789...
    ref = "GCAAAAAAAAAAGCGCT"
    hp, _, any_, _ = _masks(ref, np.ones(len(ref), dtype=bool))
    # The 10-A run (2-11) is one site, at 2; every other internal run start
    # counts once; the first run (0) and the last run (16) touch the window
    # edges, and 1 is one of the window's first two bases (no candidate's
    # context base can sit there).
    assert np.nonzero(hp)[0].tolist() == [2, 12, 13, 14, 15]
    # STR0 1bp: every position except the window's first two.
    assert np.nonzero(any_)[0].tolist() == list(range(2, len(ref)))


def test_masked_anchor_or_context_is_not_a_site():
    ref = "GCAAAAAGCGCTGC"
    antimask = np.ones(len(ref), dtype=bool)
    antimask[6] = False  # anchor of a site at 7
    hp, _, any_, _ = _masks(ref, antimask)
    assert not any_[6] and not any_[7]
    assert not hp[7]
    assert any_[8]


def test_single_str_tract_counts_by_unit_length():
    # (CA)5 at 4-13 is the window's only tract; it ends well inside.
    ref = "GTGC" + "CA" * 5 + "GTTCAGGTC"
    _, st, _, st0 = _masks(ref, np.ones(len(ref), dtype=bool), [(4, 2, 5)])
    for k in (2, 4, 5):  # whole-unit slips; pooled 5+ (6, 8, ...)
        assert np.nonzero(st[-k])[0].tolist() == [4]
        assert np.nonzero(st[k])[0].tolist() == [4]
    assert not st[-3].any() and not st[3].any()  # not a dinucleotide slip
    # STR0 insertions: anywhere (a non-matching insertion can go anywhere).
    for k in (2, 3, 5):
        assert st0[k][2:].all()
    # STR0 3bp deletions: anywhere else in the tract. Not at its first base:
    # the C before it makes deleting CAC at 4 the same as CCA at 3, which
    # left-aligns there.
    assert st0[-3][5:14].all() and not st0[-3][4] and st0[-3][3]
    # STR0 2bp deletions: not at the slip site and not inside the tract
    # (they'd left-align), but at its last base (it runs out of the tract).
    assert not st0[-2][4:13].any() and st0[-2][13]


def test_deletion_site_needs_every_deleted_base_unmasked():
    ref = "GTGC" + "CA" * 5 + "GTTCAGGTC"
    antimask = np.ones(len(ref), dtype=bool)
    antimask[6] = False  # third base of the tract
    _, st, _, _ = _masks(ref, antimask, [(4, 2, 5)])
    assert st[-2][4] and st[2][4]  # deletes 4-5, inserts before 4
    assert not st[-4][4]  # deletes 4-7, 6 masked
    assert st[4][4]


def test_deletion_longer_than_the_tract_is_not_a_slip_site():
    ref = "GTGC" + "CA" * 2 + "GTTCAGGTC"
    _, st, _, _ = _masks(ref, np.ones(len(ref), dtype=bool), [(4, 2, 2)])
    assert st[-4][4] and not st[-5][4]
    assert st[5][4]


def test_tract_reaching_the_window_end_is_not_a_site():
    ref = "GTGC" + "CA" * 5
    _, st, _, _ = _masks(ref, np.ones(len(ref), dtype=bool), [(4, 2, 5)])
    assert not st[2].any() and not st[-2].any()


def test_window_shorter_than_a_deletion_has_no_deletion_sites():
    for ref in ("GC", "GCA", "GCAT"):
        _, st, _, st0 = _masks(ref, np.ones(len(ref), dtype=bool))
        assert not any(st0[-k].any() for k in (3, 4, 5) if k >= len(ref))


def test_str_unit_multiple():
    assert str_unit_multiple(-2, 2) and str_unit_multiple(4, 2)
    assert not str_unit_multiple(3, 2) and not str_unit_multiple(-2, 3)
    assert str_unit_multiple(1, 2) and str_unit_multiple(-1, 3)


def test_channel_eeff_needs_both_strands_and_uses_length_buckets():
    depth_hpstr = np.zeros((10, 10, DEPTH_HPSTR_NCAT), dtype=np.int64)
    a_hp3 = 0 * 10 + 3 - 1
    depth_hpstr[2, 3, a_hp3] = 5
    depth_hpstr[4, 0, a_hp3] = 7  # no bottom-strand reads: not callable
    depth_hpstr[1, 1, depth_hpstr_str_bucket(2, -3)] = 11
    depth_hpstr[1, 1, depth_hpstr_str_bucket(2, 2)] = 19
    depth_hpstr[3, 3, DEPTH_HPSTR_STR0_ANY] = 13
    depth_hpstr[3, 3, depth_hpstr_str0_bucket(-3)] = 17
    depth_hpstr[3, 3, depth_hpstr_str0_bucket(7)] = 23
    init_refine_worker(np.zeros((10, 10, 64), dtype=np.int64), depth_hpstr)
    assert _channel_eeff_at_threshold("hp", (3, -1, "T")) == 5
    # Deletions and insertions of one length have their own buckets.
    assert _channel_eeff_at_threshold("str", (2, -3)) == 11
    assert _channel_eeff_at_threshold("str", (2, 3)) == 0
    assert _channel_eeff_at_threshold("str", (2, 2)) == 19
    assert _channel_eeff_at_threshold("str", (2, -2)) == 0
    assert _channel_eeff_at_threshold("str", (0, 1)) == 13
    assert _channel_eeff_at_threshold("str", (0, -3)) == 17
    assert _channel_eeff_at_threshold("str", (0, 3)) == 0
    assert _channel_eeff_at_threshold("str", (0, 5)) == 23  # 5+ pools 7


def test_buckets_are_distinct_and_in_range():
    buckets = [b * 10 + L for b in range(4) for L in range(10)]
    lengths = (-5, -4, -3, -2, 2, 3, 4, 5)
    buckets += [depth_hpstr_str_bucket(b, k) for b in range(1, 5) for k in lengths]
    buckets += [depth_hpstr_str0_bucket(k) for k in (1,) + lengths]
    assert len(set(buckets)) == len(buckets) == DEPTH_HPSTR_NCAT
    assert max(buckets) == DEPTH_HPSTR_NCAT - 1


def _solve(raw_lr, eeff, pseudocount=0.5, neg_log_lr=None):
    depth_hpstr = np.zeros((10, 10, DEPTH_HPSTR_NCAT), dtype=np.int64)
    depth_hpstr[1, 1, 0] = eeff
    init_refine_worker(np.zeros((10, 10, 64), dtype=np.int64), depth_hpstr)
    job = (
        "HP1_len-1_T",
        "hp",
        (1, -1, "T"),
        raw_lr,
        1.0,
        0.05,
        pseudocount,
        None,
        np.empty(0, dtype=np.float32) if neg_log_lr is None else neg_log_lr,
    )
    return refine_channel_task(job)[4]


def test_mu_solve_refuses_zero_pseudocount():
    # pseudocount 0 makes g(0) NaN, which brentq can't bracket.
    with pytest.raises(AssertionError, match="pseudocount"):
        _solve(np.concatenate([np.full(100, 0.01), [5.0]]), 60, pseudocount=0)


def test_cli_rejects_nonpositive_pseudocount():
    import argparse
    import importlib.util
    import os

    spec = importlib.util.spec_from_file_location(
        "dupcaller_cli",
        os.path.join(os.path.dirname(__file__), "..", "..", "src", "DupCaller.py"),
    )
    cli = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    assert cli.positive_float("0.5") == 0.5
    for bad in ("0", "-1"):
        with pytest.raises(argparse.ArgumentTypeError):
            cli.positive_float(bad)


def test_mu_solve_finds_interior_root_when_candidates_reach_eeff():
    # n == Eeff with 5 strong calls: g(1) > 0, but g dips below 0 near
    # mu ~ 0.07; the old shortcut returned 0.
    raw_lr = np.concatenate([np.full(5, 1e6), np.full(95, 1e-3)])
    mu0 = _solve(raw_lr, 100)
    assert 0.01 < mu0 < 0.2


def test_mu_solve_without_a_root_uses_fixed_rate_when_sites_exist():
    # One strong call on one site (the mock HP8 case): g ~ 1.5/mu - 1 > 0
    # on all of (0, 1), so no root; the channel has a site, so 3.5e-9.
    assert _solve(np.array([1e8]), 1) == MU0_NO_ROOT_RATE == 3.5e-9


def test_mu_solve_without_a_root_or_sites_stays_zero():
    # No sites and no coverage: g = pseudocount/mu > 0, no root, mu0 = 0.
    assert _solve(np.empty(0), 0) == 0.0


def test_mu_solve_small_eeff_takes_the_pseudocount_root():
    # No real evidence on one site: the pseudocount alone gives a root at
    # about pseudocount / Eeff (the old shortcut returned 0 here).
    assert _solve(np.full(3, 1e-3), 1) == pytest.approx(0.5, abs=0.01)


def test_mu_solve_float32_log_negatives_match_raw():
    # LR < 0 sites as float32 log10 LRs give the same mu0 as raw LRs.
    rng = np.random.default_rng(0)
    neg = rng.uniform(-8, 0, 5000)
    pos = 10.0 ** np.array([2.5, 4.0, 7.0])
    eeff = 1e7
    mu_raw = _solve(np.concatenate([pos, 10.0**neg]), eeff)
    mu_split = _solve(pos, eeff, neg_log_lr=neg.astype(np.float32))
    assert mu_split == pytest.approx(mu_raw, rel=1e-6)


def test_indel_mu_channel_routing():
    assert indel_mu_channel(3, 1, 2, 1, "A") == "STR2_len3"
    assert indel_mu_channel(-7, 4, 0, 1, "A") == "STR0_len-5"
    assert indel_mu_channel(1, 3, 0, 0, "C") == "STR0_len1"
    assert indel_mu_channel(1, 3, 0, 1, "G") == "HP3_len1_C"
    assert indel_mu_channel(-1, 10, 1, 1, "T") == "HP10_len-1_T"
    assert indel_mu_channel(-1, 0, 0, 1, "A") is None


def test_mu_override_refuses_tables_without_eeff_stamp(tmp_path):
    prefix = str(tmp_path / "s")
    pd.DataFrame(
        {
            "context": ["SBS96:A[C>A]A"],
            "n_sites": [0],
            "Eeff": [1.0],
            "mutation_rate_mle": [1e-6],
        }
    ).to_csv(prefix + "_sbs96_rate_n1.txt", sep="\t", index=False)
    old = pd.DataFrame(
        {
            "context": ["HP1_T"],
            "indel_length": [-1],
            "n_sites": [0],
            "Eeff": [1.0],
            "mutation_rate_mle": [1e-6],
        }
    )
    old.to_csv(prefix + "_indel_rate_by_hp_str.txt", sep="\t", index=False)
    with pytest.raises(ValueError, match="eeff_sites"):
        _load_mutation_rate_override(prefix)
    old["eeff_sites"] = INDEL_EEFF_SITE_SET
    old.to_csv(prefix + "_indel_rate_by_hp_str.txt", sep="\t", index=False)
    assert _load_mutation_rate_override(prefix)["HP1_len-1_T"] == 1e-6
