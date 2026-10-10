"""--skipCoveragePass power tallies: counting positions by power-table index
in round 1 and contracting with the final tables must give what the
coverage pass sums position by position (call.py's cov_mat/indel100)."""

import numpy as np
import pytest

from DupCaller_sub.funcs.misc import _compute_dbs_opportunity, build_trinuc64_order
from DupCaller_sub.funcs.power_tally import (
    PowerTally,
    contract_power_tally,
    indel_keys,
    sbs_power_index,
)

BASES = "ATCG"
TRINUC2NUM, _ = build_trinuc64_order()


def _window(rng, n):
    seq = rng.integers(0, 4, n + 2)
    seq[rng.random(n + 2) < 0.02] = 4  # N
    ref_int = seq[1:-1]
    trinuc = np.full(n, 64)
    for i in range(n):
        tri = seq[i : i + 3]
        if (tri <= 3).all():
            trinuc[i] = TRINUC2NUM["".join(BASES[b] for b in tri)]
    # Read-family counts per base (4, n): mostly constant over the
    # family, thinner at the read ends.
    top = np.zeros((4, n), dtype=int)
    bot = np.zeros((4, n), dtype=int)
    top[0] = rng.integers(0, 12)
    bot[0] = rng.integers(0, 12)
    ends = int(rng.integers(0, min(6, n)))
    top[0, :ends] = rng.integers(0, 3, ends)
    bot[0, n - ends :] = rng.integers(0, 3, ends)
    return ref_int, trinuc, top, bot


def _random_tables(rng):
    L = rng.random((10, 10, 64, 4))
    _, num2trinuc = build_trinuc64_order()
    for t, tri in enumerate(num2trinuc):
        L[:, :, t, BASES.index(tri[1])] = 0
    return L, rng.random((10, 10, 12, 4, 2)), rng.random((10, 10, 5, 11))


def test_sbs_and_dbs_tallies_match_direct_coverage():
    rng = np.random.default_rng(1)
    L, L1, Llen = _random_tables(rng)
    tally = PowerTally()
    tally.FLUSH = 50  # exercise the buffered flush
    direct_trinuc, direct_dbs = {}, {}
    direct_cov = direct_unmasked = 0.0
    for fam in range(300):
        n = int(rng.integers(2, 60))
        ref_int, trinuc, top, bot = _window(rng, n)
        valid = rng.random(n) < 0.9
        unmasked_pass = valid & (rng.random(n) < 0.95)
        pass_bool = unmasked_pass & (rng.random(n) < 0.9)
        group = f"{fam % 7}+{fam % 5}"
        # call.py / genotypeDSSnv, position by position
        n_top = np.minimum(top.sum(axis=0), 9)
        n_bot = np.minimum(bot.sum(axis=0), 9)
        cov_mat = np.zeros((n, 4))
        ok = valid & (trinuc < 64)
        cov_mat[ok] = L[n_top[ok], n_bot[ok], trinuc[ok]]
        t = direct_trinuc.setdefault(group, np.zeros((64, 4)))
        for b in range(4):
            t[:, b] += np.bincount(
                trinuc[pass_bool & ok], weights=cov_mat[pass_bool & ok, b], minlength=64
            )
        direct_dbs[group] = direct_dbs.get(group, 0) + _compute_dbs_opportunity(
            cov_mat, ref_int, pass_bool
        )
        direct_cov += cov_mat[pass_bool].sum()
        direct_unmasked += cov_mat[unmasked_pass].sum()
        # round-1 tally
        pidx = sbs_power_index(top, bot, trinuc, valid)
        tally.add_sbs(group, pidx, pass_bool, unmasked_pass, ref_int)
    trinuc_c, _, dbs_c, cov_c, unmasked_c = contract_power_tally(
        tally.summary(np.array([30.0])), L, L1, Llen
    )
    assert cov_c == pytest.approx(direct_cov, rel=1e-12)
    assert unmasked_c == pytest.approx(direct_unmasked, rel=1e-12)
    for group in direct_trinuc:
        np.testing.assert_allclose(trinuc_c[group], direct_trinuc[group], rtol=1e-12)
        np.testing.assert_allclose(dbs_c[group], direct_dbs[group], rtol=1e-12)
    assert sum(len(v[0]) for v in tally.summary(None)["dbs_other"].values()) > 0


def _direct_indel100(
    n_top,
    n_bot,
    hp_raw,
    ref_safe,
    hp_cut,
    hp_valid,
    str_unit,
    str_rep,
    str_cut,
    str_valid,
    antimask,
    next_ref,
    next_valid,
    L1,
    Llen,
):
    """call.py's indel100 accumulation, term by term."""
    v = np.zeros(100)
    to_cgta = np.array([3, 2, 0, 1])
    base = to_cgta[ref_safe]
    hp10 = np.minimum(hp_raw, 10)
    hp11 = np.minimum(hp_raw, 11)
    is_str = str_unit >= 2
    total = str_unit * str_rep
    sb = np.where(is_str, 1, 0)
    sb[is_str & (total >= 10)] = 2
    sb[is_str & (total >= 25)] = 3
    sb[is_str & (total >= 40)] = 4
    uc = np.clip(str_unit, 2, 5)

    def acc(col, w, g):
        np.add.at(v, np.broadcast_to(col, g.shape)[g], np.broadcast_to(w, g.shape)[g])

    acc(
        base * 6 + np.clip(hp10, 1, 6) - 1,
        L1[n_top, n_bot, hp11, ref_safe, 0],
        hp_valid & hp_cut,
    )
    acc(
        24 + base * 5 + np.clip(hp10, 1, 5) - 1,
        L1[n_top, n_bot, hp11, ref_safe, 1],
        (hp10 >= 1) & hp_valid & hp_cut,
    )
    g = is_str & str_valid & str_cut
    acc(
        44 + (uc - 2) * 6 + np.clip(str_rep, 1, 6) - 1,
        Llen[n_top, n_bot, sb, -uc + 5],
        g,
    )
    for u, U in enumerate([2, 3, 4, 5]):
        acc(44 + u * 6, Llen[n_top, n_bot, 0, -U + 5], antimask)
    acc(68 + (uc - 2) * 6 + np.clip(str_rep, 0, 5), Llen[n_top, n_bot, sb, uc + 5], g)
    acc(68 + (uc - 2) * 6, Llen[n_top, n_bot, sb, uc + 5], g & (str_rep == 1))
    for u, U in enumerate([2, 3, 4, 5]):
        acc(68 + u * 6, Llen[n_top, n_bot, 0, U + 5], antimask)
        acc(68 + u * 6 + 1, Llen[n_top, n_bot, 0, U + 5], antimask)
    for k in range(4):
        acc(92 + k, Llen[n_top, n_bot, 0, 3 - k], antimask)
    for b in range(4):
        acc(
            96 + b,
            Llen[n_top, n_bot, 0, 6] * next_valid * (next_ref != b) * antimask,
            np.ones_like(antimask),
        )
    return v


def test_indel_tallies_match_direct_indel100():
    rng = np.random.default_rng(2)
    _, L1, Llen = _random_tables(rng)
    L = np.zeros((10, 10, 64, 4))
    tally = PowerTally()
    tally.FLUSH = 100
    direct = {}
    for fam in range(300):
        n = int(rng.integers(1, 60))
        n_top = rng.integers(0, 10, n)
        n_bot = rng.integers(0, 10, n)
        hp_raw = rng.integers(0, 14, n)
        ref_safe = rng.integers(0, 4, n)
        hp_cut = rng.random(n) < 0.5
        hp_valid = rng.random(n) < 0.8
        str_unit = rng.integers(0, 7, n)
        str_rep = rng.integers(0, 9, n)
        str_cut = rng.random(n) < 0.5
        str_valid = rng.random(n) < 0.8
        antimask = rng.random(n) < 0.9
        next_ref = rng.integers(-1, 4, n)
        next_valid = (next_ref >= 0) & (next_ref <= 3)
        group = f"{fam % 4}+{fam % 3}"
        direct[group] = direct.get(group, 0) + _direct_indel100(
            n_top,
            n_bot,
            hp_raw,
            ref_safe,
            hp_cut,
            hp_valid,
            str_unit,
            str_rep,
            str_cut,
            str_valid,
            antimask,
            next_ref,
            next_valid,
            L1,
            Llen,
        )
        is_str = str_unit >= 2
        total = str_unit * str_rep
        sb = np.where(is_str, 1, 0)
        sb[is_str & (total >= 10)] = 2
        sb[is_str & (total >= 25)] = 3
        sb[is_str & (total >= 40)] = 4
        tally.add_indel(
            group,
            indel_keys(
                n_top,
                n_bot,
                np.minimum(hp_raw, 11),
                ref_safe,
                np.minimum(hp_raw, 10),
                hp_valid,
                hp_cut,
                is_str,
                str_valid,
                str_cut,
                sb,
                str_unit,
                str_rep,
                antimask,
                next_ref,
                next_valid,
            ),
        )
    _, indel_c, _, _, _ = contract_power_tally(tally.summary(None), L, L1, Llen)
    for group in direct:
        np.testing.assert_allclose(indel_c[group], direct[group], rtol=1e-12)
