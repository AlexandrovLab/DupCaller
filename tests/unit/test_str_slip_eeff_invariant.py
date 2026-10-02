"""STR slips are classified by sequence (str_slip_tract), str.h5 keeps one
tract per position (resolve_str_overlaps), and every indel candidate that
calling would keep sits in its own channel's Eeff site set
(indel_eeff_site_masks) -- the invariant whose break made STR3 Eeff 9 sites
against 819 candidates in rerun 6."""

import itertools

import numpy as np
import pytest

from DupCaller_sub.Index import resolve_str_overlaps
from DupCaller_sub.funcs.indels import (
    indel_context_fits_window,
    indel_context_index,
    indel_has_context,
    indel_passes_mask,
    left_align_indel,
    str_slip_tract,
)
from DupCaller_sub.funcs.misc import indel_eeff_site_masks, str_tract_valid

ACGT = "ATCG"  # base2num order


def _ref_int(seq):
    return np.array([ACGT.index(b) for b in seq])


def _str_raw(n, tracts):
    """str.h5 rows (unit, repeat count, start flag) painted from
    non-overlapping (start, end, unit, count) tracts."""
    raw = np.zeros((3, n), dtype=int)
    for start, end, unit, count in tracts:
        raw[0, start:end] = unit
        raw[1, start:end] = count
        raw[2, start] = 1
    return raw


def _hp_raw(ref_int):
    n = len(ref_int)
    cut = np.ones(n, dtype=bool)
    cut[1:] = ref_int[1:] != ref_int[:-1]
    run_id = np.cumsum(cut) - 1
    return np.vstack((np.bincount(run_id)[run_id], cut.astype(int)))


def _perf_like(seq, max_unit=4):
    """Every maximal perfect repeat of unit 2-max_unit with >= 2 units, as
    PERF reports them: each unit length scanned separately, so tracts of
    different units can overlap."""
    rows = []
    n = len(seq)
    for u in range(2, max_unit + 1):
        i = 0
        while i + 2 * u <= n:
            j = i + u
            while j < n and seq[j] == seq[j - u]:
                j += 1
            count = (j - i) // u
            motif = seq[i : i + u]
            if count >= 2 and len(set(motif)) > 1:
                rows.append((i, i + count * u, u, count))
                i += count * u
            else:
                i += 1
    return rows


# ------------------------------------------------------------ slip rule


def test_slip_rule_by_sequence():
    #      0123 4-13 (CA)5
    seq = "GTGG" + "CA" * 5 + "TTGCA"
    ref = _ref_int(seq)
    raw = _str_raw(len(seq), [(4, 14, 2, 5)])
    assert str_slip_tract(-2, "", 4, ref, raw) == (2, 10)
    assert str_slip_tract(4, "CACA", 4, ref, raw) == (2, 10)
    assert str_slip_tract(-10, "", 4, ref, raw) == (2, 10)  # whole tract
    # Right length, not the repeat: "GG" inserted into (CA)n.
    assert str_slip_tract(2, "GG", 4, ref, raw) is None
    # Out of phase: "AC" inserted before the C that starts the tract.
    assert str_slip_tract(2, "AC", 4, ref, raw) is None
    # Longer than the tract, a non-unit length, mid-tract.
    assert str_slip_tract(-12, "", 4, ref, raw) is None
    assert str_slip_tract(-3, "", 4, ref, raw) is None
    assert str_slip_tract(-2, "", 6, ref, raw) is None
    assert str_slip_tract(-1, "", 4, ref, raw) is None


# ------------------------------------------------------------ str.h5 overlaps


def test_overlapping_tracts_keep_the_longest():
    # chr22:20000006 (AGGA)2 and (AGG)2 overlap; (AGGA)2 is longer.
    assert resolve_str_overlaps([(6, 14, 4, 2), (10, 16, 3, 2)]) == [(6, 14, 4, 2)]
    # chr22:20122723 (TG)20, (GTG)2 at 37-42 and (GT)2 at 40-43.
    rows = [(0, 40, 2, 20), (37, 43, 3, 2), (40, 44, 2, 2)]
    assert resolve_str_overlaps(rows) == [(0, 40, 2, 20), (40, 44, 2, 2)]
    # Abutting tracts both stay; PERF's overextended span is clamped.
    assert resolve_str_overlaps([(0, 4, 2, 2), (4, 9, 2, 2)]) == [
        (0, 4, 2, 2),
        (4, 8, 2, 2),
    ]


def test_resolved_tracts_flag_their_own_values():
    rng = np.random.default_rng(7)
    seq = "".join(rng.choice(list("AC"), size=3000))
    tracts = resolve_str_overlaps(_perf_like(seq))
    raw = _str_raw(len(seq), tracts)
    for start, end, unit, count in tracts:
        assert raw[2, start] == 1
        assert (raw[0, start:end] == unit).all()
        assert (raw[1, start:end] == count).all()
        assert raw[2, start + 1 : end].sum() == 0


# ------------------------------------------------------------ window fit


def test_run_or_tract_cut_by_the_window_edge_is_dropped():
    seq = "GTGC" + "CA" * 5
    ref = _ref_int(seq)
    raw = _str_raw(len(seq), [(4, 14, 2, 5)])
    hp = _hp_raw(ref)
    assert not indel_context_fits_window("3:-2", 3, ref, hp, raw)
    assert indel_context_fits_window("3:-3", 3, ref, hp, raw)  # STR0
    seq = "GTCAAAAA"
    ref = _ref_int(seq)
    hp = _hp_raw(ref)
    raw = np.zeros((3, len(seq)), dtype=int)
    assert not indel_context_fits_window("2:-1", 2, ref, hp, raw)
    assert indel_context_fits_window("2:1:G", 2, ref, hp, raw)  # mismatched


# ------------------------------------------------------------ invariant


def _channel_site(indel, local_pos, ref, hp_raw, str_raw, masks):
    """The Eeff mask an indel candidate is counted against, mirroring
    genotypeDSIndel's routing and Caller.py's channels (None: 5+ pooled,
    which the site masks only approximate)."""
    hp, str_by_len, str0_any, str0_by_len = masks
    parts = indel.split(":")
    length = int(parts[1])
    ctx = indel_context_index(local_pos)
    if abs(length) == 1:
        if length == 1 and _ref_int(parts[2])[0] != ref[ctx]:
            return str0_any
        return hp
    if abs(length) >= 5:
        return None
    sign_k = length
    slip = str_slip_tract(length, parts[2] if length > 0 else "", ctx, ref, str_raw)
    return (str_by_len if slip is not None else str0_by_len)[sign_k]


@pytest.mark.parametrize("seed", range(6))
def test_every_kept_candidate_is_an_eeff_site(seed):
    rng = np.random.default_rng(seed)
    # Repeat-rich: random dinucleotide/trinucleotide blocks and noise.
    blocks = []
    while sum(map(len, blocks)) < 400:
        r = rng.random()
        if r < 0.3:
            unit = "".join(rng.choice(list(ACGT), size=rng.integers(2, 5)))
            blocks.append(unit * int(rng.integers(2, 8)))
        else:
            blocks.append("".join(rng.choice(list(ACGT), size=rng.integers(1, 6))))
    seq = "".join(blocks)[:400]
    ref = _ref_int(seq)
    n = len(seq)
    tracts = resolve_str_overlaps(_perf_like(seq))
    str_raw = _str_raw(n, tracts)
    hp_raw = _hp_raw(ref)
    antimask = rng.random(n) > 0.05
    for start, end, _, _ in tracts[::5]:
        antimask[start:end] = False  # noise intervals over some tracts
    total_len = str_raw[0] * str_raw[1]
    str_cut = str_raw[2].astype(bool)
    hp_cut = hp_raw[1].astype(bool)
    hp_valid = np.ones(n, dtype=bool)
    cuts = np.nonzero(hp_cut)[0]
    hp_valid[cuts[-1] :] = False
    hp_valid[: cuts[1]] = False
    masks = indel_eeff_site_masks(
        antimask,
        ref,
        hp_cut,
        hp_valid & antimask,
        str_raw[0],
        total_len,
        str_cut,
        str_tract_valid(str_cut, total_len) & antimask,
    )
    inserts = [
        "".join(p) for k in (1, 2, 3, 4) for p in itertools.product(ACGT, repeat=k)
    ]
    checked = 0
    for anchor in range(n):
        raw = [f"{anchor}:-{k}" for k in (1, 2, 3, 4)]
        raw += [f"{anchor}:{len(s)}:{s}" for s in rng.choice(inserts, size=12)]
        # Slip insertions at this position, in phase.
        if anchor + 1 < n and str_raw[2, anchor + 1]:
            u = str_raw[0, anchor + 1]
            unit = seq[anchor + 1 : anchor + 1 + u]
            raw += [f"{anchor}:{u * m}:{unit * m}" for m in (1, 2) if u * m <= 4]
        for indel in {left_align_indel(r, ref, 0) for r in raw}:
            pos = int(indel.split(":")[0])
            length = int(indel.split(":")[1])
            if not (
                indel_passes_mask(antimask, pos, length)
                and indel_has_context(pos, n)
                and indel_context_fits_window(indel, pos, ref, hp_raw, str_raw)
            ):
                continue
            site = _channel_site(indel, pos, ref, hp_raw, str_raw, masks)
            if site is None:
                continue
            assert site[indel_context_index(pos)], indel
            checked += 1
    assert checked > 500
