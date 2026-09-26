import numpy as np

from DupCaller_sub.funcs.learn import estimate_sbs_srd_rates
from DupCaller_sub.funcs.misc import (
    _normalize_indel_hp_mat,
    _normalize_indel_str_mat,
    _read_fallback_counts,
    apply_sbs_low_coverage_fallback,
    build_trinuc64_order,
    fallback_error_file,
)

BASE2NUM = {"A": 0, "T": 1, "C": 2, "G": 3}
_, NUM2TRINUC = build_trinuc64_order()


def _uniform_rates(alt_rate):
    mat = np.full([64, 4], alt_rate)
    for row, trinuc in enumerate(NUM2TRINUC):
        mat[row, BASE2NUM[trinuc[1]]] = 1 - 3 * alt_rate
    return mat


def _counts(ref_count, alt_count):
    mat = np.full([64, 4], float(alt_count))
    for row, trinuc in enumerate(NUM2TRINUC):
        mat[row, BASE2NUM[trinuc[1]]] = ref_count
    return mat


def test_sbs_only_zero_alt_entries_in_low_rows():
    rates = _uniform_rates(1e-6)
    fallback = _uniform_rates(1e-3)
    counts = _counts(5000, 1)
    counts[0] = [0, 1, 500, 1]  # ACA: low row, ACA>A unobserved
    counts[17] = [0, 0, 0, 0]  # CCT: low row, everything unobserved
    counts[20, 0] = 0  # CTA: A unobserved but row well sampled
    out = apply_sbs_low_coverage_fallback(rates, counts, fallback)
    assert out[0, 0] == 1e-3
    assert out[0, 1] == 1e-6 and out[0, 3] == 1e-6  # observed alts kept
    assert np.array_equal(out[17, [0, 1, 3]], fallback[17, [0, 1, 3]])
    assert np.array_equal(out[20], rates[20])
    assert np.allclose(out.sum(axis=1), 1)
    touched = np.zeros(64, dtype=bool)
    touched[[0, 17]] = True
    assert np.array_equal(out[~touched], rates[~touched])
    # Exactly 1000 sites is enough, even with a zero alt.
    counts = _counts(1000, 0)
    assert np.array_equal(
        apply_sbs_low_coverage_fallback(rates, counts, fallback), rates
    )


def test_bundled_fallback_files_found():
    for suffix in (".amp.tn.srd.txt", ".dmg.tn.txt", ".amp.hp.txt", ".dmg.str.txt"):
        assert fallback_error_file(suffix).endswith("fallback_latest" + suffix)


def test_srd_em_fallback_only_when_requested():
    hist = np.zeros([64, 4, 42])
    for row, trinuc in enumerate(NUM2TRINUC):
        hist[row, BASE2NUM[trinuc[1]], 37] = 5000
    hist[3, BASE2NUM["C"], 37] = 500  # ACG: under-sampled
    hist[3, BASE2NUM["T"], 30] = 2  # ACG>T observed
    plain = estimate_sbs_srd_rates(hist, 0.5)
    floored = estimate_sbs_srd_rates(hist, 0.5, fallback=True)
    fb = _read_fallback_counts(".amp.tn.srd.txt")
    for c in ("A", "G"):
        assert floored[3, BASE2NUM[c]] == fb[3, BASE2NUM[c]]
    assert floored[3, BASE2NUM["T"]] == plain[3, BASE2NUM["T"]]
    assert np.isclose(floored[3].sum(), 1)
    rest = np.arange(64) != 3
    assert np.array_equal(floored[rest], plain[rest])


def _hp_rates(mat, pc=0.5):
    out = np.zeros_like(mat, dtype=float)
    for g in range(4):
        b = mat[:, g * 3 : g * 3 + 3]
        out[:, g * 3 : g * 3 + 3] = (b + pc) / (b.sum(axis=1, keepdims=True) + 3 * pc)
    return out


def test_hp_fallback_chains_on_used_rate():
    fb = np.zeros([10, 12])
    fb[:, 1::3] = 1e6  # idLen=0 opportunity
    fb[:, 0::3] = 10  # deletions
    fb[:, 2::3] = 10  # insertions
    fb[3, 0] = 1000  # HP4 A-deletion: high fallback rate
    obs = fb * 0
    obs[:, 1::3] = 1e5
    obs[:, 0::3] = 1
    obs[:, 2::3] = 1
    # HP4/HP6 contexts under 1000 sites (per base group); T/C/G there
    # have every category observed, so only A's zero entries fall back.
    obs[[3, 5], 3:] = np.tile([1, 50, 1], 3)
    obs[3, 0:3] = [0, 700, 3]  # HP4 A: del unobserved, ins observed
    obs[5, 0:3] = [0, 500, 0]  # HP6 A: del+ins unobserved
    out = _normalize_indel_hp_mat(obs, 0.5, fb)
    fb_rates = _hp_rates(fb)
    obs_rates = _hp_rates(obs)
    # HP4 del takes max(fallback, HP3 used); observed ins keeps its own
    # (monotone-floored) rate, not the fallback.
    assert out[3, 0] == fb_rates[3, 0]
    assert out[3, 2] == max(obs_rates[3, 2], out[2, 2])
    # HP5 is well sampled but monotone-floored by HP4's *used* rate, and
    # HP6's unobserved del is compared against that used HP5 rate.
    assert out[4, 0] == out[3, 0]
    assert out[5, 0] == max(fb_rates[5, 0], out[4, 0])
    assert out[5, 0] > obs_rates[4, 0]
    assert out[5, 2] == max(fb_rates[5, 2], out[4, 2])
    # Opportunity column never takes the fallback.
    assert out[5, 1] == max(obs_rates[5, 1], out[4, 1])
    # C group (all observed) untouched by the A-group fallback.
    assert np.allclose(out[:, 6:9], np.maximum.accumulate(obs_rates[:, 6:9]))


def test_hp_site_total_is_per_context():
    fb = np.full([10, 12], 10.0)
    fb[:, 1::3] = 1e6
    fb_rates = _hp_rates(fb)
    obs = np.ones([10, 12])
    obs[:, 1::3] = 1e5
    # HP5 A group alone is < 1000 sites though the full row is not: A's
    # unobserved del/ins still fall back.
    obs[4, 0:3] = [0, 500, 0]
    # HP5 T group is well sampled: its unobserved del does not.
    obs[4, 3:6] = [0, 5000, 1]
    out = _normalize_indel_hp_mat(obs, 0.5, fb)
    plain = _normalize_indel_hp_mat(obs, 0.5)
    assert out[4, 0] == max(fb_rates[4, 0], out[3, 0])
    assert out[4, 2] == max(fb_rates[4, 2], out[3, 2])
    assert np.array_equal(out[:, 3:], plain[:, 3:])


def test_hp_zero_base_group_in_covered_row_takes_fallback():
    # HP9 C group has no observations at all while A/T/G are well sampled:
    # it must take the fallback, not the flat 1/3 pseudocount rate.
    fb = np.full([10, 12], 10.0)
    fb[:, 1::3] = 1e6
    fb_rates = _hp_rates(fb)
    obs = np.full([10, 12], 2.0)
    obs[:, 1::3] = 1e4
    obs[8, 6:9] = 0
    out = _normalize_indel_hp_mat(obs, 0.5, fb)
    assert out[8, 6] == max(fb_rates[8, 6], out[7, 6])
    assert out[8, 8] == max(fb_rates[8, 8], out[7, 8])
    assert out[9, 6] < 1e-3 and out[9, 8] < 1e-3


def test_hp_first_row_takes_fallback_as_is():
    fb = np.full([10, 12], 10.0)
    fb[:, 1::3] = 1e6
    obs = np.zeros([10, 12])
    obs[:, 1::3] = 100
    out = _normalize_indel_hp_mat(obs, 0.5, fb)
    fb_rates, obs_rates = _hp_rates(fb), _hp_rates(obs)
    assert np.array_equal(out[0, 0::3], fb_rates[0, 0::3])
    assert np.array_equal(out[0, 2::3], fb_rates[0, 2::3])
    assert np.array_equal(out[0, 1::3], obs_rates[0, 1::3])


def test_str_fallback_rows():
    pc = 0.5
    fb = np.full([5, 11], 5.0)
    fb[:, 5] = 1e6
    fb[3, 3] = 5000  # 25-39 bin: high -2 rate
    obs = np.full([5, 11], 1.0)
    obs[:, 5] = 1e5
    obs[1, :] = 0
    obs[1, 5] = 10  # row 1 low, all indels unobserved -> fallback as-is
    obs[3, :] = 0
    obs[3, 5] = 10
    obs[3, 7] = 2  # row 3 low, +2 observed
    obs[4, :] = 0  # row 4 low, compared to row 3's used rate
    out = _normalize_indel_str_mat(obs, pc, fb)
    fb_rates = (fb + pc) / (fb.sum(axis=1, keepdims=True) + 11 * pc)
    obs_rates = (obs + pc) / (obs.sum(axis=1, keepdims=True) + 11 * pc)
    idx = np.arange(11) != 5
    assert np.array_equal(out[0], obs_rates[0])
    assert np.array_equal(out[1, idx], fb_rates[1, idx])
    assert out[1, 5] == obs_rates[1, 5]
    assert np.array_equal(out[2], obs_rates[2])
    zero3 = idx & (np.arange(11) != 7)
    assert np.array_equal(out[3, zero3], np.maximum(fb_rates[3, zero3], out[2, zero3]))
    assert out[3, 7] == obs_rates[3, 7]
    assert np.array_equal(out[4, idx], np.maximum(fb_rates[4, idx], out[3, idx]))
    assert out[4, 3] == fb_rates[3, 3]


def test_indel_without_fallback_unchanged_behavior():
    obs = np.zeros([5, 11])
    obs[:2, 5] = 1e5
    out = _normalize_indel_str_mat(obs, 0.5)
    assert np.array_equal(out[2], out[1]) and np.array_equal(out[4], out[1])
