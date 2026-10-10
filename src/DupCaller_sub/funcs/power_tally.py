"""Round-1 power tallies for call --skipCoveragePass.

The coverage pass (round 2) weights every eligible family-position by a
detection-power lookup -- L[n_top, n_bot, trinuc] for SBS, L_indel_1bp /
L_indel_len for indels -- at the channels' final thresholds, which only
exist once round 1 has finished. Every by-duplex-group total it writes is a
sum of such lookups, so round 1 can count positions by power-table index
instead (PowerTally) and Caller.py contracts the counts with power tables
built at the final thresholds (contract_power_tally). The result equals
round 2's totals up to float summation order; only the per-locus
coverage.bed.gz needs the per-position pass.

Key spaces (d = n_top * 10 + n_bot, both capped at 9):
  sbs    d, trinuc                              -> L[d, t, :]
  dbs    d, trinuc of the first base, last base of the second base's
         trinuc (adjacent positions with the same d); other adjacent
         pairs are kept sparse as (sbs key 1, sbs key 2)
  hp     L_indel_1bp flat index (d, run 0-11, ref base, sign)
  strdel / strins  d, STR bin, unit bucket, repeat-count bucket
  rep1   d, STR bin, unit bucket (repeat count 1 STR insertions)
  anti   d (antimask_indel positions: flat novel-indel opportunity)
  ins1   d, next reference base (1bp novel-insertion opportunity)
"""

import numpy as np

from .misc import _DBS_RAW144_IDX_GRID, build_trinuc64_order

N_DEPTH = 100
SBS_SIZE = N_DEPTH * 64
DBS_SIZE = SBS_SIZE * 4
HP_SIZE = N_DEPTH * 12 * 4 * 2
STR_SIZE = N_DEPTH * 5 * 4 * 6
REP1_SIZE = N_DEPTH * 5 * 4
ANTI_SIZE = N_DEPTH
INS1_SIZE = N_DEPTH * 4

_SIZES = {
    "sbs": SBS_SIZE,
    "dbs": DBS_SIZE,
    "hp": HP_SIZE,
    "strdel": STR_SIZE,
    "strins": STR_SIZE,
    "rep1": REP1_SIZE,
    "anti": ANTI_SIZE,
    "ins1": INS1_SIZE,
}
OFFSET = {}
_off = 0
for _name, _size in _SIZES.items():
    OFFSET[_name] = _off
    _off += _size
TOTAL_SIZE = _off

_BASES = "ATCG"
_TRINUC2NUM, _NUM2TRINUC = build_trinuc64_order()
# Middle and last base (A/T/C/G = 0-3) of each trinuc, and the trinuc one
# position to the right given the next base: SHIFT[t, b].
MID_BASE = np.array([_BASES.index(t[1]) for t in _NUM2TRINUC])
LAST_BASE = np.array([_BASES.index(t[2]) for t in _NUM2TRINUC])
SHIFT = np.array(
    [[_TRINUC2NUM[t[1:] + b] for b in _BASES] for t in _NUM2TRINUC], dtype=np.int64
)
_TO_CGTA = np.array([3, 2, 0, 1])  # base2num -> CGTA order, as call.py


def sbs_power_index(F1R2_count, F2R1_count, trinuc, valid):
    """Per-position flat index (d * 64 + trinuc) of L's first three axes,
    -1 where genotypeDSSnv leaves cov_mat at 0 (valid False or trinuc >= 64).
    """
    n_top = np.minimum(F1R2_count.sum(axis=0), 9)
    n_bot = np.minimum(F2R1_count.sum(axis=0), 9)
    ok = valid & (trinuc < 64)
    return np.where(ok, (n_top * 10 + n_bot) * 64 + trinuc, -1).astype(np.int64)


class PowerTally:
    """Per-worker, per-duplex-group counts of coverage opportunity by
    power-table index. Keys are buffered per group and folded into one
    dense int64 vector (TOTAL_SIZE) per group every FLUSH keys."""

    FLUSH = 1 << 16

    def __init__(self):
        self.counts = {}
        self._buf = {}
        self._buf_n = {}
        # Sparse DBS pairs: per group, pending key arrays and their
        # compacted (keys, counts), merged every FLUSH pending keys.
        self.dbs_other = {}
        self._dbs_n = {}
        self._dbs_compact = {}
        self.sbs_unmasked = np.zeros(SBS_SIZE, dtype=np.int64)

    def _add(self, group, keys):
        if keys.size == 0:
            return
        self._buf.setdefault(group, []).append(keys)
        n = self._buf_n.get(group, 0) + keys.size
        self._buf_n[group] = n
        if n >= self.FLUSH:
            self._flush(group)

    def _flush(self, group):
        buf = self._buf.pop(group, None)
        self._buf_n.pop(group, None)
        if not buf:
            return
        if group not in self.counts:
            self.counts[group] = np.zeros(TOTAL_SIZE, dtype=np.int64)
        self.counts[group] += np.bincount(np.concatenate(buf), minlength=TOTAL_SIZE)

    def _compact_dbs(self, group):
        arrays = self.dbs_other.pop(group, None)
        self._dbs_n.pop(group, None)
        if not arrays:
            return
        keys = np.concatenate(arrays)
        counts = np.ones(keys.size, dtype=np.int64)
        if group in self._dbs_compact:
            old_keys, old_counts = self._dbs_compact[group]
            keys = np.concatenate([old_keys, keys])
            counts = np.concatenate([old_counts, counts])
        uniq, inv = np.unique(keys, return_inverse=True)
        self._dbs_compact[group] = (
            uniq,
            np.bincount(inv, weights=counts).astype(np.int64),
        )

    def add_sbs(self, group, pidx, pass_bool, unmasked_pass_bool, ref_int):
        """SBS opportunity of one family's window, as call.py's coverage /
        duplex_read_num_dict_trinuc / unmasked_coverage /
        _compute_dbs_opportunity accumulation: pidx from sbs_power_index."""
        p = pidx[pass_bool]
        self._add(group, OFFSET["sbs"] + p[p >= 0])
        u = pidx[unmasked_pass_bool]
        u = u[u >= 0]
        if u.size:
            self.sbs_unmasked += np.bincount(u, minlength=SBS_SIZE)
        if pidx.size < 2:
            return
        i1, i2 = pidx[:-1], pidx[1:]
        pair = (
            pass_bool[:-1]
            & pass_bool[1:]
            & (ref_int[:-1] <= 3)
            & (ref_int[1:] <= 3)
            & (i1 >= 0)
            & (i2 >= 0)
        )
        if not pair.any():
            return
        i1, i2 = i1[pair], i2[pair]
        d1, t1 = i1 // 64, i1 % 64
        d2, t2 = i2 // 64, i2 % 64
        last2 = LAST_BASE[t2]
        dense = (d1 == d2) & (SHIFT[t1, last2] == t2)
        self._add(group, OFFSET["dbs"] + i1[dense] * 4 + last2[dense])
        if not dense.all():
            other = i1[~dense] * SBS_SIZE + i2[~dense]
            self.dbs_other.setdefault(group, []).append(other)
            n = self._dbs_n.get(group, 0) + other.size
            self._dbs_n[group] = n
            if n >= self.FLUSH:
                self._compact_dbs(group)

    def add_indel(self, group, keys):
        """Keys from indel_keys for one family's window."""
        for name, k in keys:
            self._add(group, OFFSET[name] + k)

    def summary(self, all_quals):
        """Flush everything; return a picklable sparse summary for
        Caller.py: per group (nonzero indices, counts) of the dense vector
        and the (keys, counts) of its sparse DBS pairs."""
        for group in list(self._buf):
            self._flush(group)
        groups = {}
        for group, vec in self.counts.items():
            nz = np.flatnonzero(vec)
            groups[group] = (nz, vec[nz])
        for group in list(self.dbs_other):
            self._compact_dbs(group)
        dbs_other = dict(self._dbs_compact)
        return {
            "quals": all_quals,
            "groups": groups,
            "dbs_other": dbs_other,
            "sbs_unmasked": self.sbs_unmasked,
        }


def indel_keys(
    n_top_indel,
    n_bot_indel,
    hp_run_power,
    ref_allele_safe,
    hp_run_capped10,
    hp_repeat_valid,
    hp_cut_bool,
    is_real_str_row,
    str_repeat_valid,
    str_cut_bool,
    str_bin_row,
    str_unit_row,
    str_repeat_row,
    antimask_indel,
    next_ref_arr,
    next_ref_valid,
):
    """Tally keys for every term of call.py's indel100 accumulation for one
    family's window (same gates, same power-table indices)."""
    d = n_top_indel.astype(np.int64) * 10 + n_bot_indel
    hp_idx = ((d * 12 + hp_run_power) * 4 + ref_allele_safe) * 2
    hp_gate = hp_repeat_valid & hp_cut_bool
    ins_gate = (hp_run_capped10 >= 1) & hp_gate
    unit_bucket = np.clip(str_unit_row, 2, 5) - 2
    str_gate = is_real_str_row & str_repeat_valid & str_cut_bool
    str_base = (d * 5 + str_bin_row) * 4 + unit_bucket
    count_del = np.clip(str_repeat_row, 1, 6) - 1
    count_ins = np.clip(str_repeat_row, 0, 5)
    rep1_gate = str_gate & (str_repeat_row == 1)
    ins1_gate = antimask_indel & next_ref_valid
    return (
        ("hp", hp_idx[hp_gate]),
        ("hp", hp_idx[ins_gate] + 1),
        ("strdel", (str_base * 6 + count_del)[str_gate]),
        ("strins", (str_base * 6 + count_ins)[str_gate]),
        ("rep1", str_base[rep1_gate]),
        ("anti", d[antimask_indel]),
        ("ins1", (d * 4 + next_ref_arr)[ins1_gate]),
    )


def _hp_columns():
    """indel100 column of every L_indel_1bp flat index (per depth)."""
    cols = np.zeros((12, 4, 2), dtype=np.int64)
    for h in range(12):
        hc10 = min(h, 10)
        for r in range(4):
            base = _TO_CGTA[r]
            cols[h, r, 0] = base * 6 + np.clip(hc10, 1, 6) - 1
            cols[h, r, 1] = 24 + base * 5 + np.clip(hc10, 1, 5) - 1
    return np.tile(cols.ravel(), N_DEPTH)


_HP_COLS = _hp_columns()


def _dense(idx, counts, name):
    sel = (idx >= OFFSET[name]) & (idx < OFFSET[name] + _SIZES[name])
    return idx[sel] - OFFSET[name], counts[sel].astype(float)


def _dbs_contribution(i1, i2, counts, L3):
    """Raw-144 DBS opportunity of adjacent position pairs with SBS keys
    i1, i2 (count each), as _compute_dbs_opportunity."""
    out = np.zeros(144)
    if i1.size == 0:
        return out
    c1 = L3[i1]
    c2 = L3[i2]
    pair = c1[:, :, None] * c2[:, None, :] * counts[:, None, None]
    raw_idx = _DBS_RAW144_IDX_GRID[MID_BASE[i1 % 64], MID_BASE[i2 % 64]]
    flat = raw_idx.reshape(-1)
    keep = flat >= 0
    out += np.bincount(flat[keep], weights=pair.reshape(-1)[keep], minlength=144)
    return out


def contract_power_tally(summary, L, L_indel_1bp, L_indel_len):
    """One worker's tallies times its power tables: per duplex group the
    (64, 4) trinuc, (100,) indel100 and (144,) DBS opportunity call.py's
    coverage pass accumulates, plus the masked and unmasked SBS coverage
    totals."""
    L3 = L.reshape(SBS_SIZE, 4)
    Llen = L_indel_len.reshape(N_DEPTH, 5, 11)
    trinuc, indel, dbs = {}, {}, {}
    coverage = 0.0
    for group, (idx, cnt) in summary["groups"].items():
        k, c = _dense(idx, cnt, "sbs")
        t = np.zeros((SBS_SIZE, 4))
        t[k] = L3[k] * c[:, None]
        trinuc[group] = t.reshape(N_DEPTH, 64, 4).sum(axis=0)
        coverage += trinuc[group].sum()

        k, c = _dense(idx, cnt, "dbs")
        i1 = k // 4
        i2 = (i1 // 64) * 64 + SHIFT[i1 % 64, k % 4]
        dbs[group] = _dbs_contribution(i1, i2, c, L3)
        if group in summary["dbs_other"]:
            keys, kc = summary["dbs_other"][group]
            dbs[group] += _dbs_contribution(
                keys // SBS_SIZE, keys % SBS_SIZE, kc.astype(float), L3
            )

        v = np.zeros(100)
        k, c = _dense(idx, cnt, "hp")
        np.add.at(v, _HP_COLS[k], c * L_indel_1bp.reshape(-1)[k])
        for name, col0, sign in (("strdel", 44, -1), ("strins", 68, 1)):
            k, c = _dense(idx, cnt, name)
            cb = k % 6
            rest = k // 6
            ub = rest % 4
            sb = (rest // 4) % 5
            dd = rest // 20
            w = Llen[dd, sb, sign * (ub + 2) + 5]
            np.add.at(v, col0 + ub * 6 + cb, c * w)
        k, c = _dense(idx, cnt, "rep1")
        ub = k % 4
        sb = (k // 4) % 5
        dd = k // 20
        np.add.at(v, 68 + ub * 6, c * Llen[dd, sb, ub + 2 + 5])
        k, c = _dense(idx, cnt, "anti")
        for u in range(4):
            U = u + 2
            v[44 + u * 6] += np.sum(c * Llen[k, 0, -U + 5])
            ins = np.sum(c * Llen[k, 0, U + 5])
            v[68 + u * 6] += ins
            v[68 + u * 6 + 1] += ins
        for j in range(4):
            v[92 + j] += np.sum(c * Llen[k, 0, 3 - j])
        k, c = _dense(idx, cnt, "ins1")
        dd, nxt = k // 4, k % 4
        w = c * Llen[dd, 0, 6]
        for b in range(4):
            v[96 + b] += np.sum(w[nxt != b])
        indel[group] = v
    unmasked = float((summary["sbs_unmasked"][:, None] * L3).sum())
    return trinuc, indel, dbs, coverage, unmasked
