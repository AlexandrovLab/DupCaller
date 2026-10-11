import math
import os
import re
import sys
import time
from multiprocessing import Pool
import numpy as np
import h5py
import pysam
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
import pandas as pd
from scipy.optimize import brentq


from .prob import calculateSSPosterior, calculateDSPosterior
from .indels import STR_NCOL, STR_UNIT_GROUPS


# ======================================================================
# Generic / file / progress utilities
# ======================================================================
def _ensure_type_subdirs(prefix):
    """Create (if missing) the SBS/, INDEL/, DBS/ output subfolders under
    prefix. Type-specific outputs go in these; whole-sample files
    (_stats.txt, _coverage.bed.gz, ERROR/, ...) stay at the top level.
    Shared by Caller.py and Estimate.py so both use the same layout.

    Returns (sbs_dir, indel_dir, dbs_dir).
    """
    sbs_dir = os.path.join(prefix, "SBS")
    indel_dir = os.path.join(prefix, "INDEL")
    dbs_dir = os.path.join(prefix, "DBS")
    for d in (sbs_dir, indel_dir, dbs_dir):
        os.makedirs(d, exist_ok=True)
    return sbs_dir, indel_dir, dbs_dir


def log_progress(state, label, index, total=None, extra="", interval=5.0):
    """Print a progress line (count, percent, elapsed) at most once every
    `interval` seconds, plus always on the final item when `total` is known.

    state: {"start": <loop start time>, "last": <time of last print>},
        kept by the caller across the loop.
    total: item count if known; None prints a running count; 0 is a no-op.
    extra: optional short string appended to the line.
    """
    if total == 0:
        return
    now = time.time()
    is_last = total is not None and index + 1 == total
    if not is_last and now - state["last"] < interval:
        return
    if total is not None:
        pct = (index + 1) / total * 100
        counter = f"{index + 1}/{total} ({pct:5.1f}%)"
    else:
        counter = f"{index + 1} processed"
    suffix = f"  {extra}" if extra else ""
    print(f"{label}: {counter}{suffix}  elapsed {now - state['start']:.1f}s")
    state["last"] = now


# Labels for the 16 power-weighted indel coverage categories, in the order
# of the trailing indel columns of the coverage bed file (after the 4 SBS
# columns) and the "Effective/Unmasked Indel Coverage (<label>)" lines in
# _stats.txt. Shared by Caller.py (writer) and Estimate.py/Summarize.py
# (readers).
#
# Deletion/Insertion Length 2-5+ use a fixed no-repeat context, not the
# position's actual repeat context.
#
# Insertion A/T/C/G: rep0 opportunity to insert that base (fixed hps=1
# context), 0 where the following reference base already equals it (that
# insertion would extend a run, not start one). Independent of the
# position's own base.
INDEL_COVERAGE_CATEGORY_LABELS = [
    "Deletion of Repeat Unit",
    "Insertion of Repeat Unit",
    "Deletion Length 2",
    "Deletion Length 3",
    "Deletion Length 4",
    "Deletion Length 5+",
    "Insertion A",
    "Insertion T",
    "Insertion C",
    "Insertion G",
    "Insertion Length 2",
    "Insertion Length 3",
    "Insertion Length 4",
    "Insertion Length 5+",
    # Columns 0/1 ("...of Repeat Unit") are STR-only: one copy of the
    # position's annotated STR unit, 0 where unit_len < 2. The +-1bp
    # homopolymer-slip opportunity (always scored via the HP table, even
    # inside STRs; see prob.py's indelErrorProbs) lives in these last two
    # columns, from the self-derived homopolymer run, independent of STR
    # status. Keeping the columns non-overlapping lets them be summed
    # (e.g. indel_cov_total) without double-counting.
    "Deletion of Homopolymer",
    "Insertion of Homopolymer",
]


def parse_stats_file(path):
    """Parse a DupCaller _stats.txt file into a {label: value_string} dict,
    keyed by the text before the first tab on each line (so readers don't
    depend on line order)."""
    stats = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if "\t" not in line:
                continue
            label, value = line.split("\t", 1)
            stats[label] = value
    return stats


# Bump whenever the on-disk layout of .ref.h5/.tn.h5/.hp.h5/.str.h5 changes
# in a way ndim/shape checks can't detect. Current layout: hp.h5 has 2 rows
# (hp, cut; self-derived from the reference, homopolymers only); str.h5 has
# 3 rows (unit_len, hp, cut; from --repeatBed/--strbed, unit_len >= 2, with
# overlapping tracts resolved so each flagged start carries its own tract's
# unit/count). Index.py stamps this onto each h5 file's root attrs;
# check_h5_usable verifies it.
INDEX_SCHEMA_VERSION = 3


def check_h5_usable(h5_path, expected_ndim=None):
    """Try to open and read from an h5 file; return (ok, error_message).

    A missing or mismatched index_schema_version attribute means the index
    has a different on-disk layout, which ndim/shape checks can't catch.
    """
    try:
        with h5py.File(h5_path, "r") as f:
            keys = list(f.keys())
            if not keys:
                return False, f"{h5_path}: file contains no chromosome datasets"
            version = f.attrs.get("index_schema_version")
            if version is None:
                return False, (
                    f"{h5_path}: no index_schema_version attribute found — "
                    "this index is outdated (built by an old version of "
                    "`DupCaller.py index`)"
                )
            if version != INDEX_SCHEMA_VERSION:
                return False, (
                    f"{h5_path}: index_schema_version={version}, but this "
                    f"version of DupCaller expects {INDEX_SCHEMA_VERSION} — "
                    "the index is outdated"
                )
            data = f[keys[0]]
            ndim = len(data.shape)
            if expected_ndim is not None and ndim != expected_ndim:
                return False, (
                    f"{h5_path}: expected {expected_ndim}D data, found {ndim}D — "
                    "file may be from an incompatible version of DupCaller"
                )
            _ = data[tuple(slice(0, min(s, 10)) for s in data.shape)]
    except Exception as e:
        return False, f"{h5_path}: could not read h5 file ({e})"
    return True, None


# ======================================================================
# SBS / trinucleotide classification
# ======================================================================
def build_trinuc64_order():
    """Canonical order of the 64 real trinucleotide contexts (32 C/T-ref
    followed by 32 G/A-ref), matching Index.py's .tn.h5 trinuc_int encoding
    for indices 0-63. Index i (i<32) and i+32 are reverse complements of
    each other. Returns (trinuc2num, num2trinuc). Shared by every writer and
    reader of trinuc-indexed data.
    """
    num2trinuc = []
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                num2trinuc.append(minus_base + ref_base + plus_base)
    for plus_base in ["T", "A", "G", "C"]:
        for ref_base in ["G", "A"]:
            for minus_base in ["T", "A", "G", "C"]:
                num2trinuc.append(minus_base + ref_base + plus_base)
    trinuc2num = {t: i for i, t in enumerate(num2trinuc)}
    return trinuc2num, num2trinuc


def build_trinuc192_labels(num2trinuc):
    """192 raw (unfolded) mutation-class labels: for each of the 64
    trinucleotide contexts, its 3 non-reference alt bases in "ATCG" order,
    formatted "{trinuc}>{alt}", no reverse-complement folding. Returns
    (label2num, labels); this order is the contract between Caller.py
    (writer of _trinuc_by_duplex_group.txt) and Estimate.py (reader).
    """
    labels = []
    label2num = {}
    for trinuc in num2trinuc:
        ref_base = trinuc[1]
        for alt in "ATCG":
            if alt == ref_base:
                continue
            label = f"{trinuc}>{alt}"
            label2num[label] = len(labels)
            labels.append(label)
    return label2num, labels


_REVCOMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _build_sbs96_labels():
    """96 canonical (C/T-ref, pyrimidine-convention) SBS mutation-class
    labels, bracket notation "{minus}[{ref}>{alt}]{plus}"."""
    label2num = {}
    labels = []
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                alts = ["A", "T", "C", "G"]
                alts.remove(ref_base)
                for alt_base in alts:
                    label = f"{minus_base}[{ref_base}>{alt_base}]{plus_base}"
                    label2num[label] = len(labels)
                    labels.append(label)
    return label2num, labels


TRINUCSBS2NUM_96, NUM2TRINUCSBS_96 = _build_sbs96_labels()


def combine_raw192_to_sbs96(values_192, label2num_192):
    """Combine the 192 raw (un-folded) (trinuc, alt) mutation classes into
    the 96 canonical pyrimidine-convention SBS96 classes, by summing each
    canonical class with its reverse-complement partner.

    values_192   : array shaped (192,) or (192, n_groups).
    label2num_192: {"{trinuc}>{alt}": row_index} for values_192's rows.
    Returns an array shaped (96,) or (96, n_groups).
    """
    out_shape = (96,) if values_192.ndim == 1 else (96, values_192.shape[1])
    combined = np.zeros(out_shape)
    for k, label in enumerate(NUM2TRINUCSBS_96):
        trinuc = label[0] + label[2] + label[6]
        alt = label[4]
        rc_trinuc = _REVCOMP[trinuc[2]] + _REVCOMP[trinuc[1]] + _REVCOMP[trinuc[0]]
        rc_alt = _REVCOMP[alt]
        combined[k] = (
            values_192[label2num_192[f"{trinuc}>{alt}"]]
            + values_192[label2num_192[f"{rc_trinuc}>{rc_alt}"]]
        )
    return combined


def _safe_rate(mut, cov):
    """np.where(cov > 0, mut / cov, 0) without the divide-by-zero
    RuntimeWarning (np.where evaluates both branches)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(cov > 0, mut / cov, 0)


# ======================================================================
# DBS classification
# ======================================================================
def _num_dinuc(n):
    """Inverse of the "4*base2num[dinuc[0]]+base2num[dinuc[1]]" dinucleotide
    encoding (A=0,T=1,C=2,G=3, same convention as ref.h5)."""
    n2b = "ATCG"
    return n2b[n // 4] + n2b[n % 4]


_DBS_BASES = ["A", "T", "C", "G"]


def _dbs_alt_choices(ref):
    """The 9 possible alt dinucleotides for a reference dinucleotide (both
    positions must change)."""
    return [
        b1 + b2
        for b1 in _DBS_BASES
        if b1 != ref[0]
        for b2 in _DBS_BASES
        if b2 != ref[1]
    ]


def build_dbs_raw144_labels():
    """The 144 raw (un-folded) "REF>ALT" DBS labels (16 reference
    dinucleotides x 9 alts), the DBS analogue of build_trinuc192_labels;
    folded to DBS78 by Estimate.py's combine_raw_dbs_to_dbs78. Shared by
    the writer (_dbs_by_duplex_group.txt) and Estimate.py.
    """
    label2num = {}
    labels = []
    for ref_num in range(16):
        ref = _num_dinuc(ref_num)
        for alt in _dbs_alt_choices(ref):
            label = f"{ref}>{alt}"
            label2num[label] = len(labels)
            labels.append(label)
    assert len(labels) == 144
    return label2num, labels


def build_dbs_raw144_index_grid():
    """[4,4,4,4] lookup grid mapping (ref1, ref2, alt1, alt2), "ATCG"-order
    0-3 encoding, to a build_dbs_raw144_labels row index, for vectorized
    accumulation. Entries with alt1==ref1 or alt2==ref2 (not a DBS) are -1,
    which the caller filters out before np.bincount.
    """
    label2num, _ = build_dbs_raw144_labels()
    n2b = "ATCG"
    grid = np.full((4, 4, 4, 4), -1, dtype=np.int64)
    for r1 in range(4):
        for r2 in range(4):
            ref = n2b[r1] + n2b[r2]
            for a1 in range(4):
                if a1 == r1:
                    continue
                for a2 in range(4):
                    if a2 == r2:
                        continue
                    alt = n2b[a1] + n2b[a2]
                    grid[r1, r2, a1, a2] = label2num[f"{ref}>{alt}"]
    return grid


_DBS_RAW144_IDX_GRID = build_dbs_raw144_index_grid()


def _detect_dbs_pairs(
    mut_chrom,
    mut_positions,
    muts_ind,
    ref_int,
    alt_int,
    pass_bool,
    LR_pass_bool,
    flt_rs,
    F1R2,
    F2R1,
    setBc,
    currentStart,
    template_length,
    num2base,
    reads,
):
    """Detect DBS events within one duplex family's SNV candidates: two
    adjacent positions changed in the same molecule (the inputs all come
    from this family's own consensus calling). A pair is called only if
    both SNVs are PASS, decided from pass_bool/LR_pass_bool/flt_rs the
    same way as the caller's per-position SNV filter.

    reads: the family's reads; their fragment starts (FS) go on every DBS
    record for the depth-extraction founding-family exemption.

    Returns a list of DBS mut dicts (2-base ref/alt), shaped like the SNV
    mut dicts minus the SNV-only INFO fields.
    """

    def _flt(i):
        if pass_bool[i] and LR_pass_bool[i]:
            return flt_rs
        return "masked"

    dbs_muts = []
    fragment_starts = None
    for k in range(len(mut_positions) - 1):
        if mut_positions[k + 1] != mut_positions[k] + 1:
            continue
        i0, i1 = muts_ind[k], muts_ind[k + 1]
        if _flt(i0) != "PASS" or _flt(i1) != "PASS":
            continue
        if fragment_starts is None:
            fragment_starts = family_fragment_starts(reads)
        ref_dinuc = num2base[ref_int[i0]] + num2base[ref_int[i1]]
        alt_dinuc = num2base[alt_int[i0]] + num2base[alt_int[i1]]
        dbs_muts.append(
            {
                "chrom": mut_chrom,
                "pos": mut_positions[k],
                "ref": ref_dinuc,
                "alt": alt_dinuc,
                "filter": "PASS",
                "infos": {
                    "F1R2": F1R2,
                    "F2R1": F2R1,
                    "TAG1": setBc[0],
                    "TAG2": setBc[1],
                    "SP": currentStart,
                    "TL": template_length,
                    "FS": fragment_starts,
                },
                "formats": ["AC", "RC", "DP"],
                # Placeholder AC/RC/DP (tumor 1/0/1, normal 0/0/0): DBS has
                # no raw-BAM depth extraction like SNVs/indels do.
                "samples": [[1, 0, 1], [0, 0, 0]],
            }
        )
    return dbs_muts


def _compute_dbs_opportunity(cov_mat, ref_int, pass_bool):
    """DBS opportunity from one family's window: for each adjacent pair
    (i, i+1) with both positions in pass_bool and ACGT reference bases,
    add coverage(alt1 at i) * coverage(alt2 at i+1) for each of the 9
    non-ref (alt1, alt2) combinations. Same product-of-coverages
    approximation as Estimate.py's calculate_dbs_opportunity, but per
    family so it can be attributed to its duplex group.

    cov_mat  : (window_len, 4) L-weighted per-alt-base coverage.
    ref_int  : (window_len,) reference base, 0-3=A/T/C/G, 4=N/invalid.
    pass_bool: (window_len,) bool, calling-eligible positions.
    Returns a (144,) vector in build_dbs_raw144_labels order.
    """
    contribution = np.zeros(144)
    if cov_mat.shape[0] < 2:
        return contribution
    pair_ok = pass_bool[:-1] & pass_bool[1:]
    if not np.any(pair_ok):
        return contribution
    r1 = ref_int[:-1]
    r2 = ref_int[1:]
    valid = pair_ok & (r1 <= 3) & (r2 <= 3)
    if not np.any(valid):
        return contribution
    c1 = cov_mat[:-1][valid]  # (k, 4)
    c2 = cov_mat[1:][valid]  # (k, 4)
    pair_contrib = c1[:, :, None] * c2[:, None, :]  # (k, 4, 4)
    raw_idx = _DBS_RAW144_IDX_GRID[r1[valid], r2[valid]]  # (k, 4, 4)
    flat_idx = raw_idx.reshape(-1)
    keep = flat_idx >= 0
    contribution += np.bincount(
        flat_idx[keep], weights=pair_contrib.reshape(-1)[keep], minlength=144
    )
    return contribution


# ======================================================================
# Indel classification (indel100 / ID83)
# ======================================================================
# Bucket boundaries for the 100-class raw indel scheme (build_indel100_labels).
# Deletion bins start at 1 (there must be a run to delete from); insertion
# STR bins start at 0 (insertion into a non-repeat context). For STR
# insertion, repeat_count 1 is credited to both bucket "0" and bucket "1".
INDEL_DEL_HP_LEN_BINS = ["1", "2", "3", "4", "5", "6+"]
# ID83-resolution homopolymer-insertion bins (build_indel83_labels only),
# including the rep0 "0" bin (Cinshp0/Tinshp0).
INDEL_INS_HP_LEN_BINS = ["0", "1", "2", "3", "4", "5+"]
# Raw100 homopolymer-insertion bins: no "0" bin, since rep0 opportunity
# comes from the 1bpins{base} columns instead (Estimate.py's
# override_inshp0_with_next_base_opportunity).
INDEL100_INS_HP_LEN_BINS = ["1", "2", "3", "4", "5+"]
INDEL_STR_UNIT_BINS = ["2", "3", "4", "5+"]
INDEL_DEL_STR_COUNT_BINS = ["1", "2", "3", "4", "5", "6+"]
INDEL_INS_STR_COUNT_BINS = ["0", "1", "2", "3", "4", "5+"]
INDEL_MH_DEL_LEN_BINS = ["2", "3", "4", "5+"]


def build_indel100_labels():
    """100 raw indel classification labels, in order:
      - homopolymer deletion, 4 bases x 6 length bins (24)
      - homopolymer insertion, 4 bases x 5 length bins, no rep0 bin (20)
      - STR deletion / insertion, 4 unit sizes x 6 count bins (24 + 24)
      - microhomology deletion length (4)
      - 1bpins{base} (4): rep0 opportunity to insert that base, i.e. every
        position whose following reference base differs from it; feeds the
        ID83 Cinshp0/Tinshp0 bins (Estimate.py's
        override_inshp0_with_next_base_opportunity).
    Returns (label2num, labels); shared by writer and readers.
    """
    labels = []
    label2num = {}

    def add(label):
        label2num[label] = len(labels)
        labels.append(label)

    for base in "CGTA":
        for length in INDEL_DEL_HP_LEN_BINS:
            add(f"del{base}_hp{length}")
    for base in "CGTA":
        for length in INDEL100_INS_HP_LEN_BINS:
            add(f"ins{base}_hp{length}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_DEL_STR_COUNT_BINS:
            add(f"delSTR{unit}_rep{count}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_INS_STR_COUNT_BINS:
            add(f"insSTR{unit}_rep{count}")
    for length in INDEL_MH_DEL_LEN_BINS:
        add(f"delMH_len{length}")
    for base in "ATCG":
        add(f"1bpins{base}")

    return label2num, labels


# base2num (A=0,T=1,C=2,G=3) -> column order used by build_indel100_labels
# ("CGTA"), shared by every indel100 classification function below.
_INDEL_BASE2NUM_TO_CGTA = np.array([3, 2, 0, 1])


def indel100_reference_bucket_indices(
    hp_run_arr,
    str_unit_len_arr,
    str_repeat_count_arr,
    ref_base_arr,
    next_ref_base_arr=None,
    hp_cut_arr=None,
    str_cut_arr=None,
):
    """Genome-wide indel100 opportunity counts (build_indel100_labels order)
    from whole-chromosome index arrays: hp_run_arr = hp.h5 row 0
    (self-derived homopolymer run length), str_unit_len_arr/
    str_repeat_count_arr = str.h5 rows 0/1 (unit_len 0 where no STR is
    annotated), ref_base_arr = ref.h5 (A=0,T=1,C=2,G=3). Each eligible
    position adds 1 to each column it qualifies for.

    HP and STR opportunity are credited independently (unlike
    load_repeat_context, where STR wins): an embedded homopolymer inside an
    STR counts toward both. Every position also gets a flat, scan-free
    rep1 (deletion) / rep0+rep1 (insertion) credit for each STR unit size,
    matching call.py's fixed-context coverage credit.

    hp_cut_arr, str_cut_arr (optional): hp.h5 row 1 / str.h5 row 2
    (start-of-run flags). A repeat tract is one locus (slippage is reported
    at one left-aligned anchor), so these restrict HP/STR credit to one per
    tract; otherwise long tracts, which are also the worst covered, would
    be over-weighted. If omitted, every position is credited.

    next_ref_base_arr (optional): reference base following each position
    (out-of-range sentinel, e.g. -1, where there is none). Feeds the
    1bpins{base} columns 96-99 (A/T/C/G): 96+b counts positions whose next
    base differs from b. If omitted, those columns stay 0.

    Anchor convention (shared with call.py's indel100 accumulation):
    position P is opportunity for an event whose VCF anchor is P-1, so a
    VCF's 1-based rec.pos equals P in 0-based terms.

    Returns the (100,) count array. Contributions are accumulated in place
    rather than concatenated, to avoid many chromosome-length temporaries.
    """
    n = hp_run_arr.shape[0]
    # Non-ACGT reference positions are never opportunity (the coverage side
    # masks them too), so every column is restricted to valid_base and the
    # flat credits use n_valid.
    valid_base = (ref_base_arr >= 0) & (ref_base_arr <= 3)
    n_valid = int(np.count_nonzero(valid_base))
    ref_base_safe = np.where(valid_base, ref_base_arr, 0)
    real_str = (str_unit_len_arr >= 2) & valid_base
    base_cgta = _INDEL_BASE2NUM_TO_CGTA[ref_base_safe]
    hp_cut = (
        np.asarray(hp_cut_arr, dtype=bool)
        if hp_cut_arr is not None
        else np.ones(n, dtype=bool)
    ) & valid_base
    str_cut = (
        np.asarray(str_cut_arr, dtype=bool)
        if str_cut_arr is not None
        else np.ones(n, dtype=bool)
    )

    out = np.zeros(100)

    def add_uniform(col, mask):
        """Every masked position contributes weight 1 to the single column `col`."""
        out[col] += np.count_nonzero(mask)

    def add_varying(col_idx, mask):
        """Masked positions contribute weight 1 to their own per-position column."""
        if np.any(mask):
            out[:] += np.bincount(np.asarray(col_idx)[mask], minlength=100)

    # Homopolymer deletion (cols 0-23), one credit per run. A length-1 run
    # is always its own cut, so isolated bases are always credited.
    del_bucket = np.clip(hp_run_arr, 1, 6) - 1
    add_varying(base_cgta * 6 + del_bucket, hp_cut)

    # Homopolymer insertion (cols 24-43), length bins 1..5+, one credit per
    # run. rep0 comes from the 1bpins{base} columns instead.
    ins_bucket = np.clip(hp_run_arr, 1, 5) - 1
    add_varying(24 + base_cgta * 5 + ins_bucket, (hp_run_arr >= 1) & hp_cut)

    # STR deletion (cols 44-67), annotated tracts: one credit per tract.
    unit_bucket = np.clip(str_unit_len_arr, 2, 5) - 2
    count_bucket_del = np.clip(str_repeat_count_arr, 1, 6) - 1
    add_varying(44 + unit_bucket * 6 + count_bucket_del, real_str & str_cut)

    # Flat rep1 credit: every position is a candidate for a non-repeating
    # U-bp deletion for each unit size 2-5+. Not gated by ~real_str: inside
    # a unit-3 tract, a 2bp deletion is still an independent delSTR2_rep1
    # opportunity. Scan-free, matching call.py's coverage-side credit.
    for u_idx in range(4):
        out[44 + u_idx * 6] += n_valid

    # STR insertion (cols 68-91), annotated tracts: one credit per tract;
    # repeat_count 1 is also credited to the "0" bucket.
    real_str_cut = real_str & str_cut
    count_bucket_ins = np.clip(str_repeat_count_arr, 0, 5)
    add_varying(68 + unit_bucket * 6 + count_bucket_ins, real_str_cut)
    is_rep1 = str_repeat_count_arr == 1
    add_varying(68 + unit_bucket * 6 + 0, real_str_cut & is_rep1)

    # Flat rep0 + rep1 insertion credit for every unit size, same
    # reasoning as the deletion case.
    for u_idx in range(4):
        base_col = 68 + u_idx * 6
        out[base_col] += n_valid
        out[base_col + 1] += n_valid

    # Microhomology deletion (cols 92-95): every valid position.
    for k in range(4):
        out[92 + k] += n_valid

    # 1bpins{base} (cols 96-99): rep0 opportunity to insert base b where
    # the following base isn't b (classify_indel_channel's rep0 rule).
    if next_ref_base_arr is not None:
        next_base_valid = (next_ref_base_arr >= 0) & (next_ref_base_arr <= 3)
        for b in range(4):
            add_uniform(96 + b, valid_base & next_base_valid & (next_ref_base_arr != b))

    return out


def load_repeat_context(hp_dataset, str_dataset, start=None, end=None):
    """Merge hp.h5 (rows hp, cut) and str.h5 (rows unit_len, hp, cut;
    unit_len >= 2 only) for one chromosome, or a start:end slice, into a
    single 3-row (unit_len, cut, repeat_count) array. Inside an annotated
    STR (unit_len >= 2) the STR annotation wins; elsewhere the homopolymer
    run is used with unit_len 1.

    hp_dataset, str_dataset: h5py Datasets for one chromosome.
    start, end: optional slice bounds; defaults to the whole chromosome.

    Returns a (3, L) int16 array. Values are <= 127 (clipped by Index.py);
    int16 avoids unsigned wraparound in downstream arithmetic without the
    memory cost of int64 on whole chromosomes.
    """
    sl = slice(start, end)
    hp_run = hp_dataset[0, sl]
    hp_cut = hp_dataset[1, sl]
    str_unit = str_dataset[0, sl]
    str_hp = str_dataset[1, sl]
    str_cut = str_dataset[2, sl]

    is_str = str_unit >= 2
    unit_len_arr = np.where(is_str, str_unit, 1)
    cut_arr = np.where(is_str, str_cut, hp_cut)
    repeat_count_arr = np.where(is_str, str_hp, hp_run)
    return np.vstack((unit_len_arr, cut_arr, repeat_count_arr)).astype(np.int16)


def _count_tandem_copies(seq, unit, from_end=False):
    """Number of back-to-back copies of `unit` at the start of `seq` (or,
    with from_end=True, at its end)."""
    n = len(unit)
    k = 0
    if from_end:
        end = len(seq)
        while end - n >= 0 and seq[end - n : end] == unit:
            k += 1
            end -= n
    else:
        while (k + 1) * n <= len(seq) and seq[k * n : (k + 1) * n] == unit:
            k += 1
    return k


def classify_indel_channel(indel_seq, ref_before, ref_after, indel_len):
    """Resolve a single indel call to its ID83 channel label (internal
    syntax, see Estimate.py's build_indel83_labels), using the same rules
    as SigProfilerMatrixGenerator 1.3.6 (catalogue_generator_INDEL_single,
    limited_indel=True). Uses reference sequence only, not hp.h5/str.h5.

    indel_seq: the deleted or inserted bases (uppercase, length
        abs(indel_len)).
    ref_before: reference bases up to and including the VCF anchor base.
    ref_after: reference bases after the event (after the deleted bases,
        or after the anchor for an insertion).
    indel_len: signed length (+insertion, -deletion).
    Both flanks must reach far enough for the capped bins: 5 copies of
    the event (5*abs(indel_len) bases) on each side is enough.

    Rules (L = abs(indel_len), unit = indel_seq):
      - copies = tandem copies of unit immediately left (ending at the
        anchor) plus immediately right of the event. For a left-aligned
        event the left side is always 0.
      - Deletion, L == 1: {C|T}delhp{1+copies}, capped at 6+ (C+G -> C,
        A+T -> T by the deleted base).
      - Deletion, L >= 2, copies > 0: {L}delstr{1+copies} (6+); the unit
        is the deleted sequence itself, so ACAC from (AC)8 is 4delstr4.
      - Deletion, L >= 2, copies == 0: microhomology = the longer of the
        longest prefix of unit[:-1] matching ref_after and the longest
        suffix of unit[1:] matching ref_before; mh > 0 -> {L}delMH{mh}
        (5+ for L >= 5 and mh >= 5), mh == 0 -> {L}delstr1 (= {L}:Del:R:0,
        plotted as repeat size 1).
      - Insertion, L == 1: {C|T}inshp{copies} (5+).
      - Insertion, L >= 2: {L}insstr{copies} (5+). SPMG's insertion
        microhomology collapses to R:0 under limited_indel, so a
        non-repeat insertion is {L}insstr0 either way.
    L >= 5 is reported in the "5+" length group.
    """
    assert len(indel_seq) == abs(indel_len)
    unit = indel_seq
    L = len(unit)
    copies = _count_tandem_copies(ref_before, unit, from_end=True)
    copies += _count_tandem_copies(ref_after, unit)
    pool = "C" if unit[0] in ("C", "G") else "T"
    len_label = "5+" if L >= 5 else str(L)

    if indel_len > 0:
        count_label = "5+" if copies >= 5 else str(copies)
        if L == 1:
            return f"{pool}inshp{count_label}"
        return f"{len_label}insstr{count_label}"

    count = 1 + copies
    count_label = "6+" if count >= 6 else str(count)
    if L == 1:
        return f"{pool}delhp{count_label}"
    if copies > 0:
        return f"{len_label}delstr{count_label}"

    forward = unit[:-1]
    for_hom = 0
    for i in range(len(forward), 0, -1):
        if ref_after[:i] == forward[:i]:
            for_hom = i
            break
    reverse = unit[1:]
    rev_hom = 0
    for i in range(len(reverse), 0, -1):
        if len(ref_before) >= i and ref_before[-i:] == reverse[-i:]:
            rev_hom = i
            break
    mh = max(for_hom, rev_hom)
    if mh == 0:
        return f"{len_label}delstr1"
    mh_label = "5+" if (L >= 5 and mh >= 5) else str(mh)
    return f"{len_label}delMH{mh_label}"


def indel_coverage_category_index(indel_seq, ref_after, indel_len):
    """Column index (0-15) into INDEL_COVERAGE_CATEGORY_LABELS for one
    indel event. Only the 1bp insertion looks at sequence (does the
    inserted base match the following base?); longer indels are bucketed
    by length alone. 1bp deletions and repeat-extending 1bp insertions go
    to the homopolymer columns 14/15, never the STR-only columns 0/1,
    matching genotypeDSIndel's HP-table scoring of +-1bp events.

    indel_seq/ref_after/indel_len: same arguments as classify_indel_channel.
    """
    if indel_len < 0:
        del_len = -indel_len
        if del_len == 1:
            return 14  # Deletion of Homopolymer
        return min(del_len, 5)  # Deletion Length 2/3/4/5+ -> cols 2-5

    ins_len = indel_len
    if ins_len == 1:
        matches_after = len(ref_after) >= 1 and ref_after[:1] == indel_seq
        if matches_after:
            return 15  # Insertion of Homopolymer
        base2idx = {"A": 0, "T": 1, "C": 2, "G": 3}
        return 6 + base2idx[indel_seq[0].upper()]  # Insertion A/T/C/G
    return 8 + min(ins_len, 5)  # Insertion Length 2/3/4/5+ -> cols 10-13


# ======================================================================
# VCF output helpers
# ======================================================================
def createVcfStrings(chromDict, infoDict, formatDict, filterDict, recs):
    lines = ["##fileformat=VCFv4.2"]
    for filterr in filterDict.keys():
        lines.append(f'##FILTER=<ID={filterr},Description="{filterDict[filterr]}">')
    for info in infoDict.keys():
        lines.append(
            '##INFO=<ID={},Number={},Type={},Description="{}">'.format(
                info, infoDict[info][0], infoDict[info][1], infoDict[info][2]
            )
        )
    for form in formatDict.keys():
        lines.append(
            '##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(
                form, formatDict[form][0], formatDict[form][1], formatDict[form][2]
            )
        )
    for chrom in chromDict.keys():
        lines.append(f"##contig=<ID={chrom},length={chromDict[chrom]}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
    for rec in recs:
        chrom = rec["chrom"]
        pos = rec["pos"]
        alt = rec["alt"]
        ref = rec["ref"]
        infos = rec["infos"]
        filter = rec["filter"]
        formats = rec["formats"]
        samples = rec["samples"]
        lineEntries = [
            chrom,
            str(pos),
            ".",
            ref,
            alt,
            ".",
            filter,
            ";".join([f"{info}={infos[info]}" for info in infoDict.keys()]),
            ":".join(formats),
            ":".join([str(s) for s in samples[0]]),
            ":".join([str(s) for s in samples[1]]),
        ]
        lines.append("\t".join(lineEntries))
    return "\n".join(lines) + "\n"


# ======================================================================
# BAM / region iteration and read-family grouping
# ======================================================================
def bamReadCount(bam, chrom, start, end, ref, regionFile=None):
    if not regionFile:
        count = getAlignmentObject(bam, "rb", ref).count(chrom, start, end)
    else:
        # With -R, count only reads inside the bed intervals.
        regionBed = BED(regionFile)
        count = 0
        if chrom in regionBed.contigs:
            bamObject = getAlignmentObject(bam, "rb", ref)
            for interval in regionBed.fetch(chrom, start, end, parser=pysam.asBed()):
                interval_start = max(interval.start, start)
                interval_end = min(interval.end, end)
                if interval_start < interval_end:
                    count += bamObject.count(chrom, interval_start, interval_end)
    return count


def drop_empty_regions(regionSequence, bamObject):
    """Drop regions spanning zero bases: (contig, start, end) with
    start >= end, or an open-ended (contig, start) at or past the contig
    end. These can come from cut sites at contig boundaries; an empty
    region would write a zero-width coverage.bed.gz interval (breaking
    tabix sort order) or fetch an out-of-range span.
    """
    filtered = []
    for region in regionSequence:
        if len(region) == 3 and region[1] >= region[2]:
            continue
        if len(region) == 2 and region[1] >= bamObject.get_reference_length(region[0]):
            continue
        filtered.append(region)
    return filtered


def splitBamRegions(
    bams, num, contigs, step, ref=None, regionFile=None, min_chunk_length=10000
):
    """Balance read counts without creating chromosome pieces below the minimum.

    Short contigs remain whole. Reject cuts too close to a previous cut or
    either contig end, so the actual worker count can be less than requested.
    """
    if min_chunk_length < 1:
        raise ValueError("min_chunk_length must be a positive integer")
    bamObject = getAlignmentObject(bams[0], "rb", ref)
    contigs_set = set(contigs)
    contigs_sorted = [_ for _ in bamObject.references if _ in contigs]

    total_read_num = 0

    cutSite = [(0, 0)]
    break_flag = 0
    current_read_count = 0
    contig_lens = np.array([bamObject.get_reference_length(c) for c in contigs_sorted])
    window_nums = np.ceil(contig_lens / step).astype(int)
    window_nums_cumulative = np.cumsum(window_nums)
    total_reads_by_windows = np.zeros(window_nums_cumulative[-1])
    # One pool of `num` workers reused for every contig/bam.
    with Pool(num) as pool:
        for bam in bams:
            current_window = 0
            for c, contig in enumerate(contigs_sorted):
                print(
                    f"......Splitting bam regions into roughly equal chunks. Counting reads in {contig} of {bam}......"
                )
                contig_windows = window_nums_cumulative[c] - current_window
                current_arguments = [
                    (
                        bam,
                        contig,
                        k * step,
                        k * step + step,
                        ref,
                        regionFile,
                    )
                    for k in range(contig_windows)
                ]
                counts = np.array(pool.starmap(bamReadCount, current_arguments))
                total_reads_by_windows[
                    current_window : window_nums_cumulative[c]
                ] += counts
                current_window = window_nums_cumulative[c]
    total_reads_by_windows_cumulative = np.cumsum(total_reads_by_windows)
    read_num = total_reads_by_windows_cumulative[-1]
    chunkSize = math.ceil(read_num / num)
    cut_inds = np.searchsorted(
        total_reads_by_windows_cumulative,
        np.arange(1, num) * chunkSize,
    )
    cut_contigs = np.searchsorted(window_nums_cumulative, cut_inds)
    # Interpolate each cut within its window by the fraction of the
    # window's reads needed to reach the target, so one read-dense window
    # can hold several distinct cuts.
    targets = np.arange(1, num) * chunkSize
    # When -p exceeds what the reads can fill, cut_inds can equal
    # total_reads_by_windows.size; cap for indexing (cut_pos is clamped
    # below).
    cut_inds_capped = np.minimum(cut_inds, total_reads_by_windows.size - 1)
    window_start_reads = np.concatenate([[0], total_reads_by_windows_cumulative])[
        cut_inds_capped
    ]
    window_counts = total_reads_by_windows[cut_inds_capped]
    frac = np.where(
        window_counts > 0,
        (targets - window_start_reads) / np.maximum(window_counts, 1),
        0,
    )
    frac = np.clip(frac, 0, 1)
    cut_pos = np.round(
        (cut_inds - np.concatenate([[0], window_nums_cumulative])[cut_contigs]) * step
        + frac * step
    ).astype(int)
    # A contig's last window is usually partial, so cut_pos can overshoot
    # the contig end; clamp it so no region starts past the contig end.
    cut_pos = np.minimum(cut_pos, contig_lens[cut_contigs])
    for nn in range(cut_inds.size):
        candidate = (int(cut_contigs[nn]), int(cut_pos[nn]))
        chrom_idx, pos = candidate
        previous_pos = cutSite[-1][1] if cutSite[-1][0] == chrom_idx else 0
        if (
            pos - previous_pos < min_chunk_length
            or contig_lens[chrom_idx] - pos < min_chunk_length
        ):
            continue
        # Skip duplicate cuts: they would become zero-width regions
        # downstream. Callers tolerate fewer chunks than threads.
        if candidate == cutSite[-1]:
            continue
        cutSite.append(candidate)
    return cutSite, chunkSize, contigs_sorted


def getAlignmentObject(bam, mode, refpath=None):
    if bam.endswith(".bam"):
        bamObject = BAM(bam, mode)
    elif bam.endswith(".cram"):
        bamObject = BAM(bam, mode, reference_filename=refpath)
    else:
        raise NameError(f"{bam} should have .bam or .cram as file extension")
    return bamObject


def parse_barcode_option(value):
    """Parse the --barcode/-bc "TAG,NORMALIZE,SEP" option string into a dict
    of {"tag", "normalize", "sep"}.

    NORMALIZE "1": TAG holds a "bc1<SEP>bc2" pair to be reordered by read
    orientation (like DupCaller's DB tag). "0": TAG already holds an
    orientation-independent family barcode and SEP is unused.
    """
    parts = value.split(",")
    if len(parts) != 3:
        raise ValueError(
            f"--barcode/-bc expects 'TAG,NORMALIZE,SEP' (e.g. 'DB,1,-'), got {value!r}"
        )
    tag, normalize_str, sep = parts
    if not tag:
        raise ValueError(f"--barcode/-bc: TAG must not be empty, got {value!r}")
    if normalize_str not in ("0", "1"):
        raise ValueError(
            f"--barcode/-bc: NORMALIZE must be 0 or 1, got {normalize_str!r} in {value!r}"
        )
    return {"tag": tag, "normalize": normalize_str == "1", "sep": sep}


def get_duplex_barcode(rec, params):
    """Return the (bc1, bc2) duplex barcode pair for an aligned read, from
    params["barcodeTag"] (default "DB"). With barcodeNormalize (default
    True) the tag is split on barcodeSep (default "-"); otherwise it is
    returned as (value, value), so orientation reordering is a no-op and
    strand is tracked from read flags only.
    """
    raw = rec.get_tag(params.get("barcodeTag", "DB"))
    if not params.get("barcodeNormalize", True):
        return raw, raw
    return tuple(raw.split(params.get("barcodeSep", "-")))


def fragment_start(aln):
    """Leftmost reference start of aln's fragment: its own start when it is
    the fragment's leftmost read (template_length >= 0), else its mate's.
    Both mates of a pair, and every read of one molecule, share it."""
    if aln.template_length >= 0:
        return aln.reference_start
    return aln.next_reference_start


def family_fragment_starts(reads):
    """Sorted tuple of the fragment starts of one duplex family's reads --
    stored on every call record ("FS", not written to the VCF) so depth
    extraction can recognize the founding family's reads and their mates."""
    return tuple(sorted({fragment_start(r) for r in reads}))


def bamIterateMultipleRegion(bam, regions, ref, regionFile=None):
    bamObject = getAlignmentObject(bam, "rb", ref)
    if not regionFile:
        for region in regions:
            for rec in bamObject.fetch(*region):
                if len(region) >= 2:
                    if rec.reference_start < region[1]:
                        continue
                yield rec, region
    else:
        # With -R, fetch only reads in bed intervals overlapping the region.
        # A read can overlap several intervals, so dedup within the region;
        # cross-region dedup uses the same reference_start >= region[1]
        # rule as above.
        regionBed = BED(regionFile)
        for region in regions:
            chrom = region[0]
            if chrom not in regionBed.contigs:
                continue
            seen = set()
            for interval in regionBed.fetch(*region, parser=pysam.asBed()):
                interval_start = interval.start
                interval_end = interval.end
                if len(region) >= 2:
                    interval_start = max(interval_start, region[1])
                if len(region) == 3:
                    interval_end = min(interval_end, region[2])
                if interval_start >= interval_end:
                    continue
                for rec in bamObject.fetch(chrom, interval_start, interval_end):
                    if len(region) >= 2 and rec.reference_start < region[1]:
                        continue
                    key = (rec.query_name, rec.is_read1, rec.reference_start)
                    if key in seen:
                        continue
                    seen.add(key)
                    yield rec, region


def trim_mask(n, left, right):
    """Boolean mask over a family's n window positions marking the first
    `left` and last `right` as trimmed. right == 0 trims nothing (a bare
    `mask[-0:]` would select the whole row)."""
    mask = np.zeros(n, dtype=bool)
    if left > 0:
        mask[:left] = True
    if right > 0:
        mask[-right:] = True
    return mask


def determineTrimLength(seq, params, processed_flag):
    if seq.template_length > 0 and not processed_flag:
        overlap = 0  # Never mask overlap of forward read
        left = params["trim5"]
        right_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        right_read = params["trim3"]
        right = max(right_frag, right_read)
    else:
        # Mask overlap of reverse read
        if processed_flag:
            mate_reference_length = _cigar_ref_consumed(seq.get_tag("MC"))
            overlap = max(
                0,
                seq.reference_length
                + mate_reference_length
                - abs(seq.template_length)
                - params["trim3"],
            )
        else:
            overlap = 0
        right_frag = params["trim5"]
        right_read = params["trim3"]
        right = max(right_frag, right_read)
        left_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        left = max(left_frag, overlap, params["trim3"])
    return left, right


def nums2str(nums, num2base="ATCG"):
    bases = [num2base[_] for _ in nums]
    return "".join(bases)


def get_bed_file_for_position(
    pos,
    chrom,
    regions_start_chrom,
    regions_start_pos,
    regions_end_chrom,
    regions_end_pos,
    locus_bed,
    locus_bed_prev,
    locus_bed_next,
):
    """
    Determine which bed file to write to based on position relative to region boundaries
    Only compare positions within the same chromosome
    """
    if chrom == regions_start_chrom and pos < regions_start_pos:
        return locus_bed_prev
    elif chrom == regions_end_chrom and pos >= regions_end_pos:
        return locus_bed_next
    else:
        return locus_bed


def _filter_reason_label(
    ncov_flag, nm_flag, trim_flag, extra_flag=None, extra_label=None
):
    """Filter label for a --rescue candidate naming the mask that blocks
    it: n_cov_mask, nm_mask, trim, or (indel only, via extra_flag/
    extra_label) indel_mask. Args are bools or bool arrays over the
    candidate's position(s). Checked in fixed priority order so exactly one
    label is returned.
    """
    if extra_flag is not None and np.any(extra_flag):
        return extra_label
    if np.any(ncov_flag):
        return "cov_mask"
    if np.any(nm_flag):
        return "nm_mask"
    if np.any(trim_flag):
        return "trim_mask"
    return "other_mask"


def _compute_read_label(rec, params):
    """Duplex-family label for rec: barcode pair plus signed template_length."""
    bc1, bc2 = get_duplex_barcode(rec, params)
    if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
        return bc1 + "+" + bc2 + "+" + str(rec.template_length)
    else:
        return bc2 + "+" + bc1 + "+" + str(rec.template_length)


def _place_read_in_dict(rec, label, currentReadDict):
    """Add rec to currentReadDict under label, merging into an existing
    family entry if one is already open."""
    if label in currentReadDict:
        currentReadDict[label]["seqs"].append(rec)
        currentReadDict[label]["names"][rec.query_name] = (
            len(currentReadDict[label]["seqs"]) - 1
        )
    else:
        currentReadDict[label] = {
            "seqs": [rec],
            "F1R2": 0,
            "F2R1": 0,
            "names": {rec.query_name: 0},
        }
    if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
        currentReadDict[label]["F1R2"] += 1
    else:
        currentReadDict[label]["F2R1"] += 1


def _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool):
    """For each upstream-anchor family (template_length > 0) whose members
    disagree on next_reference_start, redirect (via rugged_reads_index)
    every read whose mate is not at the LARGEST mate position, and pre-open
    that position's rugged_reads_pool bucket.

    The largest position, not the mode: the BAM is streamed in coordinate
    order and a bucket is drained when the scan reaches it, so it must sit
    at or after every member's arrival. This affects only when the family
    is processed, not its coordinates.

    Keys are (chrom, position) so a bucket that is never drained can't
    collide with the same position on a later chromosome.
    """
    for group in currentReadDict.values():
        group_seqs = group["seqs"]
        if group_seqs[0].template_length <= 0:
            continue
        chrom = group_seqs[0].reference_name
        mate_starts = [r.next_reference_start for r in group_seqs]
        largest_mate_start = max(mate_starts)
        if any(ms != largest_mate_start for ms in mate_starts):
            pool_key = (chrom, largest_mate_start)
            rugged_reads_pool.setdefault(pool_key, [])
            for r, ms in zip(group_seqs, mate_starts):
                if ms != largest_mate_start:
                    rugged_reads_index[r.query_name] = pool_key


def _drain_rugged_pool(chrom, position, rugged_reads_pool, params, currentReadDict):
    """Merge any reads pooled for (chrom, position) into currentReadDict."""
    pooled = rugged_reads_pool.pop((chrom, position), None)
    if pooled:
        for pooled_rec in pooled:
            _place_read_in_dict(
                pooled_rec, _compute_read_label(pooled_rec, params), currentReadDict
            )


_REF_CONSUMING_CIGAR_OPS = set("MDN=X")


def _cigar_ref_consumed(cigar_str):
    """Reference length consumed by a CIGAR string (M/D/N/=/X ops)."""
    consumed = 0
    for length_str, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar_str):
        if op in _REF_CONSUMING_CIGAR_OPS:
            consumed += int(length_str)
    return consumed


def _eligible_overflow_read(rec):
    return (
        rec.is_proper_pair
        and not rec.is_supplementary
        and not rec.is_secondary
        and not rec.is_qcfail
    )


def check_mate_cigar_tags(bam_path, reference=None, sample_size=100_000):
    """Exit with an error if a proper-pair read (see _eligible_overflow_read)
    lacks the MC tag that chunk-boundary overflow relies on. Checks only
    the first `sample_size` records, so MC tags missing only deep in a
    merged BAM can go undetected.
    """
    print(
        f"Checking MC (mate CIGAR) tags on proper pairs in {bam_path} "
        f"(sampling up to {sample_size} reads)...",
        flush=True,
    )
    with getAlignmentObject(bam_path, "rb", reference) as bam:
        for n_checked, rec in enumerate(bam.fetch(until_eof=True)):
            if n_checked >= sample_size:
                break
            if _eligible_overflow_read(rec) and not rec.has_tag("MC"):
                print(
                    f"ERROR: Proper-pair read {rec.query_name!r} in {bam_path} is missing "
                    "the required MC (mate CIGAR) tag. DupCaller cannot safely "
                    "process chunk boundaries without MC tags. Add MC tags "
                    "to proper-pair reads in the tumor alignment file before rerunning.",
                    file=sys.stderr,
                )
                sys.exit(1)


def previous_region_coverage_start(bam, previous_regions, ref=None, region_file=None):
    """Return (chrom, exclusive end) of the preceding worker's coverage,
    scanning its last region with the same MC-based overflow iterator
    used for calling.
    """
    if not previous_regions or len(previous_regions[-1]) != 3:
        return None
    region = previous_regions[-1]
    end = region[2]
    for rec, _ in bamIterateMultipleRegionWithOverflow(bam, [region], ref, region_file):
        if _eligible_overflow_read(rec) and rec.reference_end is not None:
            end = max(end, rec.reference_end)
    return region[0], end


def _mate_reference_end(rec):
    """Mate's reference_end from rec's MC tag and next_reference_start (no
    mate() lookup). None if rec has no MC tag."""
    if not rec.has_tag("MC"):
        return None
    return rec.next_reference_start + _cigar_ref_consumed(rec.get_tag("MC"))


def bamIterateMultipleRegionWithOverflow(bam, regions, ref, region_file=None):
    """Like bamIterateMultipleRegion, but if regions[-1] has an explicit
    end, also yields an overflow continuation past it, up to the furthest
    mate end (from MC tags) of reads crossing the boundary, so the worker
    that saw a family's upstream anchor can reconstruct the whole family.

    Exactly-once: the overflow pass yields only reads whose mate starts
    before the boundary, and callBam's prev_boundary check makes the next
    worker skip those same reads.
    """
    boundary_chrom = None
    boundary_pos = None
    if len(regions[-1]) == 3:
        boundary_chrom = regions[-1][0]
        boundary_pos = regions[-1][2]
    last_end = boundary_pos

    for rec, region in bamIterateMultipleRegion(bam, regions, ref, region_file):
        if (
            boundary_chrom is not None
            and rec.reference_name == boundary_chrom
            and _eligible_overflow_read(rec)
        ):
            mate_end = _mate_reference_end(rec)
            if mate_end is not None and mate_end > last_end:
                last_end = mate_end
        yield rec, region

    if boundary_chrom is not None and last_end > boundary_pos:
        # Fetch ends are exclusive: a mate starting exactly at the boundary
        # belongs to this continuation, not to the preceding normal fetch.
        overflow_region = (boundary_chrom, boundary_pos, int(last_end) + 1)
        for rec, region in bamIterateMultipleRegion(
            bam, [overflow_region], ref, region_file
        ):
            if rec.next_reference_start >= boundary_pos:
                # Belongs to the next worker's normal scan.
                continue
            yield rec, region


# ======================================================================
# Error-matrix normalization and loading
# ======================================================================
def _normalize_indel_hp_mat(mat, pseudocount, fallback=None):
    """Normalize a raw (10, 12) hp.txt count matrix (rows hp run length
    1-10+, columns ref_allele*3+(idLen+1) for idLen in {-1,0,1}) into
    probabilities: each (hp_len, base) 3-column group is Dirichlet-smoothed
    by `pseudocount` and divided by its sum.

    `fallback` (bundled fallback_latest hp.txt counts): in a group with
    fewer than FALLBACK_MIN_SITES observations, a zero-count del/ins entry
    takes the fallback rate. Rows are then made non-decreasing in hp_len
    (row n >= row n-1 elementwise, using row n-1's final rates), per base
    group. Without `fallback`, an all-zero row copies the previous row.
    """
    mat_new = np.zeros_like(mat, dtype=float)
    for g in range(4):
        block = mat[:, g * 3 : g * 3 + 3]
        row_sum = block.sum(axis=1, keepdims=True)
        block_new = (block + pseudocount) / (row_sum + 3 * pseudocount)
        if fallback is not None:
            fb_block = fallback[:, g * 3 : g * 3 + 3]
            fb_new = (fb_block + pseudocount) / (
                fb_block.sum(axis=1, keepdims=True) + 3 * pseudocount
            )
            # Gate on this (hp_len, base) group's own site total.
            use_fb = (block == 0) & (row_sum < FALLBACK_MIN_SITES)
            # del (col 0) / ins (col 2) only; col 1 is the opportunity.
            use_fb[:, 1] = False
        for nn in range(block_new.shape[0]):
            if fallback is not None and use_fb[nn].any():
                block_new[nn, use_fb[nn]] = fb_new[nn, use_fb[nn]]
            if nn == 0:
                continue
            if fallback is None and row_sum[nn, 0] == 0:
                block_new[nn, :] = block_new[nn - 1, :]
                continue
            current_row = block_new[nn, :]
            smaller_entries = current_row <= block_new[nn - 1, :]
            current_row[smaller_entries] = block_new[nn - 1, :][smaller_entries]
            block_new[nn, :] = current_row
        mat_new[:, g * 3 : g * 3 + 3] = block_new
    return mat_new


def _normalize_indel_str_mat(mat, pseudocount, fallback=None):
    """Normalize a raw (5, 15) str.txt count matrix (rows STR-length bin
    0 = not a real repeat through 4 = "40+", columns a (deletion,
    reference, insertion) triple per repeat-unit group u1-u5+) into
    probabilities: each (row, unit group) triple is Dirichlet-smoothed by
    `pseudocount` and divided by its sum.

    `fallback` (bundled fallback_latest str.txt counts): in a group with
    fewer than FALLBACK_MIN_SITES observations, a zero-count deletion/
    insertion entry takes the fallback rate; in rows 2-4 it is raised to
    the previous row's final rate if higher, unless that previous row has
    no fallback opportunity (its rate would be the flat pseudocount 1/3).
    Without `fallback`, an all-zero group in rows 2-4 copies the previous
    row's rates if that row has any.
    """
    if mat.shape != (5, STR_NCOL):
        raise ValueError(
            f"str.txt error profile has shape {mat.shape}, expected (5, {STR_NCOL}) "
            "(per-repeat-unit columns, DupCaller >= 1.2.12); relearn the error "
            "profile with this DupCaller version"
        )
    mat_new = np.zeros_like(mat, dtype=float)
    for g in range(STR_UNIT_GROUPS):
        cols = slice(g * 3, g * 3 + 3)
        block = mat[:, cols]
        row_sum = block.sum(axis=1, keepdims=True)
        block_new = (block + pseudocount) / (row_sum + 3 * pseudocount)
        if fallback is not None:
            fb_block = fallback[:, cols]
            fb_new = (fb_block + pseudocount) / (
                fb_block.sum(axis=1, keepdims=True) + 3 * pseudocount
            )
            use_fb = (block == 0) & (row_sum < FALLBACK_MIN_SITES)
            use_fb[:, 1] = False
            for r in range(5):
                fb_cols = use_fb[r]
                if not fb_cols.any():
                    continue
                if r >= 2 and fb_block[r - 1, 1] > 0:
                    block_new[r, fb_cols] = np.maximum(
                        fb_new[r, fb_cols], block_new[r - 1, fb_cols]
                    )
                else:
                    block_new[r, fb_cols] = fb_new[r, fb_cols]
        else:
            has_rates = row_sum[:, 0] > 0
            for r in range(2, 5):
                if not has_rates[r] and has_rates[r - 1]:
                    block_new[r, :] = block_new[r - 1, :]
                    has_rates[r] = True
        mat_new[:, cols] = block_new
    return mat_new


def str_profile_is_current(path):
    """True if an .amp/.dmg.str.txt profile has the per-repeat-unit layout
    (STR_NCOL data columns); older profiles have 11."""
    with open(path) as fh:
        return len(fh.readline().rstrip("\n").split("\t")) - 1 == STR_NCOL


# In a context with fewer observations than this, a zero-count error
# category is treated as under-sampled and takes the bundled
# fallback_latest.* rate instead.
FALLBACK_MIN_SITES = 1000


def fallback_error_file(suffix):
    """Path of the bundled fallback_latest{suffix} error profile. Looked up
    in DupCaller_sub/ERROR (where setup.py installs src/ERROR) first, then
    src/ERROR itself (source checkout / editable install)."""
    pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates = [
        os.path.join(pkg_dir, "ERROR", "fallback_latest" + suffix),
        os.path.join(os.path.dirname(pkg_dir), "ERROR", "fallback_latest" + suffix),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Fallback error profile fallback_latest{suffix} not found in: "
        + ", ".join(candidates)
    )


def _read_fallback_counts(suffix):
    return (
        pd.read_csv(fallback_error_file(suffix), sep="\t", index_col=0)
        .to_numpy()
        .astype(float)
    )


def apply_sbs_low_coverage_fallback(
    rate_mat, counts, fallback_rate_mat, min_sites=FALLBACK_MIN_SITES
):
    """In a (64, 4) SBS rate matrix (SRD or SSM/damage), replace each alt
    (error) entry whose own count is 0 in a row with fewer than min_sites
    total observations by the bundled fallback profile's rate. Observed
    (non-zero) alt entries and well-sampled rows are left alone; a touched
    row's reference column is reset to 1 - sum(alt rates)."""
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    _, num2trinuc = build_trinuc64_order()
    counts = np.asarray(counts)
    out = np.array(rate_mat, dtype=float)
    site_totals = counts.sum(axis=1)
    for row, trinuc in enumerate(num2trinuc):
        if site_totals[row] >= min_sites:
            continue
        ref_col = base2num[trinuc[1]]
        zero = counts[row] == 0
        zero[ref_col] = False
        if not zero.any():
            continue
        out[row, zero] = fallback_rate_mat[row, zero]
        out[row, ref_col] = 1 - (out[row].sum() - out[row, ref_col])
    return out


def ref_conditional_error_rates(rates, ref_col, fallback_rates):
    """Turn a per-context error-rate matrix into the binary-model rates the
    genotyper consumes: P(error | error or reference).

    Per row, Pref = 1 - sum(other entries). A row with Pref < 0 is replaced
    by that row of `fallback_rates` and recomputed. Each non-reference
    entry P becomes P / (Pref + P); the reference column holds Pref.

    rates/fallback_rates: (rows, cols); ref_col: scalar or per-row column
    index of the reference entry.
    """
    out = np.array(rates, dtype=float)
    n = out.shape[0]
    idx = np.arange(n)
    ref_col = np.broadcast_to(np.asarray(ref_col), (n,))
    alt_mask = np.ones_like(out, dtype=bool)
    alt_mask[idx, ref_col] = False
    ref = 1 - np.where(alt_mask, out, 0).sum(axis=1)
    bad = ref < 0
    if bad.any():
        out[bad] = np.asarray(fallback_rates, dtype=float)[bad]
        ref[bad] = 1 - np.where(alt_mask[bad], out[bad], 0).sum(axis=1)
        if (ref < 0).any():
            raise ValueError(
                "Fallback error profile has a row whose error rates sum past 1"
            )
    denom = out + ref[:, None]
    cond = np.divide(out, denom, out=np.zeros_like(out), where=denom > 0)
    out[alt_mask] = cond[alt_mask]
    out[idx, ref_col] = ref
    return out


def _sbs_ref_cols():
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    _, num2trinuc = build_trinuc64_order()
    return np.array([base2num[trinuc[1]] for trinuc in num2trinuc])


def _hp_ref_conditional(rates, fallback_counts, pseudocount):
    """ref_conditional_error_rates per (hp_len, base) 3-column group of a
    normalized (10, 12) hp matrix: ref (idLen=0, col 1) = 1 - del - ins."""
    out = np.array(rates, dtype=float)
    for g in range(4):
        fb_block = fallback_counts[:, g * 3 : g * 3 + 3]
        fb_rates = (fb_block + pseudocount) / (
            fb_block.sum(axis=1, keepdims=True) + 3 * pseudocount
        )
        out[:, g * 3 : g * 3 + 3] = ref_conditional_error_rates(
            out[:, g * 3 : g * 3 + 3], 1, fb_rates
        )
    return out


def _str_ref_conditional(rates, fallback_counts, pseudocount):
    """ref_conditional_error_rates per (STR row, unit group) 3-column group
    of a normalized (5, 15) str matrix: ref = 1 - deletion - insertion."""
    out = np.array(rates, dtype=float)
    for g in range(STR_UNIT_GROUPS):
        cols = slice(g * 3, g * 3 + 3)
        fb_block = fallback_counts[:, cols]
        fb_rates = (fb_block + pseudocount) / (
            fb_block.sum(axis=1, keepdims=True) + 3 * pseudocount
        )
        out[:, cols] = ref_conditional_error_rates(out[:, cols], 1, fb_rates)
    return out


def load_error_matrices(params):
    """Load and regularize the SBS (amperr/dmgerr) and indel
    (amperri/dmgerri) error-rate matrices named in params, populating
    params in place with ampmat/ampmat_rev/dmgmat_top/dmgmat_rev_top/
    dmgmat_bot/dmgmat_rev_bot/ampmat_hp/dmgmat_hp/ampmat_str/dmgmat_str/
    trinuc_convert/trinuc2num_dict/num2trinuc_list (the inputs to callBam's
    power tables). Also used by Caller.py to re-simulate detection power
    without rescanning the BAM.

    With params["isLearn"] (default False) the amp matrices are zeros and
    dmg is uniform; neither is used for calling there.
    """
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    isLearn = params.get("isLearn", False)

    trinuc2num, num2trinuc = build_trinuc64_order()

    trinuc_convert_np = np.zeros([64, 4], dtype=np.uint8)
    for trinuc in trinuc2num.keys():
        row = np.zeros(4)
        row_num = trinuc2num[trinuc]
        row[0] = trinuc2num[trinuc[0] + "A" + trinuc[2]]
        row[1] = trinuc2num[trinuc[0] + "T" + trinuc[2]]
        row[2] = trinuc2num[trinuc[0] + "C" + trinuc[2]]
        row[3] = trinuc2num[trinuc[0] + "G" + trinuc[2]]
        trinuc_convert_np[row_num, :] = row
    params["trinuc_convert"] = trinuc_convert_np
    params["trinuc2num_dict"] = trinuc2num
    params["num2trinuc_list"] = num2trinuc
    # Load amp error matrix
    if isLearn:
        ampmat = np.zeros([64, 4])
    else:
        # amperr_file is an SRD rate matrix (learn.py's
        # estimate_sbs_srd_rates): rows are already smoothed probability
        # distributions, so no row-sum normalization here.
        ampmat = (
            pd.read_csv(params["amperr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )
    # Convert to the binary model's P(error | error or ref).
    if not isLearn:
        ampmat = ref_conditional_error_rates(
            ampmat, _sbs_ref_cols(), _read_fallback_counts(".amp.tn.srd.txt")
        )
    ampmat_min_error = ampmat.min(axis=1, keepdims=True)
    ampmat = np.concatenate([ampmat, ampmat_min_error], axis=1)
    params["ampmat"] = ampmat

    ampmat_rev = np.zeros([64, 4])
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            ampmat_rev[trinuc2num[trinuc], nn] = ampmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    ampmat_rev_min_error = ampmat_rev.min(axis=1, keepdims=True)
    ampmat_rev = np.concatenate([ampmat_rev, ampmat_rev_min_error], axis=1)
    params["ampmat_rev"] = ampmat_rev

    if isLearn:
        params["ampmat_hp"] = np.zeros([10, 12])
        params["ampmat_str"] = np.zeros([5, STR_NCOL])
    else:
        # The indel profile is stored as separate .hp.txt/.str.txt files.
        amperri_hp_file = params["amperri_file"].replace(".amp.id.txt", ".amp.hp.txt")
        amperri_str_file = params["amperri_file"].replace(".amp.id.txt", ".amp.str.txt")
        ampmat_hp = pd.read_csv(amperri_hp_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        ampmat_str = pd.read_csv(amperri_str_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        pseudocount = params.get("pseudocount", 0.5)
        amp_hp_fb = _read_fallback_counts(".amp.hp.txt")
        ampmat_hp = _normalize_indel_hp_mat(ampmat_hp, pseudocount, amp_hp_fb)
        ampmat_hp = _hp_ref_conditional(ampmat_hp, amp_hp_fb, pseudocount)
        params["ampmat_hp"] = ampmat_hp

        amp_str_fb = _read_fallback_counts(".amp.str.txt")
        ampmat_str = _normalize_indel_str_mat(ampmat_str, pseudocount, amp_str_fb)
        ampmat_str = _str_ref_conditional(ampmat_str, amp_str_fb, pseudocount)
        params["ampmat_str"] = ampmat_str

    # Load damage matrix
    if isLearn:
        dmgmat = np.zeros([64, 4])
    else:
        dmgmat = (
            pd.read_csv(params["dmgerr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )

    # Dirichlet-smoothed on the counts: a zero-observation context lands
    # at the uniform 1/4.
    pseudocount = params.get("pseudocount", 0.5)
    dmgmat_counts = dmgmat
    dmgmat_row_sum = dmgmat.sum(axis=1, keepdims=True)
    dmgmat = (dmgmat + pseudocount) / (dmgmat_row_sum + 4 * pseudocount)
    if not isLearn:
        # Unobserved error categories in under-sampled rows take the
        # fallback profile's rates, smoothed the same way as the sample's
        # own counts.
        fallback_counts = _read_fallback_counts(".dmg.tn.txt")
        fallback_rates = (fallback_counts + pseudocount) / (
            fallback_counts.sum(axis=1, keepdims=True) + 4 * pseudocount
        )
        dmgmat = apply_sbs_low_coverage_fallback(dmgmat, dmgmat_counts, fallback_rates)
        dmgmat = ref_conditional_error_rates(dmgmat, _sbs_ref_cols(), fallback_rates)
    dmgmat_ref_error = dmgmat.min(axis=1, keepdims=True)
    dmgmat = np.concatenate([dmgmat, dmgmat_ref_error], axis=1)

    params["dmgmat_top"] = dmgmat
    params["trinuc2num_dict"] = trinuc2num

    dmgmat_rev = np.zeros([64, 4])
    dmgmat_rev_ref_error = np.zeros(64)
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            dmgmat_rev[trinuc2num[trinuc], nn] = dmgmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    dmgmat_rev_ref_error = dmgmat_rev.min(axis=1, keepdims=True)
    dmgmat_rev = np.concatenate([dmgmat_rev, dmgmat_rev_ref_error], axis=1)
    params["dmgmat_rev_top"] = dmgmat_rev

    dmgmat_b = np.vstack((dmgmat[32:64, [1, 0, 3, 2]], dmgmat[:32, [1, 0, 3, 2]]))
    dmgmat_b_ref_error = dmgmat_b.min(axis=1, keepdims=True)
    dmgmat_rev_b = np.vstack(
        (dmgmat_rev[32:64, [1, 0, 3, 2]], dmgmat_rev[:32, [1, 0, 3, 2]])
    )
    dmgmat_rev_b_ref_error = dmgmat_rev_b.min(axis=1, keepdims=True)
    dmgmat_b = np.concatenate([dmgmat_b, dmgmat_b_ref_error], axis=1)
    dmgmat_rev_b = np.concatenate([dmgmat_rev_b, dmgmat_rev_b_ref_error], axis=1)
    params["dmgmat_bot"] = dmgmat_b
    params["dmgmat_rev_bot"] = dmgmat_rev_b

    if isLearn:
        params["dmgmat_hp"] = np.zeros([10, 12])
        params["dmgmat_str"] = np.zeros([5, STR_NCOL])
    else:
        dmgerri_hp_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.hp.txt")
        dmgerri_str_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.str.txt")
        dmgmat_hp = pd.read_csv(dmgerri_hp_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        dmgmat_str = pd.read_csv(dmgerri_str_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        dmg_hp_fb = _read_fallback_counts(".dmg.hp.txt")
        dmgmat_hp = _normalize_indel_hp_mat(dmgmat_hp, pseudocount, dmg_hp_fb)
        dmgmat_hp = _hp_ref_conditional(dmgmat_hp, dmg_hp_fb, pseudocount)
        params["dmgmat_hp"] = dmgmat_hp

        dmg_str_fb = _read_fallback_counts(".dmg.str.txt")
        dmgmat_str = _normalize_indel_str_mat(dmgmat_str, pseudocount, dmg_str_fb)
        dmgmat_str = _str_ref_conditional(dmgmat_str, dmg_str_fb, pseudocount)
        params["dmgmat_str"] = dmgmat_str


# ======================================================================
# FDR / detection-power calibration
# ======================================================================
def threshold_rng(base_seed, threshold):
    """RNG for one simulate_power_grid draw at LR `threshold`, seeded from
    (base_seed, threshold) so results don't depend on worker count or
    scheduling order.

    The threshold enters as round(|threshold| * 1e4) mod 1e5 (4 decimals,
    one integer digit), mixed with base_seed via SeedSequence. An infinite
    threshold (_refine_channel's reject-all sentinel) uses a fixed value;
    the seed is irrelevant there since the grid comparison is constant.
    """
    if np.isfinite(threshold):
        threshold_digits = int(round(abs(threshold) * 10000)) % 100000
    else:
        threshold_digits = 99999
    return np.random.default_rng((int(base_seed), threshold_digits))


def simulate_power_grid(
    Pamp_c,
    Pamp_rev_c,
    Pdmg_c,
    Pdmg_rev_c,
    Pdmg_bot_c,
    Pdmg_rev_bot_c,
    threshold,
    all_quals,
    rng=None,
    N_SIM=100,
):
    """Monte Carlo (10, 10) detection-power grid: for every (top-strand,
    bottom-strand) read-family count composition, the fraction of N_SIM
    simulated quality draws where calculateDSPosterior's log-likelihood
    ratio (LL_B1 - LL_B2) clears `threshold`, given this single context's
    fixed amplification/damage probabilities.

    Shared by callBam's power tables (L, L_indel_1bp, L_indel_len) and
    Caller.py's re-simulation at a refined threshold, so both use
    identical math.
    """
    if rng is None:
        rng = np.random.default_rng()
    ln10 = np.log(10)
    grid = np.zeros([10, 10])
    P_arr = np.full(N_SIM, Pamp_c)
    P_rev_arr = np.full(N_SIM, Pamp_rev_c)
    Pt_arr = np.full(N_SIM, Pdmg_c)
    Prev_t_arr = np.full(N_SIM, Pdmg_rev_c)
    Pb_arr = np.full(N_SIM, Pdmg_bot_c)
    Prev_b_arr = np.full(N_SIM, Pdmg_rev_bot_c)
    for i in range(10):
        for j in range(10):
            if i == 0 and j == 0:
                continue
            if i > 0:
                q_top = rng.choice(all_quals, size=(i, N_SIM), replace=True)
                F1R2_Pseq = -q_top / 10 * ln10
                F1R2_bin = np.ones([i, N_SIM], dtype=bool)
            else:
                F1R2_Pseq = np.zeros([0, N_SIM])
                F1R2_bin = np.zeros([0, N_SIM], dtype=bool)
            if j > 0:
                q_bot = rng.choice(all_quals, size=(j, N_SIM), replace=True)
                F2R1_Pseq = -q_bot / 10 * ln10
                F2R1_bin = np.ones([j, N_SIM], dtype=bool)
            else:
                F2R1_Pseq = np.zeros([0, N_SIM])
                F2R1_bin = np.zeros([0, N_SIM], dtype=bool)
            F1R2_b1, F1R2_b2 = calculateSSPosterior(
                P_arr, P_rev_arr, F1R2_bin, F1R2_Pseq
            )
            F2R1_b1, F2R1_b2 = calculateSSPosterior(
                P_arr, P_rev_arr, F2R1_bin, F2R1_Pseq
            )
            LL_B1, LL_B2 = calculateDSPosterior(
                Pt_arr,
                Prev_t_arr,
                Pb_arr,
                Prev_b_arr,
                F1R2_b1,
                F2R1_b1,
                F1R2_b2,
                F2R1_b2,
            )
            grid[i, j] = ((LL_B1 - LL_B2) >= threshold).mean()
    return grid


# depth_by_hpstr's context axis:
#   0-39   homopolymer: base * 10 + run length (1-10+) - 1
#   40-71  STR bin 1-4 x repeat unit 2-5+ x deletion/insertion: slip sites
#          (depth_hpstr_str_bucket)
#   72     STR0 1bp (mismatched insertion): every position
#   73-80  STR0 indel length 2-5+ x deletion/insertion: sites of a non-slip
#          indel of that length (depth_hpstr_str0_bucket)
# Deletions and insertions are split because their sites and masks differ:
# a deletion needs every deleted base unmasked and can only sit where it
# can't left-align further; an insertion that isn't a slip can sit
# anywhere.
DEPTH_HPSTR_STR_BASE = 40
DEPTH_HPSTR_STR0_ANY = 72
DEPTH_HPSTR_STR0_LEN_BASE = 73
DEPTH_HPSTR_NCAT = 81
STR_EEFF_LENGTHS = (2, 3, 4, 5)


def _signed_len_offset(indel_len):
    """Offset of a signed indel length (STR0) or repeat unit (STR bins 1-4)
    2-5+ within its 8 buckets: pairs, deletion first."""
    return (min(abs(indel_len), 5) - 2) * 2 + (1 if indel_len > 0 else 0)


def depth_hpstr_str_bucket(str_bin, indel_len):
    """depth_by_hpstr bucket of STR bin 1-4 and signed repeat unit 2-5+
    (negative: deletion)."""
    return DEPTH_HPSTR_STR_BASE + (str_bin - 1) * 8 + _signed_len_offset(indel_len)


def depth_hpstr_str0_bucket(indel_len):
    """depth_by_hpstr bucket of an STR0 channel of this signed indel length."""
    if abs(indel_len) == 1:
        return DEPTH_HPSTR_STR0_ANY
    return DEPTH_HPSTR_STR0_LEN_BASE + _signed_len_offset(indel_len)


def str_tract_valid(str_cut_bool, total_len):
    """STR tract starts of one family window whose tract can be scored: its
    anchor (the base before the start) and the base after the tract
    (start + unit * repeat count) both fall inside the window."""
    p = np.arange(str_cut_bool.shape[0])
    return str_cut_bool & (p >= 1) & (p + total_len < str_cut_bool.shape[0])


def indel_eeff_site_masks(
    antimask_indel,
    ref_int,
    hp_cut_bool,
    hp_repeat_valid,
    str_unit,
    str_total_len,
    str_cut_bool,
    str_repeat_valid,
):
    """Positions of one family window that count toward each indel
    channel's Eeff (depth_by_hpstr): exactly where a candidate of that
    channel can be called. A site is the context base (anchor + 1). Every
    site needs the anchor and the context base unmasked (antimask_indel, the
    indel masks, no snp_mask), an ACGT context base and an anchor that isn't
    the window's first base (indel_has_context); a deletion of k also needs
    all k deleted bases unmasked and inside the window (the same span as
    indel_mask_span), and its last deleted base to differ from the anchor
    (otherwise it left-aligns further).

    HP: one site per run, its first base, skipping runs cut by the window
    edge (hp_repeat_valid). STR(+-u): tract starts of repeat unit u (any
    number of units slips, str_slip_tract), whose tract fits the window
    (str_repeat_valid); a deletion site needs the one-unit deletion's span.
    STR0(+k): every site (an insertion that isn't copies of the tract's
    unit can go anywhere). STR0(-k): every deletion site that isn't a slip
    (a tract start whose unit divides k). The capped 5+ unit/length pools
    5, 6, 7, ...; it is approximated with k = 5's mask span, every tract
    start of unit 5+ for STR and no left-alignment test.

    str_unit: the STR unit length at each position (0/1 outside an STR);
    str_total_len: the tract length (unit * repeat count).
    Returns (hp, str_by_unit, str0_any, str0_by_len), str_by_unit keyed by
    signed repeat unit and str0_by_len by signed length, both
    +-STR_EEFF_LENGTHS."""
    n = antimask_indel.shape[0]
    site_ok = antimask_indel & (ref_int < 4)
    site_ok[1:] &= antimask_indel[:-1]
    site_ok[:2] = False
    p = np.arange(n)
    # masked_before[i]: masked positions in [0, i).
    masked_before = np.concatenate(([0], np.cumsum(~antimask_indel)))
    is_str = str_unit >= 2
    unit = np.maximum(str_unit, 1)
    # str_length_bin: a tract needs two units to be a repeat.
    str_start = str_cut_bool & is_str & (str_total_len >= 2 * unit)
    str_by_unit = {}
    str0_by_len = {}
    for k in STR_EEFF_LENGTHS:
        end = np.minimum(p + k, n)
        del_ok = site_ok & (p + k <= n) & (masked_before[end] - masked_before[p] == 0)
        if k == STR_EEFF_LENGTHS[-1]:
            unit_start = str_start & (unit >= k)
            str_by_unit[-k] = unit_start & str_repeat_valid & del_ok
            str_by_unit[k] = unit_start & str_repeat_valid & site_ok
            str0_by_len[-k] = del_ok
            str0_by_len[k] = site_ok
            continue
        # Last deleted base of a k-deletion at p (ref_int[p + k - 1]),
        # compared with the anchor.
        after = np.full(n, -1)
        m = max(n - k, 0)
        after[1 : 1 + m] = ref_int[k : k + m]
        before = np.full(n, -2)
        before[1:] = ref_int[:-1]
        left_aligned = after != before
        unit_start = str_start & (unit == k)
        str_by_unit[-k] = unit_start & left_aligned & str_repeat_valid & del_ok
        str_by_unit[k] = unit_start & str_repeat_valid & site_ok
        unit_divides = str_start & (k % unit == 0)
        slip_del = unit_divides & (str_total_len >= k) & left_aligned
        str0_by_len[-k] = del_ok & left_aligned & ~slip_del
        str0_by_len[k] = site_ok
    return site_ok & hp_repeat_valid & hp_cut_bool, str_by_unit, site_ok, str0_by_len


def str_channel_name(str_bin, signed_unit):
    """Name of an STR mu-solve channel: STR{bin}_u{U}_{del|ins}. U is the
    repeat unit 2-5 (5 = 5+) in STR bins 1-4 and the indel length 1-5
    (5 = 5+) in STR0 (not a slip; STR0_u1_ins is the mismatched 1bp
    insertion)."""
    sign = "del" if signed_unit < 0 else "ins"
    return f"STR{str_bin}_u{min(abs(int(signed_unit)), 5)}_{sign}"


def parse_str_channel_name(name):
    """(str_bin, signed unit) of a str_channel_name name."""
    bin_str, unit_str, sign = name[3:].split("_")
    unit = int(unit_str[1:])
    return int(bin_str), -unit if sign == "del" else unit


def indel_mu_channel(id_len, hp_len, str_bin, hm, base_char, unit_len=0):
    """Name of the mu-solve channel an indel candidate belongs to, or None
    (an HP-routed 1bp event outside any homopolymer). Same branch selection
    as genotypeDSIndel/indelErrorProbs: |id_len| >= 2 is STR-context, by
    the tract's repeat unit unit_len in STR bins 1-4 and by the length in
    STR0; a 1bp mismatched insertion (hm == 0) is STR0_u1_ins; any other
    1bp event is HP, pooled by base_char (the deleted or inserted base):
    "C" for C/G, "T" for A/T."""
    if abs(id_len) >= 2:
        unit = unit_len if str_bin >= 1 else abs(id_len)
        return str_channel_name(str_bin, unit if id_len > 0 else -unit)
    if hm == 0:
        return str_channel_name(0, 1)
    if hp_len >= 1:
        pool = "C" if base_char in ("C", "G") else "T"
        return f"HP{hp_len}_len{id_len}_{pool}"
    return None


_REFINE_WORKER_CTX = {}


def init_refine_worker(depth_by_trinuc, depth_by_hpstr):
    """Pool initializer for refine_channel_task: stores the read-only depth
    matrices once per worker instead of pickling them into every
    per-channel task.
    """
    global _REFINE_WORKER_CTX
    _REFINE_WORKER_CTX = dict(
        depth_by_trinuc=depth_by_trinuc,
        depth_by_hpstr=depth_by_hpstr,
        n1_mask=np.minimum(*np.indices((10, 10))) >= 1,
    )


def _channel_eeff_at_threshold(kind, ctx_key):
    """Raw opportunity (Eeff) for one channel: duplex depth (both strands
    present) in the channel's context, with no detection-power weighting.
    Reads the worker context set by init_refine_worker.
    """
    c = _REFINE_WORKER_CTX
    if kind == "sbs96":
        t_fwd, b_fwd, t_rc = ctx_key
        n1_mask = c["n1_mask"]
        eeff = float(np.sum(c["depth_by_trinuc"][:, :, t_fwd] * n1_mask))
        eeff += float(np.sum(c["depth_by_trinuc"][:, :, t_rc] * n1_mask))
        return eeff
    # Indel channels, like SBS, need reads on both strands (n1_mask).
    if kind == "hp":
        # HP channels (id_len +-1) sum the depth of their pool's two bases:
        # "T" pool = A/T, "C" pool = C/G.
        hp_len, id_len, pool = ctx_key
        n1_mask = c["n1_mask"]
        total = 0.0
        for base in (0, 1) if pool == "T" else (2, 3):
            total += float(
                np.sum(c["depth_by_hpstr"][:, :, base * 10 + hp_len - 1] * n1_mask)
            )
        return total
    if kind == "str":
        # str_bin 1-4: STR contexts (whole-unit slips, id_len the signed
        # repeat unit); 0: the mismatched 1bp insertion channel and the
        # >=2bp indels that aren't a slip of an STR (id_len the length).
        str_bin, id_len = ctx_key
        if str_bin == 0:
            bucket = depth_hpstr_str0_bucket(id_len)
        else:
            bucket = depth_hpstr_str_bucket(str_bin, id_len)
        return float(np.sum(c["depth_by_hpstr"][:, :, bucket] * c["n1_mask"]))
    raise ValueError(f"unknown channel kind {kind!r}")


def mu_at_min_lr(fdr_thr, min_log10_lr):
    """The mutation rate whose local-fdr threshold (see _refine_channel)
    equals min_log10_lr: (1-f)/((1-f) + f*10**min_log10_lr), about 1.9e-4
    at f = 0.05, min_log10_lr = 5. The threshold falls as mu rises, so a
    channel's final threshold max(min_log10_lr, threshold(mu0)) is exactly
    threshold(min(mu0, this)); per-call FDRs use that capped mu, so a call
    at its channel's threshold has local fdr fdr_thr."""
    return (1.0 - fdr_thr) / ((1.0 - fdr_thr) + fdr_thr * 10.0**min_log10_lr)


def _refine_channel(mu0, threshold0, fdr_thr, min_log10_lr=0.0):
    """FDR-controlled log10 LR threshold for one channel, taking mu0 as the
    channel's mutation rate: the LR where the local fdr
    (1-mu0)/(LR*mu0+1-mu0) equals fdr_thr, i.e.
    LR = (1-fdr_thr)*(1-mu0)/(fdr_thr*mu0).

    Floored at min_log10_lr (--minLR): mu includes the channel's own
    candidates, so in a low-coverage channel one call could otherwise
    raise mu enough to pass itself. Also floored at threshold0, since
    round 1 never records candidates below it. The same threshold drives
    the coverage pass's power simulation.

    mu0 == 0 (no coverage and no sites, or a -mr entry of 0) returns +inf,
    so nothing in the channel passes.
    """
    if mu0 == 0:
        return float("inf")
    lr_threshold = (1.0 - fdr_thr) / (fdr_thr * mu0 / (1 - mu0))
    return max(threshold0, min_log10_lr, float(np.log10(lr_threshold)))


# Upper end of the mu0 bracket: mu0 = 1 would make _refine_channel's
# mu0/(1-mu0) divide by zero.
MU0_UPPER_GAP = 1e-12


def refine_channel_task(job):
    """Per-channel worker task: compute the channel's raw Eeff, take mu0
    from mu0_override (-mr/--muterateprefix) or solve for it, and return
    the refined threshold (_refine_channel). Eeff0 is always computed from
    this run's coverage, even with an override.

    raw_lr_list: raw (10**LR) LRs of the channel's LR >= 0 sites;
    neg_log_lr: log10 LRs (float32) of its LR < 0 sites.

    Without an override, mu0 is the maximum-likelihood mixture weight, the
    root of the log-likelihood derivative
        score(mu) = sum((r_i-1)/(1-mu+mu*r_i)) - (Eeff0-n)/(1-mu),
    which is strictly decreasing (at most one root in (0, 1)).
    score(0) = sum(r_i) - Eeff0; if that is <= 0 (no more evidence than
    noise, since E[LR] = 1 under the null) the likelihood peaks at 0 and
    mu0 = 1/(Eeff0+1) instead. Its threshold log10(Eeff0*(1-f)/f) exceeds
    every site's LR when f < 0.5, so there are still no PASS calls but the
    coverage pass gets finite power. No coverage and no sites gives
    mu0 = 0 (threshold +inf). n <= Eeff0 should hold (each site is one
    covered family-position); if not, a warning is printed and Eeff0 is
    raised to n. If score stays positive up to the bracket end (possible
    when n == Eeff0), mu0 is just below 1.
    """
    (
        name,
        kind,
        ctx_key,
        raw_lr_list,
        threshold0,
        fdr_thr,
        mu0_override,
        neg_log_lr,
        min_log10_lr,
    ) = job
    Eeff0 = _channel_eeff_at_threshold(kind, ctx_key)

    if mu0_override is not None:
        mu0 = mu0_override
    else:
        raw_lr = np.concatenate(
            [
                np.asarray(raw_lr_list, dtype=float),
                10.0 ** np.asarray(neg_log_lr, dtype=float),
            ]
        )
        n = len(raw_lr)
        if not np.isfinite(raw_lr).all():
            raise ValueError(f"channel {name}: non-finite LR among its mu sites")
        eeff_solve = Eeff0
        if n > Eeff0:
            print(
                f"WARNING: channel {name} has {n} mu sites but Eeff "
                f"{Eeff0:g}; solving its mutation rate with Eeff = {n}"
            )
            eeff_solve = float(n)
        if eeff_solve <= 0:
            mu0 = 0.0
        elif n == 0 or raw_lr.sum() <= eeff_solve:
            mu0 = 1.0 / (eeff_solve + 1.0)
        else:

            def score(mu):
                return np.sum((raw_lr - 1.0) / (1.0 - mu + mu * raw_lr)) - (
                    eeff_solve - n
                ) / (1.0 - mu)

            mu_hi = 1.0 - MU0_UPPER_GAP
            # Relative tolerance only: mu can be ~1e-12, below brentq's
            # default absolute xtol.
            mu0 = brentq(score, 0.0, mu_hi, xtol=1e-300) if score(mu_hi) < 0 else mu_hi
    new_threshold = _refine_channel(mu0, threshold0, fdr_thr, min_log10_lr)
    return name, kind, ctx_key, Eeff0, mu0, new_threshold


# ======================================================================
# Depth-matrix accumulation
# ======================================================================
def _accumulate_depth_matrix(depth_mat, n_top, n_bot, category, valid, n_cat):
    """
    Add one count per valid position to depth_mat[n_top, n_bot, category]
    via a single bincount.

    depth_mat: (10, 10, n_cat) running total, updated in place.
    n_top, n_bot: (window_len,) family top/bottom-strand read counts,
        already capped at 9.
    category: (window_len,) context bucket in [0, n_cat); out-of-range
        entries are fine where `valid` is False.
    valid: (window_len,) bool mask of positions to count.
    """
    if not np.any(valid):
        return
    idx = (
        n_top[valid].astype(np.int64) * 10 * n_cat
        + n_bot[valid].astype(np.int64) * n_cat
        + category[valid].astype(np.int64)
    )
    depth_mat += np.bincount(idx, minlength=10 * 10 * n_cat).reshape(10, 10, n_cat)
