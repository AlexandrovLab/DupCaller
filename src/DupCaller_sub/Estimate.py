#!/usr/bin/env python3
import os
import sys
import time
import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sigProfilerPlotting as sigPlt
from pysam import VariantFile as VCF, TabixFile
from scipy.stats import gamma, barnard_exact
import pysam
from . import __version__
from .funcs.misc import check_h5_usable
from .funcs.misc import log_progress
from .funcs.misc import _ensure_type_subdirs
from .funcs.misc import parse_stats_file
from .funcs.misc import build_trinuc64_order, build_trinuc192_labels
from .funcs.misc import (
    _REVCOMP,
    _build_sbs96_labels,
    TRINUCSBS2NUM_96,
    NUM2TRINUCSBS_96,
    combine_raw192_to_sbs96,
    _safe_rate,
)
from .funcs.misc import _dbs_alt_choices, build_dbs_raw144_labels
from .funcs.misc import (
    build_indel100_labels,
    indel100_reference_bucket_indices,
    classify_indel_channel,
    indel_coverage_category_index,
    INDEL_DEL_HP_LEN_BINS,
    INDEL_INS_HP_LEN_BINS,
    INDEL100_INS_HP_LEN_BINS,
    INDEL_STR_UNIT_BINS,
    INDEL_DEL_STR_COUNT_BINS,
    INDEL_INS_STR_COUNT_BINS,
    INDEL_MH_DEL_LEN_BINS,
)


def calculate_ref_trinuc(args):
    """Reference trinucleotide counts over the 64 raw (un-folded) contexts;
    reverse-complement folding happens in estimate_96.

    -ft/--refTrinuc loads a table written by -ot/--outTrinuc instead of
    rescanning .tn.h5; -ot/--outTrinuc writes the table used.
    """
    _, num2trinuc = build_trinuc64_order()
    if getattr(args, "refTrinuc", None):
        table = pd.read_csv(args.refTrinuc, sep="\t", index_col=0, comment="#")
        if list(table.index) != num2trinuc:
            raise ValueError(
                f"-ft/--refTrinuc file {args.refTrinuc} does not hold the 64 "
                "trinucleotide contexts written by -ot/--outTrinuc"
            )
        with open(args.refTrinuc) as f:
            header = f.readline().strip()
        if header.startswith("# regions:"):
            file_regions = header.split(":", 1)[1].split()
            if file_regions != list(args.regions):
                print(
                    f"WARNING: -ft/--refTrinuc was computed over {file_regions}, "
                    f"but -r/--regions is {list(args.regions)}"
                )
        trinuc_count = table["count"].to_numpy(dtype=float)
    else:
        trinuc_count = _scan_ref_trinuc(args)
    if getattr(args, "outTrinuc", None):
        with open(args.outTrinuc, "w") as f:
            f.write("# regions: " + " ".join(args.regions) + "\n")
            pd.DataFrame(
                {"count": trinuc_count}, index=pd.Index(num2trinuc, name="trinuc")
            ).to_csv(f, sep="\t")
    return trinuc_count


def _scan_ref_trinuc(args):
    tn_int = h5py.File(args.reference + ".tn.h5", "r")
    trinuc_count = np.zeros(97)
    n_chroms = len(args.regions)
    progress = {"start": time.time(), "last": time.time()}
    for i, chrom in enumerate(args.regions):
        log_progress(
            progress,
            "Computing reference trinucleotide composition",
            i,
            n_chroms,
            extra=chrom,
        )
        tn_vals = tn_int[chrom][:]
        trinuc_count += np.bincount(tn_vals, minlength=97)
    return trinuc_count[0:64]


# base2num convention (A=0,T=1,C=2,G=3); values outside 0-3 decode to "N",
# which never matches in classify_indel_channel's sequence comparisons.
_NUM2BASE = "ATCG"


def _decode_ref_seq(int_arr):
    """Decode a ref.h5 base2num-encoded array into an uppercase base string."""
    return "".join(_NUM2BASE[b] if 0 <= b <= 3 else "N" for b in int_arr)


# Minimum reference flank on each side of an indel for ID83 classification;
# widened to 6 event lengths when longer, so the repeat-count scan always
# reaches the capped 5+/6+ bins.
INDEL_CONTEXT_WINDOW = 50


def _indel83_record_label(ref_h5, rec):
    """ID83 label of one indel VCF record, plus the (indel_seq, ref_after,
    indel_len) it was resolved from. rec.pos (1-based anchor) is the 0-based
    index of the base after the anchor, so ref_before ends at rec.pos - 1."""
    indel_len = len(rec.alts[0]) - len(rec.ref)
    chrom_ref = ref_h5[rec.chrom]
    chrom_len = chrom_ref.shape[0]
    flank = max(INDEL_CONTEXT_WINDOW, 6 * abs(indel_len))
    if indel_len < 0:
        indel_seq = _decode_ref_seq(chrom_ref[rec.pos : rec.pos - indel_len])
        after_start = rec.pos - indel_len
    else:
        indel_seq = rec.alts[0][1:].upper()
        after_start = rec.pos
    ref_before = _decode_ref_seq(chrom_ref[max(0, rec.pos - flank) : rec.pos])
    ref_after = _decode_ref_seq(
        chrom_ref[after_start : min(chrom_len, after_start + flank)]
    )
    label = classify_indel_channel(indel_seq, ref_before, ref_after, indel_len)
    return label, indel_seq, ref_after, indel_len


def calculate_ref_indel100(args):
    """Genome-wide indel100 reference opportunity (unweighted position
    counts per category, per indel100_reference_bucket_indices), the indel
    analogue of calculate_ref_trinuc.

    hp.h5 and str.h5 are read raw (not via load_repeat_context's
    STR-priority merge) so homopolymer and STR opportunity stay independent;
    a position can validly contribute to both.
    """
    hp_h5 = h5py.File(args.reference + ".hp.h5", "r")
    str_h5 = h5py.File(args.reference + ".str.h5", "r")
    ref_h5 = h5py.File(args.reference + ".ref.h5", "r")
    ref_indel100 = np.zeros(100)
    n_chroms = len(args.regions)
    progress = {"start": time.time(), "last": time.time()}
    for i, chrom in enumerate(args.regions):
        log_progress(
            progress,
            "Computing reference indel100 composition",
            i,
            n_chroms,
            extra=chrom,
        )
        hp_run_arr = hp_h5[chrom][0, :]
        hp_cut_arr = hp_h5[chrom][1, :]
        str_unit_len_arr = str_h5[chrom][0, :]
        str_repeat_count_arr = str_h5[chrom][1, :]
        str_cut_arr = str_h5[chrom][2, :]
        ref_base_arr = ref_h5[chrom][:]
        # Next reference base per position; -1 sentinel for the last
        # position (hence the signed dtype).
        next_ref_base_arr = np.empty(ref_base_arr.shape[0], dtype=np.int16)
        next_ref_base_arr[:-1] = ref_base_arr[1:]
        next_ref_base_arr[-1] = -1
        ref_indel100 += indel100_reference_bucket_indices(
            hp_run_arr,
            str_unit_len_arr,
            str_repeat_count_arr,
            ref_base_arr,
            next_ref_base_arr,
            hp_cut_arr,
            str_cut_arr,
        )
    return ref_indel100


def redistribute_indel_cov16_to_100(
    cov16_row, hp_run, str_unit_len, str_repeat_count, ref_base, label2num_100, out100
):
    """Add one locus's 16 coverage.bed.gz indel opportunity values
    (INDEL_COVERAGE_CATEGORY_LABELS order) into their indel100 bin(s),
    chosen by this locus's own hp.h5/str.h5/ref.h5 annotation. The stored
    values are already gated as in the main pipeline, so no eligibility
    check is needed here. Used by the -rb/--reestimatebed path.

    Column -> destination:
      0      del_unit       -> delSTR{unit}_rep{count}
      1      ins_unit       -> insSTR{unit}_rep{count}, plus _rep0 if count==1
      2-5    del_2..5+      -> both delMH_len{L} and delSTR{L}_rep1
      6-9    ins_A/T/C/G    -> 1bpins{base}
      10-13  ins_len2..5+   -> both insSTR{L}_rep0 and _rep1
      14     del_hp         -> del{ref_base}_hp{run}
      15     ins_hp         -> ins{ref_base}_hp{run} (capped 5+, no hp0 bin)

    hp_run, str_unit_len, str_repeat_count, ref_base: plain ints from
    hp.h5 row 0, str.h5 rows 0/1, ref.h5. Mutates out100 (100,) in place.
    """
    (
        del_unit,
        ins_unit,
        del2,
        del3,
        del4,
        del5p,
        insA,
        insT,
        insC,
        insG,
        insl2,
        insl3,
        insl4,
        insl5p,
        del_hp,
        ins_hp,
    ) = cov16_row

    unit_label = "5+" if str_unit_len >= 5 else str(max(int(str_unit_len), 2))
    del_count_label = (
        "6+" if str_repeat_count >= 6 else str(max(int(str_repeat_count), 1))
    )
    ins_count_label = (
        "5+" if str_repeat_count >= 5 else str(max(int(str_repeat_count), 0))
    )

    out100[label2num_100[f"delSTR{unit_label}_rep{del_count_label}"]] += del_unit
    out100[label2num_100[f"insSTR{unit_label}_rep{ins_count_label}"]] += ins_unit
    if str_repeat_count == 1:
        out100[label2num_100[f"insSTR{unit_label}_rep0"]] += ins_unit

    # INDEL_STR_UNIT_BINS and INDEL_MH_DEL_LEN_BINS share the same
    # ["2","3","4","5+"] bins, so one zip covers both destinations.
    for unit, v in zip(INDEL_STR_UNIT_BINS, (del2, del3, del4, del5p)):
        out100[label2num_100[f"delMH_len{unit}"]] += v
        out100[label2num_100[f"delSTR{unit}_rep1"]] += v

    for base, v in zip("ATCG", (insA, insT, insC, insG)):
        out100[label2num_100[f"1bpins{base}"]] += v

    for unit, v in zip(INDEL_STR_UNIT_BINS, (insl2, insl3, insl4, insl5p)):
        out100[label2num_100[f"insSTR{unit}_rep0"]] += v
        out100[label2num_100[f"insSTR{unit}_rep1"]] += v

    if 0 <= ref_base <= 3:
        pool_base = _NUM2BASE[ref_base]
        del_hp_label = "6+" if hp_run >= 6 else str(max(int(hp_run), 1))
        ins_hp_label = "5+" if hp_run >= 5 else str(max(int(hp_run), 1))
        out100[label2num_100[f"del{pool_base}_hp{del_hp_label}"]] += del_hp
        out100[label2num_100[f"ins{pool_base}_hp{ins_hp_label}"]] += ins_hp


def fay_feuer_confint(mut_counts, correction_ratio, cov, alpha=0.05):
    """Fay-Feuer (1997) gamma CI for a correction-ratio-weighted burden.

    Each channel count is its own Poisson MLE, so Y = sum(w_i*count_i) has
    mean y = sum(w*count) and variance v = sum(w^2*count); the bounds come
    from gamma(shape=y^2/v, scale=v/y), with the upper bound using one extra
    event in the largest-weight channel (y+w_max, v+w_max^2) so it stays
    nonzero at y == 0. With all weights 1 this reduces exactly to the
    classic chi-square Poisson interval, so callers pass (k, 1.0, cov) for
    an unweighted burden.

    mut_counts: (n_channels,) counts, or a scalar total k.
    correction_ratio: (n_channels,) weights, or scalar 1.0.
    cov: coverage denominator of the point estimate; bounds are divided by it.
    """
    if cov <= 0:
        return float("nan"), float("nan")
    mut_counts = np.asarray(mut_counts, dtype=float)
    correction_ratio = np.asarray(correction_ratio, dtype=float)
    y = np.sum(correction_ratio * mut_counts)
    v = np.sum(correction_ratio**2 * mut_counts)
    w_max = float(correction_ratio.max()) if correction_ratio.size else 0.0

    low = gamma.ppf(alpha / 2, a=y**2 / v, scale=v / y) if y > 0 else 0.0

    y_star = y + w_max
    v_star = v + w_max**2
    high = (
        gamma.ppf(1 - alpha / 2, a=y_star**2 / v_star, scale=v_star / y_star)
        if y_star > 0
        else 0.0
    )

    return float(low / cov), float(high / cov)


def _dinuc_num(dinuc):
    """Encode a 2-base string as 4*base2num[b0] + base2num[b1] (A=0,T=1,C=2,G=3)."""
    b2n = {"A": 0, "T": 1, "C": 2, "G": 3}
    return 4 * b2n[dinuc[0]] + b2n[dinuc[1]]


def _dinuc_revcomp(dinuc):
    return _REVCOMP[dinuc[1]] + _REVCOMP[dinuc[0]]


# The 10 canonical DBS78 reference dinucleotides: a fixed COSMIC/
# SigProfiler convention, not derivable from a simple rule (e.g. the AA/TT
# pair's representative is "TT").
DBS_CANONICAL_REFS = ["AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"]


def _build_dbs_raw_to_canonical_ref():
    """All 16 raw dinucleotides -> their canonical-ref representative."""
    mapping = {}
    for ref in DBS_CANONICAL_REFS:
        mapping[ref] = ref
        mapping[_dinuc_revcomp(ref)] = ref
    assert len(mapping) == 16
    return mapping


DBS_RAW_TO_CANONICAL_REF = _build_dbs_raw_to_canonical_ref()

# For each self-complementary ref (AT, TA, CG, GC), the 6 of its 9 raw alts
# that are the canonical representatives. Taken verbatim from the COSMIC
# DBS78 catalog (as used by SigProfilerPlotting); there is no simple rule,
# and matching it keeps our labels interoperable with COSMIC signatures.
_DBS78_SELFCOMP_CANONICAL_ALTS = {
    "AT": ["CA", "CC", "CG", "GA", "GC", "TA"],
    "TA": ["AT", "CG", "CT", "GC", "GG", "GT"],
    "CG": ["AT", "GC", "GT", "TA", "TC", "TT"],
    "GC": ["AA", "AG", "AT", "CA", "CG", "TA"],
}


def build_dbs78_labels():
    """78 canonical DBS labels "REF>ALT" (deterministic order), plus the
    (canonical_ref, raw_alt) -> canonical_alt map.

    78 rather than 10*9: AT, TA, CG and GC are their own reverse complement,
    so their 9 alts pair up among themselves (e.g. AT>CA == AT>TG viewed
    from the other strand), leaving 3 self-inverse + 3 merged = 6 channels
    each: 4*6 + 6*9 = 78 (COSMIC DBS78).
    """
    label2num = {}
    labels = []
    alt_canonical = {}

    for ref in DBS_CANONICAL_REFS:
        alts = _dbs_alt_choices(ref)
        if _dinuc_revcomp(ref) == ref:
            # Self-complementary ref: exactly one of alt / RC(alt) is in
            # the canonical set.
            canonical_alts = _DBS78_SELFCOMP_CANONICAL_ALTS[ref]
            canonical_set = set(canonical_alts)
            for alt in alts:
                if alt in canonical_set:
                    alt_canonical[(ref, alt)] = alt
                else:
                    rc_alt = _dinuc_revcomp(alt)
                    assert rc_alt in canonical_set
                    alt_canonical[(ref, alt)] = rc_alt
            for alt in canonical_alts:
                label = f"{ref}>{alt}"
                label2num[label] = len(labels)
                labels.append(label)
        else:
            for alt in alts:
                alt_canonical[(ref, alt)] = alt
                label = f"{ref}>{alt}"
                label2num[label] = len(labels)
                labels.append(label)

    assert len(labels) == 78
    return label2num, labels, alt_canonical


DBS78_LABEL2NUM, DBS78_LABELS, _DBS78_ALT_CANONICAL = build_dbs78_labels()


def dbs_raw_event_to_label(raw_ref, raw_alt):
    """Fold one raw (reference-strand) dinucleotide substitution onto its
    DBS78 label."""
    canonical_ref = DBS_RAW_TO_CANONICAL_REF[raw_ref]
    eff_alt = raw_alt if canonical_ref == raw_ref else _dinuc_revcomp(raw_alt)
    canonical_alt = _DBS78_ALT_CANONICAL[(canonical_ref, eff_alt)]
    return f"{canonical_ref}>{canonical_alt}"


DBS_RAW144_LABEL2NUM, DBS_RAW144_LABELS = build_dbs_raw144_labels()


def combine_raw_dbs_to_dbs78(values_144, label2num_144):
    """Fold a raw-144 DBS array to the 78 DBS78 classes, summing each raw
    (ref,alt) row into dbs_raw_event_to_label's class. The number of raw
    rows per class varies (self-complementary refs), hence the per-label
    lookup.

    values_144   : array shaped (144,) or (144, n_groups).
    label2num_144: {"{ref}>{alt}": row_index} for values_144's rows.
    Returns an array shaped (78,) or (78, n_groups).
    """
    out_shape = (78,) if values_144.ndim == 1 else (78, values_144.shape[1])
    combined = np.zeros(out_shape)
    for raw_label, raw_idx in label2num_144.items():
        raw_ref, raw_alt = raw_label.split(">")
        canonical_label = dbs_raw_event_to_label(raw_ref, raw_alt)
        combined[DBS78_LABEL2NUM[canonical_label]] += values_144[raw_idx]
    return combined


def calculate_ref_dbs(args):
    """Genome-wide reference dinucleotide counts (16 raw contexts), computed
    from ref.h5 by pairing each position with its successor. Folding to
    DBS78 happens in _dbs78_ref_weighted.
    """
    ref_h5 = h5py.File(args.reference + ".ref.h5", "r")
    dinuc_count = np.zeros(16)
    n_chroms = len(args.regions)
    progress = {"start": time.time(), "last": time.time()}
    for i, chrom in enumerate(args.regions):
        log_progress(
            progress,
            "Computing reference dinucleotide composition",
            i,
            n_chroms,
            extra=chrom,
        )
        ref_base_arr = ref_h5[chrom][:]
        if ref_base_arr.shape[0] < 2:
            continue
        base1 = ref_base_arr[:-1]
        base2 = ref_base_arr[1:]
        # Exclude any pair touching an N/out-of-range base (encoded as 4).
        valid = (base1 <= 3) & (base2 <= 3)
        dinuc_idx = 4 * base1[valid] + base2[valid]
        dinuc_count += np.bincount(dinuc_idx, minlength=16)
    return dinuc_count


def calculate_dbs_opportunity(coverage_bed_path, ref_h5, regions):
    """Raw-144 DBS opportunity approximated from the per-locus SBS coverage
    bed. Not used by do_estimate (which reads _dbs_by_duplex_group.txt);
    kept as a standalone cross-check.

    For each adjacent covered pair (start, start+1) with A/T/C/G reference
    at both, opportunity for alt pair (b1, b2) is cov(b1 at pos1) *
    cov(b2 at pos2), i.e. assumes the two positions are independent. Masks
    baked into the coverage bed are inherited.

    coverage_bed_path: path to {sample}_coverage.bed.gz.
    ref_h5: open h5py.File for {reference}.ref.h5.
    regions: chromosome names to process (args.regions).
    """
    tbx = TabixFile(coverage_bed_path)
    bases = "ATCG"
    opp144 = np.zeros(144)

    n_chroms = len(regions)
    progress = {"start": time.time(), "last": time.time()}
    for ci, chrom in enumerate(regions):
        log_progress(
            progress,
            "Computing DBS opportunity from coverage",
            ci,
            n_chroms,
            extra=chrom,
        )
        if chrom not in tbx.contigs:
            continue
        starts = []
        covs = []
        for line in tbx.fetch(chrom):
            parts = line.split("\t", 7)
            if len(parts) < 7:
                continue
            starts.append(int(parts[1]))
            covs.append(
                (float(parts[3]), float(parts[4]), float(parts[5]), float(parts[6]))
            )
        if len(starts) < 2:
            continue
        starts = np.array(starts, dtype=np.int64)
        covs = np.array(covs, dtype=float)  # (n, 4): A,T,C,G coverage columns

        adjacent = starts[1:] == starts[:-1] + 1
        if not np.any(adjacent):
            continue

        ref_arr = ref_h5[chrom][:]
        pos1 = starts[:-1][adjacent]
        ref1 = ref_arr[pos1]
        ref2 = ref_arr[pos1 + 1]
        cov1 = covs[:-1][adjacent]  # (k, 4)
        cov2 = covs[1:][adjacent]  # (k, 4)

        # 4x4 ref combos (N, encoded as 4, never matches) x 9 alt combos,
        # each vectorized over all adjacent pairs on the chromosome.
        for ref1_num in range(4):
            m1 = ref1 == ref1_num
            if not np.any(m1):
                continue
            for ref2_num in range(4):
                m2 = m1 & (ref2 == ref2_num)
                if not np.any(m2):
                    continue
                ref_dinuc = bases[ref1_num] + bases[ref2_num]
                c1 = cov1[m2]
                c2 = cov2[m2]
                for b1_num in range(4):
                    if b1_num == ref1_num:
                        continue
                    for b2_num in range(4):
                        if b2_num == ref2_num:
                            continue
                        alt_dinuc = bases[b1_num] + bases[b2_num]
                        contribution = float(np.sum(c1[:, b1_num] * c2[:, b2_num]))
                        raw_idx = DBS_RAW144_LABEL2NUM[f"{ref_dinuc}>{alt_dinuc}"]
                        opp144[raw_idx] += contribution
    return opp144


_INDEL_HP_BASE_POOLS = (("CG", ("C", "G")), ("AT", ("A", "T")))


def build_indel76_labels():
    """74 strand-folded indel labels (the name "indel76" is historical),
    the folded counterpart of build_indel100_labels.

    Homopolymer channels pool C+G and A+T (a C run on one strand is a G run
    on the other); STR and microhomology channels pass through. Insertion
    hp has no rep0 bin (INDEL100_INS_HP_LEN_BINS).
    6+6 (del hp) + 5+5 (ins hp) + 24 (STR del) + 24 (STR ins) + 4 (MH) = 74.
    Deterministic order; must not be built independently elsewhere.
    """
    labels = []
    label2num = {}

    def add(label):
        label2num[label] = len(labels)
        labels.append(label)

    for pool, _ in _INDEL_HP_BASE_POOLS:
        for length in INDEL_DEL_HP_LEN_BINS:
            add(f"del{pool}_hp{length}")
    for pool, _ in _INDEL_HP_BASE_POOLS:
        for length in INDEL100_INS_HP_LEN_BINS:
            add(f"ins{pool}_hp{length}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_DEL_STR_COUNT_BINS:
            add(f"delSTR{unit}_rep{count}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_INS_STR_COUNT_BINS:
            add(f"insSTR{unit}_rep{count}")
    for length in INDEL_MH_DEL_LEN_BINS:
        add(f"delMH_len{length}")

    return label2num, labels


def combine_indel100_to_indel76(values_100, label2num_100):
    """Fold the 100 raw indel classes into the 74 indel76 classes by
    summing homopolymer C+G / A+T pairs; STR and MH pass through.

    values_100    : array shaped (100,) or (100, n_groups).
    label2num_100 : {"label": row_index} for values_100's rows, as written
        by Caller.py / build_indel100_labels.
    Returns an array shaped (74,) or (74, n_groups).
    """
    _, labels_76 = build_indel76_labels()
    out_shape = (74,) if values_100.ndim == 1 else (74, values_100.shape[1])
    combined = np.zeros(out_shape)
    idx = 0
    for prefix, len_bins in (
        ("del", INDEL_DEL_HP_LEN_BINS),
        ("ins", INDEL100_INS_HP_LEN_BINS),
    ):
        for _, bases in _INDEL_HP_BASE_POOLS:
            for length in len_bins:
                for base in bases:
                    combined[idx] += values_100[
                        label2num_100[f"{prefix}{base}_hp{length}"]
                    ]
                idx += 1
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_DEL_STR_COUNT_BINS:
            combined[idx] = values_100[label2num_100[f"delSTR{unit}_rep{count}"]]
            idx += 1
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_INS_STR_COUNT_BINS:
            combined[idx] = values_100[label2num_100[f"insSTR{unit}_rep{count}"]]
            idx += 1
    for length in INDEL_MH_DEL_LEN_BINS:
        combined[idx] = values_100[label2num_100[f"delMH_len{length}"]]
        idx += 1
    assert idx == 74
    return combined


def _safe_correction_ratio(ref_counts, obs_counts):
    """Per-class (ref fraction) / (observed coverage fraction). Classes with
    zero observed coverage get 0 rather than Inf, since Inf * 0 mutations
    would be NaN and poison the corrected-burden sum; an all-zero stratum
    gets all zeros."""
    obs_sum = obs_counts.sum()
    if obs_sum <= 0:
        return np.zeros_like(obs_counts, dtype=float)
    ref_frac = ref_counts / ref_counts.sum()
    obs_frac = obs_counts / obs_sum
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(obs_frac > 0, ref_frac / obs_frac, 0.0)


def estimate_96(trinuc_cov_96_by_rf, trinuc_mut_by_rf, ref_trinuc_64, n):
    """Estimate SBS96 mutation rates with per-class correction ratios.

    trinuc_cov_96_by_rf : (96, n_groups) L-weighted coverage per SBS96 class.
    trinuc_mut_by_rf : (96, n_groups) observed mutation counts per SBS96 class.
    ref_trinuc_64    : (64,) raw reference trinucleotide counts; folded to 32
        contexts here via the i/i+32 reverse-complement pairing.

    Returns burden curves both cumulative (nmin >= 1..5) and exact
    (nmin == 1..4, plus pooled nmin >= 5; the _exact return values), where
    nmin = min(F1R2, F2R1) of the duplex group.
    """
    print("........Estimating mutation rate for each trinucleotide context.......")
    ref_trinuc_32 = ref_trinuc_64[:32] + ref_trinuc_64[32:64]

    ref_trinuc_96 = np.repeat(
        ref_trinuc_32, 3
    )  # (96,) same ref count for the 3 alts of each trinuc

    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    burden_uncorrected = np.zeros(5)
    burden_uncorrected_ub = np.zeros(5)
    burden_uncorrected_lb = np.zeros(5)
    burden_corrected = np.zeros(5)
    burden_corrected_ub = np.zeros(5)
    burden_corrected_lb = np.zeros(5)
    covs = np.zeros(5)
    hap_trinuc = np.zeros([96, 5])

    ## Calculate burden when minimum strand depth >= 5
    trinuc_mut = trinuc_mut_by_rf[:, nmin >= 5].sum(axis=1)  # (96,)
    trinuc_cov = trinuc_cov_96_by_rf[:, nmin >= 5].sum(
        axis=1
    )  # (96,) alt-base-specific
    mutnum = trinuc_mut.sum()
    cov = trinuc_cov.sum() / 3  # divide by 3 to get per-locus-equivalent coverage
    burden_uncorrected[4] = mutnum / cov if cov > 0 else 0.0
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = fay_feuer_confint(
        mutnum, 1.0, cov
    )
    trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov)
    correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov)
    mutnum_corrected = correction_ratio * trinuc_mut
    burden_corrected[4] = mutnum_corrected.sum() / cov if cov > 0 else 0.0
    burden_corrected_lb[4], burden_corrected_ub[4] = fay_feuer_confint(
        trinuc_mut, correction_ratio, cov
    )
    hap_trinuc[:, 4] = trinuc_rate * ref_trinuc_96
    covs[4] = cov

    for nn in range(4, 0, -1):
        trinuc_mut = trinuc_mut + trinuc_mut_by_rf[:, nmin == nn].sum(axis=1)
        trinuc_cov = trinuc_cov + trinuc_cov_96_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = trinuc_mut.sum()
        cov = trinuc_cov.sum() / 3
        burden_uncorrected[nn - 1] = mutnum / cov if cov > 0 else 0.0
        (
            burden_uncorrected_lb[nn - 1],
            burden_uncorrected_ub[nn - 1],
        ) = fay_feuer_confint(mutnum, 1.0, cov)
        trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov)
        covs[nn - 1] = cov
        correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov)
        mutnum_corrected = trinuc_mut * correction_ratio
        burden_corrected[nn - 1] = mutnum_corrected.sum() / cov if cov > 0 else 0.0
        (
            burden_corrected_lb[nn - 1],
            burden_corrected_ub[nn - 1],
        ) = fay_feuer_confint(trinuc_mut, correction_ratio, cov)
        hap_trinuc[:, nn - 1] = trinuc_rate * ref_trinuc_96

    # Exact (non-cumulative) group-size bins: nmin == 1..4, plus pooled
    # nmin >= 5; each bin is computed independently.
    burden_uncorrected_exact = np.zeros(5)
    burden_uncorrected_exact_lb = np.zeros(5)
    burden_uncorrected_exact_ub = np.zeros(5)
    burden_corrected_exact = np.zeros(5)
    burden_corrected_exact_lb = np.zeros(5)
    burden_corrected_exact_ub = np.zeros(5)
    covs_exact = np.zeros(5)
    for nn in range(1, 6):
        sel = (nmin == nn) if nn < 5 else (nmin >= 5)
        trinuc_mut_e = trinuc_mut_by_rf[:, sel].sum(axis=1)
        trinuc_cov_e = trinuc_cov_96_by_rf[:, sel].sum(axis=1)
        mutnum_e = trinuc_mut_e.sum()
        cov_e = trinuc_cov_e.sum() / 3
        covs_exact[nn - 1] = cov_e
        burden_uncorrected_exact[nn - 1] = mutnum_e / cov_e if cov_e > 0 else 0.0
        (
            burden_uncorrected_exact_lb[nn - 1],
            burden_uncorrected_exact_ub[nn - 1],
        ) = fay_feuer_confint(mutnum_e, 1.0, cov_e)
        correction_ratio_e = _safe_correction_ratio(ref_trinuc_96, trinuc_cov_e)
        mutnum_corrected_e = (trinuc_mut_e * correction_ratio_e).sum()
        burden_corrected_exact[nn - 1] = (
            mutnum_corrected_e / cov_e if cov_e > 0 else 0.0
        )
        (
            burden_corrected_exact_lb[nn - 1],
            burden_corrected_exact_ub[nn - 1],
        ) = fay_feuer_confint(trinuc_mut_e, correction_ratio_e, cov_e)

    # Zero the count of any channel with zero coverage (its correction_ratio
    # is already 0), so mutation_number_uncorrected never reports mutations
    # without coverage. Fill with int 0, not 0.0: a float fill would upcast
    # the int array and break consumers that parse an integer total (e.g.
    # Summarize.py's "Uncorrected mutation number").
    trinuc_mut = np.where(trinuc_cov > 0, trinuc_mut, 0)

    return (
        trinuc_mut,
        mutnum_corrected,
        correction_ratio,
        hap_trinuc[:, 0],
        burden_corrected,
        burden_corrected_lb,
        burden_corrected_ub,
        burden_uncorrected,
        burden_uncorrected_lb,
        burden_uncorrected_ub,
        ref_trinuc_64.sum(),
        covs,
        trinuc_cov,
        burden_corrected_exact,
        burden_corrected_exact_lb,
        burden_corrected_exact_ub,
        burden_uncorrected_exact,
        burden_uncorrected_exact_lb,
        burden_uncorrected_exact_ub,
        covs_exact,
    )


# Microhomology-length sub-channels per ID83 deletion-length group
# ("2","3","4","5+"): mh_len < del_len, so {1}, {1,2}, {1,2,3},
# {1,2,3,4,5+} -> 11 channels.
INDEL83_MH_GROUP_SIZES = [1, 2, 3, 5]


def build_indel83_labels():
    """83-channel COSMIC ID83 labels (internal syntax), the label set
    classify_indel_channel resolves calls to:
      Cdelhp1-6+, Tdelhp1-6+ (6+6): homopolymer deletion, base pooled
        C+G->C / A+T->T, by existing homopolymer length.
      Cinshp0-5+, Tinshp0-5+ (6+6): homopolymer insertion, same pooling,
        by existing homopolymer length (0 = no repeat there yet).
      {2,3,4,5+}delstr{1-6+} (24): STR deletion by unit size x existing
        repeat count.
      {2,3,4,5+}insstr{0-5+} (24): STR insertion by unit size x existing
        repeat count.
      {2,3,4,5+}delMH{1..} (11): microhomology deletion by deletion length
        x microhomology length (INDEL83_MH_GROUP_SIZES sub-channels per
        deletion-length group).
    6+6+6+6+24+24+11 = 83. Deterministic order — must not be built
    independently elsewhere.
    """
    labels = []
    label2num = {}

    def add(label):
        label2num[label] = len(labels)
        labels.append(label)

    for base in ("C", "T"):
        for length in INDEL_DEL_HP_LEN_BINS:
            add(f"{base}delhp{length}")
    for base in ("C", "T"):
        for length in INDEL_INS_HP_LEN_BINS:
            add(f"{base}inshp{length}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_DEL_STR_COUNT_BINS:
            add(f"{unit}delstr{count}")
    for unit in INDEL_STR_UNIT_BINS:
        for count in INDEL_INS_STR_COUNT_BINS:
            add(f"{unit}insstr{count}")
    for del_group, n_mh in zip(["2", "3", "4", "5+"], INDEL83_MH_GROUP_SIZES):
        for mh in range(1, n_mh + 1):
            mh_label = "5+" if (del_group == "5+" and mh == 5) else str(mh)
            add(f"{del_group}delMH{mh_label}")

    assert len(labels) == 83
    return label2num, labels


def build_indel83_sigprofiler_labels():
    """Map internal ID83 labels (e.g. "Cdelhp1", "2delstr3") to
    SigProfilerPlotting's labels (e.g. "1:Del:C:0", "2:Del:R:2"); same
    channel order. Mirrors build_indel83_labels' loops rather than
    hardcoding a table.

    Returns {our_label: sigprofiler_label}, 83 entries.
    """
    mapping = {}
    for base in ("C", "T"):
        for i, length in enumerate(INDEL_DEL_HP_LEN_BINS):
            mapping[f"{base}delhp{length}"] = f"1:Del:{base}:{i}"
    for base in ("C", "T"):
        for i, length in enumerate(INDEL_INS_HP_LEN_BINS):
            mapping[f"{base}inshp{length}"] = f"1:Ins:{base}:{i}"
    for unit in INDEL_STR_UNIT_BINS:
        unit_num = unit.rstrip("+")
        for i, count in enumerate(INDEL_DEL_STR_COUNT_BINS):
            mapping[f"{unit}delstr{count}"] = f"{unit_num}:Del:R:{i}"
    for unit in INDEL_STR_UNIT_BINS:
        unit_num = unit.rstrip("+")
        for i, count in enumerate(INDEL_INS_STR_COUNT_BINS):
            mapping[f"{unit}insstr{count}"] = f"{unit_num}:Ins:R:{i}"
    for del_group, n_mh in zip(["2", "3", "4", "5+"], INDEL83_MH_GROUP_SIZES):
        group_num = del_group.rstrip("+")
        for mh in range(1, n_mh + 1):
            mh_label = "5+" if (del_group == "5+" and mh == 5) else str(mh)
            mapping[f"{del_group}delMH{mh_label}"] = f"{group_num}:Del:M:{mh}"
    assert len(mapping) == 83
    return mapping


INDEL83_TO_SIGPROFILER_LABELS = build_indel83_sigprofiler_labels()


def expand_indel76_to_indel83(values_76):
    """Expand a 74-class indel76 array to the ID83 layout.

    Coverage has no microhomology-length resolution, so each
    delMH_len{2,3,4,5+} value is replicated across its 1/2/3/5 MH
    sub-channels. Del-hp, STR del and STR ins pass through. indel76 has no
    ins-hp rep0 bin, so Cinshp0 (12) and Tinshp0 (18) are left at 0; the
    caller must fill them with override_inshp0_with_next_base_opportunity.

    values_76: array shaped (74,) or (74, n_groups).
    Returns an array shaped (83,) or (83, n_groups), with indices 12 and
    18 left at 0.
    """
    out_shape = (83,) if values_76.ndim == 1 else (83, values_76.shape[1])
    out = np.zeros(out_shape)
    out[:12] = values_76[:12]  # del hp (CG/AT pools x len 1-6+)
    out[13:18] = values_76[12:17]  # ins hp, CG pool len 1..5+ -> Cinshp1..5+
    out[19:24] = values_76[17:22]  # ins hp, AT pool len 1..5+ -> Tinshp1..5+
    out[24:72] = values_76[22:70]  # STR del (24) + STR ins (24), unchanged
    idx = 72
    for group_idx, n_mh in enumerate(INDEL83_MH_GROUP_SIZES):
        out[idx : idx + n_mh] = values_76[70 + group_idx]
        idx += n_mh
    assert idx == 83
    return out


def override_inshp0_with_next_base_opportunity(values_83, values_100):
    """Fill ID83 Cinshp0/Tinshp0 (indices 12/18) from the raw indel100
    1bpins{A,T,C,G} columns (96-99), pooled C+G / A+T. A rep0 insertion of
    base b only needs the next reference base to differ from b, which is
    exactly what those columns count (for both reference composition and
    coverage).

    values_83  : (83,) or (83, n_groups); not modified in place.
    values_100 : (100,) or (100, n_groups) raw array at the same resolution.
    Returns a new array shaped like values_83.
    """
    out = values_83.copy()
    out[12] = values_100[98] + values_100[99]  # Cinshp0 = Cins + Gins
    out[18] = values_100[96] + values_100[97]  # Tinshp0 = Ains + Tins
    return out


def _indel83_mh_group_slices():
    slices = []
    idx = 72
    for n_mh in INDEL83_MH_GROUP_SIZES:
        slices.append(slice(idx, idx + n_mh))
        idx += n_mh
    return slices


def _grouped_indel83_correction_ratio(ref_frac, obs_frac):
    """ID83 correction ratio ref_frac/obs_frac: per channel for 0-71, but
    one shared ratio per delMH deletion-length group for 72-82, matching
    the coverage resolution. Zero observed fraction gives ratio 0 (see
    _safe_correction_ratio).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(obs_frac > 0, ref_frac / obs_frac, 0.0)
    for sl in _indel83_mh_group_slices():
        obs_group = obs_frac[sl].sum()
        ratio[sl] = ref_frac[sl].sum() / obs_group if obs_group > 0 else 0.0
    return ratio


def estimate_indel83(indel83_cov_by_rf, indel83_mut_by_rf, ref_indel83, n):
    """ID83 counterpart of estimate_96 (same cumulative/exact nmin curves),
    using _grouped_indel83_correction_ratio. Coverage is per opportunity
    channel, not per locus; callers rescale by indel_locus_multiplier.
    """
    print("........Estimating mutation rate for each indel context (ID83).......")
    ref_frac = ref_indel83 / ref_indel83.sum()

    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    burden_uncorrected = np.zeros(5)
    burden_uncorrected_ub = np.zeros(5)
    burden_uncorrected_lb = np.zeros(5)
    burden_corrected = np.zeros(5)
    burden_corrected_ub = np.zeros(5)
    burden_corrected_lb = np.zeros(5)
    covs = np.zeros(5)
    hap_indel83 = np.zeros([83, 5])

    indel_mut = indel83_mut_by_rf[:, nmin >= 5].sum(axis=1)
    indel_cov = indel83_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = indel_mut.sum()
    cov = indel_cov.sum()
    burden_uncorrected[4] = mutnum / cov if cov > 0 else float("nan")
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = fay_feuer_confint(
        mutnum, 1.0, cov
    )
    indel_rate = _safe_rate(indel_mut, indel_cov)
    with np.errstate(divide="ignore", invalid="ignore"):
        obs_frac = indel_cov / indel_cov.sum()
    correction_ratio = _grouped_indel83_correction_ratio(ref_frac, obs_frac)
    mutnum_corrected = correction_ratio * indel_mut
    burden_corrected[4] = mutnum_corrected.sum() / cov if cov > 0 else float("nan")
    burden_corrected_lb[4], burden_corrected_ub[4] = fay_feuer_confint(
        indel_mut, correction_ratio, cov
    )
    hap_indel83[:, 4] = indel_rate * ref_indel83
    covs[4] = cov

    for nn in range(4, 0, -1):
        indel_mut = indel_mut + indel83_mut_by_rf[:, nmin == nn].sum(axis=1)
        indel_cov = indel_cov + indel83_cov_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = indel_mut.sum()
        cov = indel_cov.sum()
        burden_uncorrected[nn - 1] = mutnum / cov if cov > 0 else float("nan")
        (
            burden_uncorrected_lb[nn - 1],
            burden_uncorrected_ub[nn - 1],
        ) = fay_feuer_confint(mutnum, 1.0, cov)
        indel_rate = _safe_rate(indel_mut, indel_cov)
        covs[nn - 1] = cov
        with np.errstate(divide="ignore", invalid="ignore"):
            obs_frac = indel_cov / indel_cov.sum()
        correction_ratio = _grouped_indel83_correction_ratio(ref_frac, obs_frac)
        mutnum_corrected = indel_mut * correction_ratio
        burden_corrected[nn - 1] = (
            mutnum_corrected.sum() / cov if cov > 0 else float("nan")
        )
        (
            burden_corrected_lb[nn - 1],
            burden_corrected_ub[nn - 1],
        ) = fay_feuer_confint(indel_mut, correction_ratio, cov)
        hap_indel83[:, nn - 1] = indel_rate * ref_indel83

    # Exact (non-cumulative) group-size bins, as in estimate_96.
    burden_uncorrected_exact = np.zeros(5)
    burden_uncorrected_exact_lb = np.zeros(5)
    burden_uncorrected_exact_ub = np.zeros(5)
    burden_corrected_exact = np.zeros(5)
    burden_corrected_exact_lb = np.zeros(5)
    burden_corrected_exact_ub = np.zeros(5)
    covs_exact = np.zeros(5)
    for nn in range(1, 6):
        sel = (nmin == nn) if nn < 5 else (nmin >= 5)
        indel_mut_e = indel83_mut_by_rf[:, sel].sum(axis=1)
        indel_cov_e = indel83_cov_by_rf[:, sel].sum(axis=1)
        mutnum_e = indel_mut_e.sum()
        cov_e = indel_cov_e.sum()
        covs_exact[nn - 1] = cov_e
        burden_uncorrected_exact[nn - 1] = mutnum_e / cov_e if cov_e > 0 else 0.0
        (
            burden_uncorrected_exact_lb[nn - 1],
            burden_uncorrected_exact_ub[nn - 1],
        ) = fay_feuer_confint(mutnum_e, 1.0, cov_e)
        with np.errstate(divide="ignore", invalid="ignore"):
            obs_frac_e = indel_cov_e / indel_cov_e.sum()
        correction_ratio_e = _grouped_indel83_correction_ratio(ref_frac, obs_frac_e)
        mutnum_corrected_e = (indel_mut_e * correction_ratio_e).sum()
        burden_corrected_exact[nn - 1] = (
            mutnum_corrected_e / cov_e if cov_e > 0 else 0.0
        )
        (
            burden_corrected_exact_lb[nn - 1],
            burden_corrected_exact_ub[nn - 1],
        ) = fay_feuer_confint(indel_mut_e, correction_ratio_e, cov_e)

    # Zero counts in zero-coverage channels; int fill keeps int dtype
    # (see estimate_96).
    indel_mut = np.where(indel_cov > 0, indel_mut, 0)

    return (
        indel_mut,
        mutnum_corrected,
        correction_ratio,
        hap_indel83[:, 0],
        burden_corrected,
        burden_corrected_lb,
        burden_corrected_ub,
        burden_uncorrected,
        burden_uncorrected_lb,
        burden_uncorrected_ub,
        ref_indel83.sum(),
        covs,
        indel_cov,
        burden_corrected_exact,
        burden_corrected_exact_lb,
        burden_corrected_exact_ub,
        burden_uncorrected_exact,
        burden_uncorrected_exact_lb,
        burden_uncorrected_exact_ub,
        covs_exact,
    )


def _dbs78_ref_weighted(dinuc_count_16):
    """Reference opportunity for the 78 DBS78 channels.

    Each raw ref's genome count is credited to all 9 of its raw alts, then
    folded with combine_raw_dbs_to_dbs78, the same fold used for observed
    opportunity and mutations. A flat per-ref broadcast would be wrong for
    the self-complementary refs, whose merged channels receive 2 raw alts
    but self-inverse ones only 1, and would skew correction_ratio ~2x
    there. No /9 normalization: it would cancel in ref_frac.

    dinuc_count_16: (16,) raw counts from calculate_ref_dbs.
    Returns (78,) canonical-channel reference opportunity.
    """
    raw144 = np.zeros(144)
    for label, idx in DBS_RAW144_LABEL2NUM.items():
        raw_ref = label.split(">")[0]
        raw144[idx] = dinuc_count_16[_dinuc_num(raw_ref)]
    return combine_raw_dbs_to_dbs78(raw144, DBS_RAW144_LABEL2NUM)


def estimate_dbs78(dbs_mut_78, dbs_opp_78, ref_dbs_78):
    """Whole-sample DBS78 burden and per-class correction ratios (see
    estimate_dbs78_by_group for the group-size breakdown).

    dbs_mut_78 : (78,) PASS DBS counts.
    dbs_opp_78 : (78,) opportunity summed over all duplex groups.
    ref_dbs_78 : (78,) reference composition from _dbs78_ref_weighted.
    """
    print("........Estimating mutation rate for each DBS78 context.......")
    # Each raw ref dinucleotide has 9 alt pairs and the 144->78 fold
    # preserves the total, so /9 gives per-locus-pair coverage.
    cov = dbs_opp_78.sum() / 9
    mutnum = dbs_mut_78.sum()
    burden_uncorrected = mutnum / cov if cov > 0 else float("nan")
    burden_uncorrected_lb, burden_uncorrected_ub = fay_feuer_confint(mutnum, 1.0, cov)

    correction_ratio = _safe_correction_ratio(ref_dbs_78, dbs_opp_78)
    mutnum_corrected = correction_ratio * dbs_mut_78
    burden_corrected = mutnum_corrected.sum() / cov if cov > 0 else float("nan")
    burden_corrected_lb, burden_corrected_ub = fay_feuer_confint(
        dbs_mut_78, correction_ratio, cov
    )

    # Zero counts in zero-opportunity channels; int fill keeps int dtype
    # (see estimate_96).
    dbs_mut_78 = np.where(dbs_opp_78 > 0, dbs_mut_78, 0)

    return (
        dbs_mut_78,
        mutnum_corrected,
        correction_ratio,
        burden_corrected,
        burden_corrected_lb,
        burden_corrected_ub,
        burden_uncorrected,
        burden_uncorrected_lb,
        burden_uncorrected_ub,
        cov,
    )


def estimate_dbs78_by_group(dbs_mut_by_rf, dbs_opp144_by_rf, ref_dbs_78, n):
    """DBS78 cumulative and exact group-size burden curves, as in
    estimate_96. Returns the cumulative 7-tuple followed by the exact one.

    dbs_mut_by_rf    : (78, n_groups) PASS DBS counts per duplex group.
    dbs_opp144_by_rf : (144, n_groups) raw opportunity per duplex group;
        each stratum is summed, then folded to 78.
    ref_dbs_78       : (78,) reference composition.
    n                : duplex-group labels ("F1R2+F2R1") in column order,
        which must match trinuc_by_rf.columns.
    """
    print("........Estimating DBS78 mutation rate by duplex group size.......")
    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)

    def _burden_for(sel):
        opp144 = dbs_opp144_by_rf[:, sel].sum(axis=1)
        opp78 = combine_raw_dbs_to_dbs78(opp144, DBS_RAW144_LABEL2NUM)
        mut78 = dbs_mut_by_rf[:, sel].sum(axis=1)
        cov = opp78.sum() / 9
        mutnum = mut78.sum()
        uncorrected = mutnum / cov if cov > 0 else 0.0
        u_lb, u_ub = fay_feuer_confint(mutnum, 1.0, cov)
        correction_ratio = _safe_correction_ratio(ref_dbs_78, opp78)
        mutnum_corrected = (correction_ratio * mut78).sum()
        corrected = mutnum_corrected / cov if cov > 0 else 0.0
        c_lb, c_ub = fay_feuer_confint(mut78, correction_ratio, cov)
        return uncorrected, u_lb, u_ub, corrected, c_lb, c_ub, cov

    def _five_point_curve(sel_fn):
        uncorrected = np.zeros(5)
        uncorrected_lb = np.zeros(5)
        uncorrected_ub = np.zeros(5)
        corrected = np.zeros(5)
        corrected_lb = np.zeros(5)
        corrected_ub = np.zeros(5)
        covs = np.zeros(5)
        for nn in range(1, 6):
            (
                uncorrected[nn - 1],
                uncorrected_lb[nn - 1],
                uncorrected_ub[nn - 1],
                corrected[nn - 1],
                corrected_lb[nn - 1],
                corrected_ub[nn - 1],
                covs[nn - 1],
            ) = _burden_for(sel_fn(nn))
        return (
            uncorrected,
            uncorrected_lb,
            uncorrected_ub,
            corrected,
            corrected_lb,
            corrected_ub,
            covs,
        )

    min_curve = _five_point_curve(lambda nn: nmin >= nn)
    exact_curve = _five_point_curve(lambda nn: (nmin == nn) if nn < 5 else (nmin >= 5))
    return min_curve + exact_curve


def estimate_id(trinuc_cov_by_rf, muts_by_rf, n):
    print("........Estimating indel rate")
    trinuc_mut_cov_by_rf = np.repeat(trinuc_cov_by_rf, 3, axis=0)
    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    burden_indel = np.zeros(5)
    burden_indel_ub = np.zeros(5)
    burden_indel_lb = np.zeros(5)
    trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = muts_by_rf[:, nmin >= 5].sum(axis=1).sum()
    cov = trinuc_cov.sum() / 3
    burden_indel[4] = mutnum / cov
    burden_indel_lb[4], burden_indel_ub[4] = fay_feuer_confint(mutnum, 1.0, cov)
    for nn in range(4, 0, -1):
        trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= nn].sum(axis=1)
        mutnum = muts_by_rf[:, nmin >= nn].sum(axis=1).sum()
        cov = trinuc_cov.sum() / 3
        burden_indel[nn - 1] = mutnum / cov
        burden_indel_lb[nn - 1], burden_indel_ub[nn - 1] = fay_feuer_confint(
            mutnum, 1.0, cov
        )

    return (burden_indel, burden_indel_lb, burden_indel_ub, mutnum)


def write_burden_by_group_size(
    outdir,
    sample,
    mutation_type,
    min_uncorrected,
    min_uncorrected_lb,
    min_uncorrected_ub,
    min_corrected,
    min_corrected_lb,
    min_corrected_ub,
    min_cov,
    exact_uncorrected,
    exact_uncorrected_lb,
    exact_uncorrected_ub,
    exact_corrected,
    exact_corrected_lb,
    exact_corrected_ub,
    exact_cov,
):
    """Write {outdir}/{sample}_{mutation_type}_burden_by_group_size.{txt,pdf}:
    a table with cumulative (min_group_size, nmin >= N) and exact
    (group_size, nmin == N) column blocks, and a 2-panel plot with 95% CIs.
    Every burden/cov argument is a (5,) array over N = 1,2,3,4,5+.

    mutation_type: lowercase file-name infix ("sbs"/"indel"/"dbs").
    """
    read_numbers = np.arange(1, 6, dtype=np.int16)
    table = pd.DataFrame(
        {
            "read_number": read_numbers,
            "min_group_size_Uncorrected_burden": min_uncorrected,
            "min_group_size_Uncorrected_burden_lower": min_uncorrected_lb,
            "min_group_size_Uncorrected_burden_upper": min_uncorrected_ub,
            "min_group_size_Coverage_base": min_cov,
            "min_group_size_Corrected_burden": min_corrected,
            "min_group_size_Corrected_burden_lower": min_corrected_lb,
            "min_group_size_Corrected_burden_upper": min_corrected_ub,
            "group_size_Uncorrected_burden": exact_uncorrected,
            "group_size_Uncorrected_burden_lower": exact_uncorrected_lb,
            "group_size_Uncorrected_burden_upper": exact_uncorrected_ub,
            "group_size_Coverage_base": exact_cov,
            "group_size_Corrected_burden": exact_corrected,
            "group_size_Corrected_burden_lower": exact_corrected_lb,
            "group_size_Corrected_burden_upper": exact_corrected_ub,
        }
    )
    table.to_csv(
        f"{outdir}/{sample}_{mutation_type}_burden_by_group_size.txt",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    panels = (
        (
            axes[0],
            "Burden by minimum group size (min(F1R2,F2R1) >= N)",
            min_uncorrected,
            min_uncorrected_lb,
            min_uncorrected_ub,
            min_corrected,
            min_corrected_lb,
            min_corrected_ub,
        ),
        (
            axes[1],
            "Burden by exact group size (min(F1R2,F2R1) == N)",
            exact_uncorrected,
            exact_uncorrected_lb,
            exact_uncorrected_ub,
            exact_corrected,
            exact_corrected_lb,
            exact_corrected_ub,
        ),
    )
    for (
        ax,
        title,
        uncorrected,
        uncorrected_lb,
        uncorrected_ub,
        corrected,
        corrected_lb,
        corrected_ub,
    ) in panels:
        ax.plot(
            read_numbers,
            uncorrected,
            marker="o",
            color="tab:blue",
            label="Uncorrected burden",
        )
        ax.plot(
            read_numbers, uncorrected_lb, linestyle="--", color="tab:blue", alpha=0.5
        )
        ax.plot(
            read_numbers, uncorrected_ub, linestyle="--", color="tab:blue", alpha=0.5
        )
        ax.plot(
            read_numbers,
            corrected,
            marker="o",
            color="tab:orange",
            label="Corrected burden",
        )
        ax.plot(
            read_numbers, corrected_lb, linestyle="--", color="tab:orange", alpha=0.5
        )
        ax.plot(
            read_numbers, corrected_ub, linestyle="--", color="tab:orange", alpha=0.5
        )
        ax.set_yscale("log")
        ax.set_xlabel("Read number (N)")
        ax.set_ylabel("Mutations per base")
        ax.set_title(title)
        ax.legend()
    fig.suptitle(f"{sample} {mutation_type.upper()} burden by duplex group size")
    fig.tight_layout()
    fig.savefig(f"{outdir}/{sample}_{mutation_type}_burden_by_group_size.pdf")
    plt.close(fig)


def write_genome_extrapolation(f, point, lb, ub, scale):
    """Write the "Mutation number per genome" lines (point, 95% lower/upper)
    as burden * scale.

    point/lb/ub: corrected burden and CI bounds.
    scale: reference base count to extrapolate to, including any unit
        multiplier the burden needs (e.g. indel_locus_multiplier).
    """
    f.write(f"Mutation number per genome\t{point * scale}\n")
    f.write(f"Mutation number per genome 95% lower\t{lb * scale}\n")
    f.write(f"Mutation number per genome 95% upper\t{ub * scale}\n")


def call_flag_set(prefix, sample, flag):
    """True when the call run that wrote prefix had store_true option flag
    set (its _call_params.log)."""
    path = os.path.join(prefix, f"{sample}_call_params.log")
    if not os.path.exists(path):
        return False
    with open(path) as f:
        for line in f:
            key, sep, value = line.strip().partition(":")
            if sep and key == flag:
                return value.strip() == "True"
    return False


def call_skipped_coverage_pass(prefix, sample):
    """True when the call run that wrote prefix used --skipCoveragePass, so
    there is no {sample}_coverage.bed.gz."""
    return call_flag_set(prefix, sample, "skipCoveragePass")


def do_estimate(args):
    log_path = args.prefix + "_estimate_params.log"
    with open(log_path, "w") as log:
        log.write("DupCaller estimate — parameter log\n")
        log.write(f"Version:  {__version__}\n")
        log.write(f"Run time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"Command:  {' '.join(sys.argv)}\n")
        log.write("\n--- All parameters (resolved values) ---\n")
        for key, value in sorted(vars(args).items()):
            log.write(f"  {key}: {value}\n")

    if args.reference:
        bad_h5_files = []
        for suffix, expected_ndim in (
            (".ref.h5", 1),
            (".tn.h5", 1),
            (".hp.h5", 2),
            (".str.h5", 2),
        ):
            ok, msg = check_h5_usable(
                args.reference + suffix, expected_ndim=expected_ndim
            )
            if not ok:
                bad_h5_files.append(msg)
        if bad_h5_files:
            print(
                "ERROR: One or more reference index files are unreadable or incompatible:"
            )
            for msg in bad_h5_files:
                print(f"  - {msg}")
            print(
                "\nPlease re-index the reference genome with the current version of DupCaller:\n"
                f"  DupCaller.py index -f {args.reference} -rt <repeats.tsv>"
            )
            sys.exit(1)

    if not args.reestimatebed:
        ref_trinuc = calculate_ref_trinuc(args)
        prefix = args.prefix
        if len(prefix.split("/")[-1]) == 0:
            sample = prefix.split("/")[-2]
        else:
            sample = prefix.split("/")[-1]
        sbs_dir, indel_dir, dbs_dir = _ensure_type_subdirs(prefix)
        trinuc_by_rf = pd.read_csv(
            sbs_dir + "/" + sample + "_trinuc_by_duplex_group.txt",
            sep="\t",
            index_col=0,
        )
        # Raw-144 DBS opportunity per duplex group; columns are in the same
        # order as trinuc_by_rf.columns (Caller.py writes both from one list).
        dbs144_by_rf = pd.read_csv(
            dbs_dir + "/" + sample + "_dbs_by_duplex_group.txt", sep="\t", index_col=0
        )
        dbs144_by_rf_np = dbs144_by_rf.to_numpy().astype(float)
        ###Estimate SNV burden
        vcf = VCF(sbs_dir + "/" + sample + "_sbs.vcf", "r")
        mut_by_rf = dict()
        # trinuc_by_rf has 192 raw rows "{trinuc}>{alt}" (64 contexts x 3
        # alts), folded to SBS96 later by combine_raw192_to_sbs96.
        label2num_192 = {label: i for i, label in enumerate(trinuc_by_rf.index)}
        trinuc_by_rf_np = trinuc_by_rf.to_numpy().astype(float)
        trinuc_mut_np = np.zeros([96, len(trinuc_by_rf.columns)], dtype=int)
        duplex_no_dict = dict()
        revcomp = _REVCOMP
        for nn, duplex_no in enumerate(trinuc_by_rf.columns):
            duplex_no_dict[duplex_no] = nn

        # (chrom, pos, ref, alt) -> [duplex count, TAC, TDP(, indel cov column)]
        unique_mutations = {}

        if args.dilute:
            vcf_out = VCF(
                sbs_dir + "/" + sample + "_sbs_flt.vcf", "w", header=vcf.header
            )
        stats = parse_stats_file(f"{prefix}/{sample}_stats.txt")
        # "Unmasked Coverage" is already per-locus (call divides by 3).
        # With --doNotOutputMaskedSite there are no masked records, so the
        # Unmasked burden is N/A.
        no_masked = call_flag_set(prefix, sample, "doNotOutputMaskedSite")
        unmasked_cov = None if no_masked else float(stats["Unmasked Coverage"])
        ###Classify PASS DBS events into DBS78 channels. Their two bases are
        ### labeled dbs_member rather than PASS, so they are absent from
        ### _sbs.vcf and SBS/DBS counts are mutually exclusive.
        vcf_dbs = VCF(dbs_dir + "/" + sample + "_dbs.vcf", "r")
        dbs_mut_78 = np.zeros(78, dtype=int)
        # Per-duplex-group counts, columns indexed by duplex_no_dict.
        dbs_mut_by_rf = np.zeros([78, len(trinuc_by_rf.columns)], dtype=int)
        for rec in vcf_dbs.fetch():
            if "PASS" not in rec.filter:
                continue
            label = dbs_raw_event_to_label(rec.ref, rec.alts[0])
            dbs_mut_78[DBS78_LABEL2NUM[label]] += 1
            duplex_no = f"{rec.info['F1R2']}+{rec.info['F2R1']}"
            dbs_mut_by_rf[DBS78_LABEL2NUM[label], duplex_no_dict[duplex_no]] += 1
        vcf_dbs.close()

        ###Calculate masked mutation rate
        # Unmasked burden = all _sbs.vcf records + _sbs_fail.vcf records
        # with filter "masked", SNPM/NOISEM set and DP > 0 (snp/noise-masked
        # sites that cleared LR and had depth extracted). Other fail records
        # have no real depth and must not count.
        if no_masked:
            unmasked_sbs_burden = "N/A"
            unmasked_sbs_burden_lb = unmasked_sbs_burden_ub = "N/A"
        else:
            vcf_fail = VCF(sbs_dir + "/" + sample + "_sbs_fail.vcf", "r")
            unmasked_mut_count = 0
            count_progress = {"start": time.time(), "last": time.time()}
            for i, rec in enumerate(vcf.fetch()):
                log_progress(
                    count_progress,
                    "Counting unmasked SNVs",
                    i,
                    extra=f"{rec.chrom}:{rec.pos}",
                )
                unmasked_mut_count += 1
            for i, rec in enumerate(vcf_fail.fetch()):
                log_progress(
                    count_progress,
                    "Counting unmasked SNVs (masked, rescued by LR)",
                    i,
                    extra=f"{rec.chrom}:{rec.pos}",
                )
                if not (
                    "masked" in rec.filter
                    and (rec.info.get("SNPM") or rec.info.get("NOISEM"))
                    and rec.samples["TUMOR"]["DP"] > 0
                ):
                    continue
                unmasked_mut_count += 1
            vcf_fail.close()
            unmasked_sbs_burden = (
                float(unmasked_mut_count) / float(unmasked_cov)
                if unmasked_cov > 0
                else float("nan")
            )
            unmasked_sbs_burden_lb, unmasked_sbs_burden_ub = fay_feuer_confint(
                unmasked_mut_count, 1.0, unmasked_cov
            )
        snv_progress = {"start": time.time(), "last": time.time()}
        for i, rec in enumerate(vcf.fetch()):
            log_progress(
                snv_progress,
                "Processing SNV mutations",
                i,
                extra=f"{rec.chrom}:{rec.pos}",
            )
            if "PASS" not in rec.filter:
                continue
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = rec.alts[0]
            TAC = rec.samples["TUMOR"]["AC"]
            TDP = rec.samples["TUMOR"]["DP"]
            mutation_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if not unique_mutations.get(mutation_key):
                unique_mutations[mutation_key] = [0, TAC, TDP]
            unique_mutations[mutation_key][0] += 1
            F1R2 = rec.info["F1R2"]
            F2R1 = rec.info["F2R1"]
            if TAC > 1 and args.dilute:
                TDP = rec.samples["TUMOR"]["DP"]
                NAC = rec.samples["NORMAL"]["AC"]
                NDP = rec.samples["NORMAL"]["DP"]
                barnard_p = barnard_exact([[TAC, TDP - TAC], [NAC, NDP - NAC]]).pvalue
                if barnard_p <= 0.05:
                    continue
            if args.dilute:
                vcf_out.write(rec)
            duplex_no = str(F1R2) + "+" + str(F2R1)
            ref = rec.ref
            raw_trinuc = rec.info["TN"]
            raw_alt = rec.alts[0]
            if ref == "C" or ref == "T":
                trinuc = raw_trinuc
                alt = raw_alt
            else:
                trinuc = (
                    revcomp[raw_trinuc[2]]
                    + revcomp[raw_trinuc[1]]
                    + revcomp[raw_trinuc[0]]
                )
                alt = revcomp[raw_alt]
            trinucSbs = trinuc[0] + "[" + trinuc[1] + ">" + alt + "]" + trinuc[2]
            # Mutant loci stay in the coverage denominator: coverage depends
            # only on depth and reference context, not the observed base.
            trinuc_mut_np[TRINUCSBS2NUM_96[trinucSbs], duplex_no_dict[duplex_no]] += 1
        if args.dilute:
            vcf_out.close()

        print("......Estimating mutational burden and SBS96 profile........")
        trinuc_cov_96_by_rf = combine_raw192_to_sbs96(trinuc_by_rf_np, label2num_192)
        (
            mutnum,
            mutnum_corrected,
            correction_ratio,
            corrected_trinuc_num,
            burden,
            burden_lb,
            burden_ub,
            uburden,
            uburden_lb,
            uburden_ub,
            genome_cov,
            cov_by_minread,
            trinuc_opportunity,
            burden_exact,
            burden_exact_lb,
            burden_exact_ub,
            uburden_exact,
            uburden_exact_lb,
            uburden_exact_ub,
            cov_by_exact_group,
        ) = estimate_96(
            trinuc_cov_96_by_rf,
            trinuc_mut_np,
            ref_trinuc,
            trinuc_by_rf.columns,
        )
        # Reference bases considered across args.regions (non-N,
        # ACGT-flanked; no noise BED): each is in exactly one of the 64
        # trinuc classes. Shared genome-size scale for SBS/DBS/indel
        # (ref_indel83.sum() is not a base count: positions overlap columns).
        reference_base_number = genome_cov
        corrected_trinuc_pd = pd.DataFrame(
            np.stack(
                [
                    mutnum,
                    mutnum_corrected,
                    correction_ratio,
                    corrected_trinuc_num,
                    np.repeat(ref_trinuc[:32] + ref_trinuc[32:64], 3),
                    _safe_rate(mutnum_corrected, trinuc_opportunity),
                ],
                axis=1,
            ),
            index=NUM2TRINUCSBS_96,
            columns=[
                "mutation_number_uncorrected",
                "mutation_number_corrected",
                "correction_ratio",
                "mutation_number_genome",
                "trinuc_number_genome",
                "mutations_per_opportunity",
            ],
        )
        corrected_trinuc_pd.to_csv(
            sbs_dir + "/" + sample + "_sbs_96_corrected.txt", sep="\t"
        )
        # One plot page per column (SBS_96_plots_{sample}.pdf). The index
        # must be named "MutationType", otherwise SigProfilerPlotting
        # discards it and uses the first column as the index.
        sigPlt.plotSBS(
            corrected_trinuc_pd[
                ["mutation_number_uncorrected", "mutation_number_corrected"]
            ].rename_axis("MutationType"),
            sbs_dir,
            sample,
            "96",
        )
        write_burden_by_group_size(
            sbs_dir,
            sample,
            "sbs",
            uburden,
            uburden_lb,
            uburden_ub,
            burden,
            burden_lb,
            burden_ub,
            cov_by_minread,
            uburden_exact,
            uburden_exact_lb,
            uburden_exact_ub,
            burden_exact,
            burden_exact_lb,
            burden_exact_ub,
            cov_by_exact_group,
        )

        with open(sbs_dir + "/" + sample + "_sbs_burden.txt", "w") as f:
            f.write(f"Uncorrected burden\t{uburden[0]}\n")
            f.write(f"Uncorrected burden 95% lower\t{uburden_lb[0]}\n")
            f.write(f"Uncorrected burden 95% upper\t{uburden_ub[0]}\n")
            f.write(f"Uncorrected mutation number\t{mutnum.sum()}\n")
            f.write(f"Corrected burden\t{burden[0]}\n")
            f.write(f"Corrected burden 95% lower\t{burden_lb[0]}\n")
            f.write(f"Corrected burden 95% upper\t{burden_ub[0]}\n")
            f.write(f"Corrected mutation number\t{mutnum_corrected.sum()}\n")
            write_genome_extrapolation(
                f, burden[0], burden_lb[0], burden_ub[0], reference_base_number
            )
            # Effective coverage at min_group_size=1.
            f.write(f"Duplex coverage\t{cov_by_minread[0]}\n")
            f.write(f"Unmasked burden\t{unmasked_sbs_burden}\n")
            f.write(f"Unmasked burden 95% lower\t{unmasked_sbs_burden_lb}\n")
            f.write(f"Unmasked burden 95% upper\t{unmasked_sbs_burden_ub}\n")
            f.write(f"Reference base number\t{reference_base_number}\n")

        vcf.close()
        ###Calculate indel burdens
        print("......Estimating Indel Burden.......")
        vcf_indel = VCF(f"{indel_dir}/{sample}_indel.vcf", "r")
        # indel100_by_rf: raw 100-class coverage, folded to indel76 and
        # expanded to ID83 (no microhomology-length resolution), while
        # observed mutations are classified at full ID83 resolution.
        indel100_by_rf = pd.read_csv(
            indel_dir + "/" + sample + "_indel_by_duplex_group.txt",
            sep="\t",
            index_col=0,
        )
        label2num_100 = {label: i for i, label in enumerate(indel100_by_rf.index)}
        indel100_by_rf_np = indel100_by_rf.to_numpy().astype(float)
        label2num_83, labels_indel83 = build_indel83_labels()
        indel83_mut_np = np.zeros([83, len(indel100_by_rf.columns)], dtype=int)
        ref_h5 = h5py.File(args.reference + ".ref.h5", "r")

        indel_count = 0
        indel_progress = {"start": time.time(), "last": time.time()}
        for i, rec in enumerate(vcf_indel.fetch()):
            log_progress(
                indel_progress,
                "Processing indel mutations",
                i,
                extra=f"{rec.chrom}:{rec.pos}",
            )
            mutation_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            indel_count += 1
            TAC = rec.samples["TUMOR"]["AC"]
            TDP = rec.samples["TUMOR"]["DP"]

            F1R2 = rec.info["F1R2"]
            F2R1 = rec.info["F2R1"]
            duplex_no = str(F1R2) + "+" + str(F2R1)
            label, indel_seq, ref_after, indel_len = _indel83_record_label(ref_h5, rec)
            duplex_idx = duplex_no_dict[duplex_no]
            indel83_mut_np[label2num_83[label], duplex_idx] += 1
            if not unique_mutations.get(mutation_key):
                # 4th element: this event's coverage.bed.gz indel category
                # (0-15), used for its duplex depth.
                unique_mutations[mutation_key] = [
                    0,
                    TAC,
                    TDP,
                    indel_coverage_category_index(indel_seq, ref_after, indel_len),
                ]
            unique_mutations[mutation_key][0] += 1

        print("......Estimating indel83 profile and corrected burden........")
        indel76_cov_by_rf = combine_indel100_to_indel76(
            indel100_by_rf_np, label2num_100
        )
        ref_indel100 = calculate_ref_indel100(args)
        ref_indel76 = combine_indel100_to_indel76(ref_indel100, label2num_100)
        indel83_cov_by_rf = expand_indel76_to_indel83(indel76_cov_by_rf)
        ref_indel83 = expand_indel76_to_indel83(ref_indel76)
        indel83_cov_by_rf = override_inshp0_with_next_base_opportunity(
            indel83_cov_by_rf, indel100_by_rf_np
        )
        ref_indel83 = override_inshp0_with_next_base_opportunity(
            ref_indel83, ref_indel100
        )
        (
            indel_mut,
            indel_mut_corrected,
            indel_correction_ratio,
            corrected_indel_num,
            indel_burden_corrected,
            indel_burden_corrected_lb,
            indel_burden_corrected_ub,
            indel_burden_uncorrected,
            indel_burden_uncorrected_lb,
            indel_burden_uncorrected_ub,
            genome_indel_cov,
            indel_cov_by_minread,
            indel_opportunity,
            indel_burden_corrected_exact,
            indel_burden_corrected_exact_lb,
            indel_burden_corrected_exact_ub,
            indel_burden_uncorrected_exact,
            indel_burden_uncorrected_exact_lb,
            indel_burden_uncorrected_exact_ub,
            indel_cov_by_exact_group,
        ) = estimate_indel83(
            indel83_cov_by_rf,
            indel83_mut_np,
            ref_indel83,
            indel100_by_rf.columns,
        )
        indel83_corrected_pd = pd.DataFrame(
            np.stack(
                [
                    indel_mut,
                    indel_mut_corrected,
                    indel_correction_ratio,
                    corrected_indel_num,
                    ref_indel83,
                    _safe_rate(indel_mut_corrected, indel_opportunity),
                ],
                axis=1,
            ),
            index=labels_indel83,
            columns=[
                "mutation_number_uncorrected",
                "mutation_number_corrected",
                "correction_ratio",
                "mutation_number_genome",
                "indel83_number_genome",
                "mutations_per_opportunity",
            ],
        )
        indel83_corrected_pd.to_csv(
            indel_dir + "/" + sample + "_indel_83_corrected.txt", sep="\t"
        )
        # Translate to SigProfilerPlotting's ID83 label syntax
        # (ID_83_plots_{sample}.pdf).
        sigPlt.plotID(
            indel83_corrected_pd[
                ["mutation_number_uncorrected", "mutation_number_corrected"]
            ]
            .rename(index=INDEL83_TO_SIGPROFILER_LABELS)
            .rename_axis("MutationType"),
            indel_dir,
            sample,
            "83",
        )
        # SBS/DBS coverage has a fixed per-locus multiplicity (3 alts, 9 alt
        # pairs) already divided out, but indel coverage sums overlapping,
        # locus-varying opportunity columns, so estimate_indel83's burdens
        # are per opportunity column. indel_locus_multiplier is the
        # reference genome's columns-per-base ratio; multiplying by it
        # (assuming covered loci share the genome-average multiplicity)
        # gives per-base burdens. Coverage is divided by it.
        indel_locus_multiplier = genome_indel_cov / reference_base_number
        write_burden_by_group_size(
            indel_dir,
            sample,
            "indel",
            indel_burden_uncorrected * indel_locus_multiplier,
            indel_burden_uncorrected_lb * indel_locus_multiplier,
            indel_burden_uncorrected_ub * indel_locus_multiplier,
            indel_burden_corrected * indel_locus_multiplier,
            indel_burden_corrected_lb * indel_locus_multiplier,
            indel_burden_corrected_ub * indel_locus_multiplier,
            indel_cov_by_minread / indel_locus_multiplier,
            indel_burden_uncorrected_exact * indel_locus_multiplier,
            indel_burden_uncorrected_exact_lb * indel_locus_multiplier,
            indel_burden_uncorrected_exact_ub * indel_locus_multiplier,
            indel_burden_corrected_exact * indel_locus_multiplier,
            indel_burden_corrected_exact_lb * indel_locus_multiplier,
            indel_burden_corrected_exact_ub * indel_locus_multiplier,
            indel_cov_by_exact_group / indel_locus_multiplier,
        )
        # Same layout as _sbs_burden.txt, minus the "Unmasked burden" lines
        # (no per-base unmasked indel coverage to divide by). Uncorrected and
        # corrected burden share estimate_indel83's ID83 coverage denominator.
        with open(indel_dir + "/" + sample + "_indel_burden.txt", "w") as f:
            f.write(
                f"Uncorrected burden\t"
                f"{indel_burden_uncorrected[0] * indel_locus_multiplier}\n"
            )
            f.write(
                f"Uncorrected burden 95% lower\t"
                f"{indel_burden_uncorrected_lb[0] * indel_locus_multiplier}\n"
            )
            f.write(
                f"Uncorrected burden 95% upper\t"
                f"{indel_burden_uncorrected_ub[0] * indel_locus_multiplier}\n"
            )
            f.write(f"Uncorrected mutation number\t{indel_count}\n")
            f.write(
                f"Corrected burden\t"
                f"{indel_burden_corrected[0] * indel_locus_multiplier}\n"
            )
            f.write(
                f"Corrected burden 95% lower\t"
                f"{indel_burden_corrected_lb[0] * indel_locus_multiplier}\n"
            )
            f.write(
                f"Corrected burden 95% upper\t"
                f"{indel_burden_corrected_ub[0] * indel_locus_multiplier}\n"
            )
            f.write(f"Corrected mutation number\t{indel_mut_corrected.sum()}\n")
            write_genome_extrapolation(
                f,
                indel_burden_corrected[0],
                indel_burden_corrected_lb[0],
                indel_burden_corrected_ub[0],
                indel_locus_multiplier * reference_base_number,
            )
            # Effective coverage at min_group_size=1, in per-base units.
            f.write(
                f"Duplex coverage\t"
                f"{indel_cov_by_minread[0] / indel_locus_multiplier}\n"
            )
            f.write(f"Reference base number\t{reference_base_number}\n")

        ###Calculate DBS burden
        print("......Estimating DBS78 Burden.......")
        dbs_opp144 = dbs144_by_rf_np.sum(axis=1)
        dbs_opp_78 = combine_raw_dbs_to_dbs78(dbs_opp144, DBS_RAW144_LABEL2NUM)
        ref_dbs_16 = calculate_ref_dbs(args)
        ref_dbs_78 = _dbs78_ref_weighted(ref_dbs_16)
        (
            dbs_mut,
            dbs_mut_corrected,
            dbs_correction_ratio,
            dbs_burden_corrected,
            dbs_burden_corrected_lb,
            dbs_burden_corrected_ub,
            dbs_burden_uncorrected,
            dbs_burden_uncorrected_lb,
            dbs_burden_uncorrected_ub,
            dbs_cov,
        ) = estimate_dbs78(dbs_mut_78, dbs_opp_78, ref_dbs_78)

        dbs78_corrected_pd = pd.DataFrame(
            np.stack(
                [
                    dbs_mut,
                    dbs_mut_corrected,
                    dbs_correction_ratio,
                    ref_dbs_78,
                    _safe_rate(dbs_mut_corrected, dbs_opp_78),
                ],
                axis=1,
            ),
            index=DBS78_LABELS,
            columns=[
                "mutation_number_uncorrected",
                "mutation_number_corrected",
                "correction_ratio",
                "dbs78_number_genome",
                "mutations_per_opportunity",
            ],
        )
        dbs78_corrected_pd.to_csv(
            dbs_dir + "/" + sample + "_dbs_78_corrected.txt", sep="\t"
        )
        # DBS78 labels already match SigProfilerPlotting's
        # (DBS_78_plots_{sample}.pdf).
        sigPlt.plotDBS(
            dbs78_corrected_pd[
                ["mutation_number_uncorrected", "mutation_number_corrected"]
            ].rename_axis("MutationType"),
            dbs_dir,
            sample,
            "78",
        )

        (
            dbs_burden_uncorrected_min,
            dbs_burden_uncorrected_min_lb,
            dbs_burden_uncorrected_min_ub,
            dbs_burden_corrected_min,
            dbs_burden_corrected_min_lb,
            dbs_burden_corrected_min_ub,
            dbs_cov_by_minread,
            dbs_burden_uncorrected_exact,
            dbs_burden_uncorrected_exact_lb,
            dbs_burden_uncorrected_exact_ub,
            dbs_burden_corrected_exact,
            dbs_burden_corrected_exact_lb,
            dbs_burden_corrected_exact_ub,
            dbs_cov_by_exact_group,
        ) = estimate_dbs78_by_group(
            dbs_mut_by_rf, dbs144_by_rf_np, ref_dbs_78, dbs144_by_rf.columns
        )
        write_burden_by_group_size(
            dbs_dir,
            sample,
            "dbs",
            dbs_burden_uncorrected_min,
            dbs_burden_uncorrected_min_lb,
            dbs_burden_uncorrected_min_ub,
            dbs_burden_corrected_min,
            dbs_burden_corrected_min_lb,
            dbs_burden_corrected_min_ub,
            dbs_cov_by_minread,
            dbs_burden_uncorrected_exact,
            dbs_burden_uncorrected_exact_lb,
            dbs_burden_uncorrected_exact_ub,
            dbs_burden_corrected_exact,
            dbs_burden_corrected_exact_lb,
            dbs_burden_corrected_exact_ub,
            dbs_cov_by_exact_group,
        )

        # Same layout as _sbs_burden.txt, minus "Unmasked burden" (DBS calls
        # have no masked tier).
        with open(dbs_dir + "/" + sample + "_dbs_burden.txt", "w") as f:
            f.write(f"Uncorrected burden\t{dbs_burden_uncorrected}\n")
            f.write(f"Uncorrected burden 95% lower\t{dbs_burden_uncorrected_lb}\n")
            f.write(f"Uncorrected burden 95% upper\t{dbs_burden_uncorrected_ub}\n")
            f.write(f"Uncorrected mutation number\t{int(dbs_mut.sum())}\n")
            f.write(f"Corrected burden\t{dbs_burden_corrected}\n")
            f.write(f"Corrected burden 95% lower\t{dbs_burden_corrected_lb}\n")
            f.write(f"Corrected burden 95% upper\t{dbs_burden_corrected_ub}\n")
            f.write(f"Corrected mutation number\t{dbs_mut_corrected.sum()}\n")
            write_genome_extrapolation(
                f,
                dbs_burden_corrected,
                dbs_burden_corrected_lb,
                dbs_burden_corrected_ub,
                reference_base_number,
            )
            # Effective coverage at min_group_size=1.
            f.write(f"Duplex coverage\t{dbs_cov_by_minread[0]}\n")
            f.write(f"Reference base number\t{reference_base_number}\n")

        # Append per-base duplex coverage (same values as each burden file's
        # "Duplex coverage") to the stats.txt written by call.
        with open(f"{prefix}/{sample}_stats.txt", "a") as f:
            f.write(f"SBS Base Coverage\t{cov_by_minread[0]}\n")
            f.write(
                f"Indel Base Coverage\t"
                f"{indel_cov_by_minread[0] / indel_locus_multiplier}\n"
            )
            f.write(f"DBS Base Coverage\t{dbs_cov_by_minread[0]}\n")

        # A call run with --skipCoveragePass writes no coverage.bed.gz, which
        # the duplex allele counts and per-gene coverage read per locus.
        no_coverage_file = call_skipped_coverage_pass(prefix, sample)
        if no_coverage_file:
            print("Skipping duplex allele count due to no coverage file.")
        else:
            print("......Processing unique mutations for duplex allele counts.......")
            coverage_file = f"{prefix}/{sample}_coverage.bed.gz"
            allele_counts_data = []

            tbx = TabixFile(coverage_file)

            n_unique_muts = len(unique_mutations)
            allele_counts_progress = {"start": time.time(), "last": time.time()}
            for i, ((chrom, pos, ref, alt), counts) in enumerate(
                unique_mutations.items()
            ):
                log_progress(
                    allele_counts_progress,
                    "Extracting duplex allele counts",
                    i,
                    n_unique_muts,
                    extra=f"{chrom}:{pos}",
                )
                duplex_depth = 0

                try:
                    # Coverage bed format: chrom start end A_cov T_cov C_cov G_cov
                    # del_unit ins_unit del_2 del_3 del_4 del_5plus
                    # ins_A ins_T ins_C ins_G ins_len2 ins_len3 ins_len4 ins_len5plus
                    # del_hp ins_hp
                    # Columns 3-6: per-alt-base L-weighted SNV coverage;
                    # columns 7-22: power-weighted indel coverage per
                    # INDEL_COVERAGE_CATEGORY_LABELS category. Depth is kept
                    # as a float (3dp): int truncation would zero out, and so
                    # drop, real sub-1.0 single-family sites.
                    alt_col = {"A": 3, "T": 4, "C": 5, "G": 6}
                    is_snv = len(ref) == 1 and len(alt) == 1 and alt in alt_col
                    # coverage.bed.gz is keyed at the SNV position, or for an
                    # indel at the base after its anchor (0-based index = POS).
                    query_pos = pos - 1 if is_snv else pos
                    for row in tbx.fetch(chrom, query_pos, query_pos + 1):
                        parts = row.split("\t")
                        if len(parts) >= 23:
                            if is_snv:
                                duplex_depth = round(float(parts[alt_col[alt]]), 3)
                            else:
                                # Indel: only this event's own category column
                                # (counts[3]).
                                duplex_depth = round(float(parts[7 + counts[3]]), 3)
                            break
                    gene_name = "."
                    if args.genebed:
                        for g in TabixFile(args.genebed).fetch(chrom, pos - 1, pos):
                            gene_name = g.split("\t")[3].split("_")[0]
                            break

                except Exception as e:
                    print(f"Warning: Could not query coverage for {chrom}:{pos}: {e}")
                    duplex_depth = 0
                if duplex_depth == 0:
                    print(
                        f"duplex depth for location {chrom}:{pos} is 0, skipping the mutation..."
                    )
                    continue
                allele_counts_data.append(
                    {
                        "chromosome": chrom,
                        "position_start": pos,
                        "ref": ref,
                        "alt": alt,
                        "count": counts[0],
                        "duplex_depth": duplex_depth,
                        "bam_alt_count": counts[1],
                        "bam_depth": counts[2],
                        "duplex_vaf": float(counts[0]) / float(duplex_depth),
                        "bam_vaf": float(counts[1]) / float(counts[2]),
                        "gene": gene_name,
                    }
                )
            tbx.close()

            allele_counts_df = pd.DataFrame(allele_counts_data)
            allele_counts_file = (
                args.prefix + "/" + sample + "_duplex_allele_counts.txt"
            )
            allele_counts_df.to_csv(allele_counts_file, sep="\t", index=False)
            print(f"Duplex allele counts written to: {allele_counts_file}")
        if args.genebed and no_coverage_file:
            print("Skipping per-gene duplex coverage due to no coverage file.")
        elif args.genebed:
            print(f"Calculating per-gene duplex coverage")
            gene_sbs_cov = dict()
            gene_indel_cov = dict()
            cds_len = dict()
            coverage_tbx = TabixFile(f"{prefix}/{sample}_coverage.bed.gz")
            for rec in TabixFile(args.genebed).fetch():
                (chrom, start, end, gene_exon) = rec.split("\t")[0:4]
                start = int(start)
                end = int(end)
                gene = gene_exon.split("_")[0]
                exon_sbs_cov = 0.0
                exon_indel_cov = 0.0
                for loc in coverage_tbx.fetch(chrom, start, end):
                    parts = loc.split("\t")
                    # Columns 3-6: per-alt-base SNV coverage (ref-base column
                    # is 0), so sum / 3 below is per-base depth.
                    exon_sbs_cov += sum(float(x) for x in parts[3:7])
                    # Columns 7-22: indel coverage, converted to per-base
                    # depth below via indel_locus_multiplier.
                    exon_indel_cov += sum(float(x) for x in parts[7:23])
                if gene in cds_len:
                    gene_sbs_cov[gene] += exon_sbs_cov
                    gene_indel_cov[gene] += exon_indel_cov
                    cds_len[gene] += end - start
                else:
                    gene_sbs_cov[gene] = exon_sbs_cov
                    gene_indel_cov[gene] = exon_indel_cov
                    cds_len[gene] = end - start
            coverage_tbx.close()
            with open(args.prefix + "/" + sample + "_gene_coverage.txt", "w") as f:
                f.write("gene\tsbs_duplex_depth\tindel_duplex_depth\n")
                for gene in cds_len.keys():
                    sbs_depth = gene_sbs_cov[gene] / 3 / cds_len[gene]
                    indel_depth = (
                        gene_indel_cov[gene] / indel_locus_multiplier / cds_len[gene]
                    )
                    f.write(f"{gene}\t{sbs_depth}\t{indel_depth}\n")

    else:
        print("Re-estimating mutational burden from provided bed file")
        ref_trinuc = calculate_ref_trinuc(args)
        prefix = args.prefix
        if len(prefix.split("/")[-1]) == 0:
            sample = prefix.split("/")[-2]
        else:
            sample = prefix.split("/")[-1]
        if call_skipped_coverage_pass(prefix, sample):
            print(
                "Error: -rb/--reestimatebed needs the per-locus coverage.bed.gz, "
                "which call --skipCoveragePass does not write."
            )
            sys.exit(1)
        sbs_dir, indel_dir, dbs_dir = _ensure_type_subdirs(prefix)
        tn_int = h5py.File(args.reference + ".tn.h5", "r")
        # trinuc_sbs_cov_64[t, b]: L-weighted coverage for raw trinuc
        # context t (0-63) and alt base b; folded to SBS96 at the end.
        trinuc2num_64, num2trinuc_64 = build_trinuc64_order()
        label2num_192, labels_192 = build_trinuc192_labels(num2trinuc_64)
        base2num_re = {"A": 0, "T": 1, "C": 2, "G": 3}
        row_idx_192 = np.array(
            [trinuc2num_64[label.split(">")[0]] for label in labels_192]
        )
        col_idx_192 = np.array(
            [base2num_re[label.split(">")[1]] for label in labels_192]
        )
        trinuc_sbs_cov_64 = np.zeros((64, 4))
        trinucSbs_count = np.zeros(96)
        # Indel coverage redistributed from coverage.bed.gz's 16 columns to
        # indel100 (redistribute_indel_cov16_to_100), then folded to ID83 so
        # the same correction scheme as estimate_indel83 applies.
        label2num_100_re, _ = build_indel100_labels()
        label2num_83, labels_indel83 = build_indel83_labels()
        indel_cov100_re = np.zeros(100)
        indel_mut83_re = np.zeros(83)
        revcomp = _REVCOMP
        vcf = VCF(sbs_dir + "/" + sample + "_sbs.vcf", "r")
        vcf_indel = VCF(indel_dir + "/" + sample + "_indel.vcf", "r")
        ref_h5 = h5py.File(args.reference + ".ref.h5", "r")
        hp_h5 = h5py.File(args.reference + ".hp.h5", "r")
        str_h5 = h5py.File(args.reference + ".str.h5", "r")
        indel_count = 0
        for interval in TabixFile(args.reestimatebed, parser=pysam.asBed()).fetch():
            # Annotation slices fetched once per interval (h5 reads per row
            # are slow), indexed by pos - interval.start.
            hp_run_slice = hp_h5[interval.contig][0, interval.start : interval.end]
            str_slice = str_h5[interval.contig][:, interval.start : interval.end]
            ref_slice = ref_h5[interval.contig][interval.start : interval.end]
            for loc in TabixFile(f"{prefix}/{sample}_coverage.bed.gz").fetch(
                interval.contig, interval.start, interval.end
            ):
                # Coverage bed format: chrom start end A_cov T_cov C_cov G_cov
                # del_unit ins_unit del_2 del_3 del_4 del_5plus
                # ins_A ins_T ins_C ins_G ins_len2 ins_len3 ins_len4 ins_len5plus
                # del_hp ins_hp
                parts = loc.split("\t")
                a_cov = float(parts[3])
                t_cov = float(parts[4])
                c_cov = float(parts[5])
                g_cov = float(parts[6])
                pos = int(parts[1])
                local_idx = pos - interval.start
                redistribute_indel_cov16_to_100(
                    [float(v) for v in parts[7:23]],
                    int(hp_run_slice[local_idx]),
                    int(str_slice[0, local_idx]),
                    int(str_slice[1, local_idx]),
                    int(ref_slice[local_idx]),
                    label2num_100_re,
                    indel_cov100_re,
                )
                trinuc_idx = int(tn_int[parts[0]][int(parts[1])])
                if trinuc_idx < 64:
                    trinuc_sbs_cov_64[trinuc_idx, 0] += a_cov
                    trinuc_sbs_cov_64[trinuc_idx, 1] += t_cov
                    trinuc_sbs_cov_64[trinuc_idx, 2] += c_cov
                    trinuc_sbs_cov_64[trinuc_idx, 3] += g_cov
            for rec in vcf_indel.fetch():
                if "PASS" not in rec.filter:
                    continue
                if rec.chrom != interval.contig:
                    continue
                if rec.pos <= interval.start or rec.pos > interval.end:
                    continue
                indel_count += 1
                label = _indel83_record_label(ref_h5, rec)[0]
                indel_mut83_re[label2num_83[label]] += 1
            for rec in vcf.fetch():
                if "PASS" not in rec.filter:
                    continue
                if rec.chrom != interval.contig:
                    continue
                if rec.pos <= interval.start or rec.pos > interval.end:
                    continue
                chrom = rec.chrom
                pos = rec.pos
                ref = rec.ref
                alt = rec.alts[0]
                F1R2 = rec.info["F1R2"]
                F2R1 = rec.info["F2R1"]
                TAC = rec.samples["TUMOR"]["AC"]
                if TAC > 1 and args.dilute:
                    TDP = rec.samples["TUMOR"]["DP"]
                    NAC = rec.samples["NORMAL"]["AC"]
                    NDP = rec.samples["NORMAL"]["DP"]
                    barnard_p = barnard_exact(
                        [[TAC, TDP - TAC], [NAC, NDP - NAC]]
                    ).pvalue
                    if barnard_p <= 0.05:
                        continue
                if ref == "C" or ref == "T":
                    trinuc = rec.info["TN"]
                    alt = rec.alts[0]
                else:
                    trinuc_revcomp = rec.info["TN"]
                    trinuc = (
                        revcomp[trinuc_revcomp[2]]
                        + revcomp[trinuc_revcomp[1]]
                        + revcomp[trinuc_revcomp[0]]
                    )
                    alt = revcomp[rec.alts[0]]
                trinucSbs = trinuc[0] + "[" + trinuc[1] + ">" + alt + "]" + trinuc[2]
                trinucSbs_count[TRINUCSBS2NUM_96[trinucSbs]] += 1

        # Select the 192 raw cells, then fold to SBS96.
        trinuc_cov_192 = trinuc_sbs_cov_64[row_idx_192, col_idx_192]
        trinuc_cov_96 = combine_raw192_to_sbs96(trinuc_cov_192, label2num_192)
        ref_trinuc_96 = np.repeat(ref_trinuc[:32] + ref_trinuc[32:64], 3)  # (96,)
        trinucSbs_count = trinucSbs_count.astype(float)
        ref_trinuc_sum = ref_trinuc.sum()
        genome_cov = ref_trinuc_sum

        trinuc_mut = trinucSbs_count
        mutnum = trinuc_mut.sum()
        cov = trinuc_cov_96.sum() / 3  # per-locus-equivalent coverage
        burden_uncorrected = mutnum / cov
        burden_uncorrected_lb, burden_uncorrected_ub = fay_feuer_confint(
            mutnum, 1.0, cov
        )
        trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov_96)
        correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov_96)
        mutnum_corrected = correction_ratio * trinuc_mut
        burden_corrected = mutnum_corrected.sum() / cov
        burden_corrected_lb, burden_corrected_ub = fay_feuer_confint(
            trinuc_mut, correction_ratio, cov
        )
        hap_trinuc = trinuc_rate * ref_trinuc_96

        with open(sbs_dir + "/" + sample + "_sbs_burden_re_estimate.txt", "w") as f:
            f.write(f"Uncorrected burden\t{burden_uncorrected}\n")
            f.write(f"Uncorrected burden 95% lower\t{burden_uncorrected_lb}\n")
            f.write(f"Uncorrected burden 95% upper\t{burden_uncorrected_ub}\n")
            f.write(f"Uncorrected mutation number\t{int(mutnum)}\n")
            f.write(f"Corrected burden\t{burden_corrected}\n")
            f.write(f"Corrected burden 95% lower\t{burden_corrected_lb}\n")
            f.write(f"Corrected burden 95% upper\t{burden_corrected_ub}\n")
            f.write(f"Corrected mutation number\t{mutnum_corrected.sum()}\n")
            write_genome_extrapolation(
                f,
                burden_corrected,
                burden_corrected_lb,
                burden_corrected_ub,
                genome_cov,
            )
            # Per-locus duplex coverage (the burden denominator).
            f.write(f"Duplex coverage\t{cov}\n")
            f.write(f"Reference base number\t{genome_cov}\n")
        # indel100 -> indel76 -> ID83 (+ rep0 fill), as in the main path.
        ref_indel100_re = calculate_ref_indel100(args)
        ref_indel76_re = combine_indel100_to_indel76(ref_indel100_re, label2num_100_re)
        ref_indel83_re = expand_indel76_to_indel83(ref_indel76_re)
        ref_indel83_re = override_inshp0_with_next_base_opportunity(
            ref_indel83_re, ref_indel100_re
        )
        indel_cov76_re = combine_indel100_to_indel76(indel_cov100_re, label2num_100_re)
        indel_cov83_re = expand_indel76_to_indel83(indel_cov76_re)
        indel_cov83_re = override_inshp0_with_next_base_opportunity(
            indel_cov83_re, indel_cov100_re
        )
        # Coverage is per opportunity channel; rescale to per-base units
        # (see indel_locus_multiplier in the main path).
        indel_locus_multiplier_re = ref_indel83_re.sum() / ref_trinuc.sum()
        indel_cov_total = indel_cov83_re.sum()
        indel_burden_re_uncorrected = (
            indel_count / float(indel_cov_total) * indel_locus_multiplier_re
            if indel_cov_total > 0
            else float("nan")
        )
        (
            indel_burden_re_uncorrected_lb,
            indel_burden_re_uncorrected_ub,
        ) = fay_feuer_confint(indel_count, 1.0, indel_cov_total)
        indel_burden_re_uncorrected_lb *= indel_locus_multiplier_re
        indel_burden_re_uncorrected_ub *= indel_locus_multiplier_re

        ref_frac_re = ref_indel83_re / ref_indel83_re.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            obs_frac_re = indel_cov83_re / indel_cov83_re.sum()
        indel_correction_ratio83_re = _grouped_indel83_correction_ratio(
            ref_frac_re, obs_frac_re
        )
        indel_mutnum_corrected = (indel_correction_ratio83_re * indel_mut83_re).sum()
        indel_burden_re_corrected = (
            indel_mutnum_corrected / float(indel_cov_total) * indel_locus_multiplier_re
            if indel_cov_total > 0
            else float("nan")
        )
        (
            indel_burden_re_corrected_lb,
            indel_burden_re_corrected_ub,
        ) = fay_feuer_confint(
            indel_mut83_re, indel_correction_ratio83_re, indel_cov_total
        )
        indel_burden_re_corrected_lb *= indel_locus_multiplier_re
        indel_burden_re_corrected_ub *= indel_locus_multiplier_re
        with open(indel_dir + "/" + sample + "_indel_burden_re_estimate.txt", "w") as f:
            f.write(f"Uncorrected indel burden\t{indel_burden_re_uncorrected}\n")
            f.write(
                f"Uncorrected indel burden 95% lower\t{indel_burden_re_uncorrected_lb}\n"
            )
            f.write(
                f"Uncorrected indel burden 95% upper\t{indel_burden_re_uncorrected_ub}\n"
            )
            f.write(f"Uncorrected indel number\t{indel_count}\n")
            f.write(f"Corrected indel burden\t{indel_burden_re_corrected}\n")
            f.write(
                f"Corrected indel burden 95% lower\t{indel_burden_re_corrected_lb}\n"
            )
            f.write(
                f"Corrected indel burden 95% upper\t{indel_burden_re_corrected_ub}\n"
            )
            f.write(f"Corrected indel number\t{indel_mutnum_corrected}\n")
            # Burdens are already per-base (indel_locus_multiplier_re).
            write_genome_extrapolation(
                f,
                indel_burden_re_corrected,
                indel_burden_re_corrected_lb,
                indel_burden_re_corrected_ub,
                genome_cov,
            )
            f.write(f"Indel coverage\t{indel_cov_total}\n")
        corrected_trinuc_pd = pd.DataFrame(
            np.stack(
                [
                    trinuc_mut,
                    mutnum_corrected,
                    correction_ratio,
                    hap_trinuc,
                    np.repeat(ref_trinuc[:32] + ref_trinuc[32:64], 3),
                ],
                axis=1,
            ),
            index=NUM2TRINUCSBS_96,
            columns=[
                "mutation_number_uncorrected",
                "mutation_number_corrected",
                "correction_ratio",
                "mutation_number_genome",
                "trinuc_number_genome",
            ],
        )
        corrected_trinuc_pd.to_csv(
            sbs_dir + "/" + sample + "_sbs_96_corrected_re_estimate.txt", sep="\t"
        )
        sigPlt.plotSBS(
            corrected_trinuc_pd[
                ["mutation_number_uncorrected", "mutation_number_corrected"]
            ].rename_axis("MutationType"),
            sbs_dir,
            sample + "_re_estimate",
            "96",
        )
