#!/usr/bin/env python3
import sys
import time
import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sigProfilerPlotting as sigPlt
from pysam import VariantFile as VCF, TabixFile
from scipy.stats import chi2, barnard_exact
import pysam
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
    """Reference-genome trinucleotide composition at full 64-context
    resolution (un-folded by reverse complement). Reverse-complement pairs
    are only combined at estimation time, inside estimate_96.
    """
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


# base2num convention shared with Index.py/funcs/call.py's nums2str
# (A=0,T=1,C=2,G=3); anything outside 0-3 (N bases, contig edges) decodes
# to "N", which simply never matches in classify_indel_channel's sequence
# comparisons.
_NUM2BASE = "ATCG"


def _decode_ref_seq(int_arr):
    """Decode a ref.h5 base2num-encoded integer array into an uppercase
    base string, for classify_indel_channel (funcs/misc.py), which needs
    actual sequence rather than pre-computed repeat annotation."""
    return "".join(_NUM2BASE[b] if 0 <= b <= 3 else "N" for b in int_arr)


# How many reference bases to fetch after an indel when resolving its
# ID83 channel: classify_indel_record's MH-vs-repeat-unit check only needs
# up to abs(indel_len) bases (typical deletions are a handful of bp), and
# repeat_count itself now comes straight from hp.h5/str.h5 (classify_
# indel_channel's `anno`), not from scanning ref_after -- this window only
# needs to cover the deletion/insertion length itself.
INDEL_CONTEXT_WINDOW = 50


def calculate_ref_indel100(args):
    """Genome-wide indel100 "opportunity" composition (unweighted position
    counts, one contribution per eligible position per category — mirrors
    indel100_reference_bucket_indices' overlap conventions), analogous to
    calculate_ref_trinuc for SBS. Read directly from the
    .hp.h5/.str.h5/.ref.h5 index files rather than the coverage bed, since
    this is genome composition, not observed coverage.

    Reads hp.h5 and str.h5 raw (not through load_repeat_context's merged,
    STR-priority view) so homopolymer and STR opportunity stay
    independent — see indel100_reference_bucket_indices' docstring for
    why a position can validly contribute to both.
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
        # ref.h5 is already uint8 on disk (base2num-encoded); no need to
        # widen it just to compute a shifted copy.
        ref_base_arr = ref_h5[chrom][:]
        # Reference base immediately following each position — the last
        # position in the chromosome has no valid "next" base, hence -1
        # (needs a signed dtype for that sentinel, unlike ref_base_arr).
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


def poisson_confint(k, cov, alpha=0.05):
    if cov == 0:
        return float("nan"), float("nan")
    low = chi2.ppf(alpha / 2, 2 * k) / 2
    high = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    if k == 0:
        low = 0
    return low / cov, high / cov


def _dinuc_num(dinuc):
    """Encode a 2-base string as an int 0-15: 4*base2num[dinuc[0]] +
    base2num[dinuc[1]] (A=0,T=1,C=2,G=3, same convention as ref.h5)."""
    b2n = {"A": 0, "T": 1, "C": 2, "G": 3}
    return 4 * b2n[dinuc[0]] + b2n[dinuc[1]]


def _dinuc_revcomp(dinuc):
    return _REVCOMP[dinuc[1]] + _REVCOMP[dinuc[0]]


# The 10 canonical DBS78 reference dinucleotides. Not derivable from a
# simple rule (e.g. NOT "alphabetically smaller of {ref, RC(ref)}" — the
# AA/TT pair's representative is "TT", not "AA") — this is a fixed,
# external convention (COSMIC/SigProfiler), extracted from
# SigProfilerPlotting's plotDBS (sigProfilerPlotting.py:10444-10455,
# LabelEncoder over these 10 strings, already alphabetically sorted) so
# our labels line up with the same palette plot_dbs78 uses.
DBS_CANONICAL_REFS = ["AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"]


def _build_dbs_raw_to_canonical_ref():
    """All 16 raw dinucleotides -> their canonical-ref representative.
    The 6 non-canonical raw refs are each some canonical ref's reverse
    complement; derived here rather than hardcoded so it can't drift out
    of sync with DBS_CANONICAL_REFS."""
    mapping = {}
    for ref in DBS_CANONICAL_REFS:
        mapping[ref] = ref
        mapping[_dinuc_revcomp(ref)] = ref
    assert len(mapping) == 16
    return mapping


DBS_RAW_TO_CANONICAL_REF = _build_dbs_raw_to_canonical_ref()

# For each self-complementary ref (AT, TA, CG, GC — see build_dbs78_labels'
# docstring), which 6 of its 9 raw alts are the canonical-direction
# representative. Not derivable from a simple rule (e.g. not "the
# alphabetically smaller of {alt, RC(alt)}" — for AT, "CA" is kept over its
# partner "TG" despite "CA" > "TG"... sorted the other way for other pairs
# too), so this is taken verbatim from SigProfilerPlotting's plotDBS
# (sigProfilerPlotting.py:10330-10409, the `dinucs` literal), which is
# itself just the fixed COSMIC DBS78 reference catalog — matching it
# exactly (rather than picking an internally-consistent but arbitrary
# direction) is what makes our DBS78 channel labels interoperable with
# COSMIC/SigProfiler DBS signature references.
_DBS78_SELFCOMP_CANONICAL_ALTS = {
    "AT": ["CA", "CC", "CG", "GA", "GC", "TA"],
    "TA": ["AT", "CG", "CT", "GC", "GG", "GT"],
    "CG": ["AT", "GC", "GT", "TA", "TC", "TT"],
    "GC": ["AA", "AG", "AT", "CA", "CG", "TA"],
}


def build_dbs78_labels():
    """78 canonical DBS labels "REF>ALT", deterministic order, plus the
    raw-alt -> canonical-alt map needed to fold observed/raw events onto
    them.

    Why 78 and not 10*9=90: dinucleotide contexts (unlike trinucleotide
    ones) can be their OWN reverse complement — AT, TA, CG, and GC each
    satisfy RC(ref)==ref (e.g. RC("AT") = comp("T")+comp("A") = "A"+"T" =
    "AT"). This is impossible for odd-length contexts like trinucleotides
    (the middle base would have to equal its own complement, and no base
    does), so it doesn't show up anywhere in the SBS96/ID83 folding logic
    above — DBS is the first place this matters. For those 4
    "self-complementary" refs, the ref doesn't change identity under
    reverse-complementing, but its 9 alts still pair up among
    THEMSELVES: e.g. for ref="AT", observing alt="CA" and observing
    alt="TG" are the same physical event viewed from opposite strands
    (RC("AT","CA") = ("AT","TG")), so they must be folded into one
    channel. Working this out per self-complementary ref leaves 6 distinct
    alt channels (3 self-complementary alts that don't pair with anything
    + 3 merged pairs), not 9 — 4 refs * 6 + 6 refs * 9 = 24 + 54 = 78,
    matching the standard COSMIC DBS78 scheme. Which member of each merged
    pair is the "canonical" direction is an external convention with no
    derivable rule, taken from _DBS78_SELFCOMP_CANONICAL_ALTS.
    """
    label2num = {}
    labels = []
    # (canonical_ref, raw_alt) -> canonical_alt, for every raw alt of
    # every canonical ref — used by dbs_raw_event_to_label below.
    alt_canonical = {}

    for ref in DBS_CANONICAL_REFS:
        alts = _dbs_alt_choices(ref)
        if _dinuc_revcomp(ref) == ref:
            # Self-complementary ref: each alt is either already one of
            # the 6 official canonical alts, or its RC is (exactly one of
            # the two always holds, since the 3 self-inverse alts are
            # trivially canonical and the remaining 6 form 3 RC-pairs with
            # exactly one representative each in the official set).
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
    """Fold one observed (raw_ref, raw_alt) dinucleotide substitution —
    as actually read off the reference/alt strand, no canonicalization
    applied yet — onto its DBS78 channel label. Mirrors how SBS defers
    reverse-complement folding to estimation time (combine_raw192_to_sbs96)
    rather than canonicalizing at call time."""
    canonical_ref = DBS_RAW_TO_CANONICAL_REF[raw_ref]
    eff_alt = raw_alt if canonical_ref == raw_ref else _dinuc_revcomp(raw_alt)
    canonical_alt = _DBS78_ALT_CANONICAL[(canonical_ref, eff_alt)]
    return f"{canonical_ref}>{canonical_alt}"


DBS_RAW144_LABEL2NUM, DBS_RAW144_LABELS = build_dbs_raw144_labels()


def combine_raw_dbs_to_dbs78(values_144, label2num_144):
    """Fold a raw-144-space DBS array down to the 78 canonical DBS78
    classes, by summing each canonical class with whatever raw (ref,alt)
    combination(s) map onto it via dbs_raw_event_to_label — unlike SBS96's
    uniform 2-raw-classes-per-canonical-class fold, this varies (2 for the
    54 channels belonging to non-self-complementary refs, but the 24
    channels belonging to the 4 self-complementary refs can fold from
    either a same-ref alt-pair or a different-ref mirror pair depending on
    the raw event), so this looks up the mapping per raw label rather than
    assuming a fixed pairing.

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
    """Genome-wide reference dinucleotide composition, raw (un-folded, 16
    contexts) resolution — analogous to calculate_ref_trinuc for SBS.
    Folding to the 78 canonical DBS78 channels is deferred to estimation
    time (_dbs78_ref_weighted below), mirroring how calculate_ref_trinuc's
    64 raw contexts only get folded to 32 inside estimate_96.

    Derived on the fly from ref.h5 by pairing each position with its
    successor (same trick as calculate_ref_indel100's next_ref_base_arr)
    — there's no separate dinuc.h5 index file, since this is cheap enough
    to compute per chromosome at estimate time.
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
        # Exclude any pair touching an N/out-of-range base (base2num
        # encodes those as 4) — same convention as calculate_ref_indel100's
        # ref_base_safe and indel100_reference_bucket_indices.
        valid = (base1 <= 3) & (base2 <= 3)
        dinuc_idx = 4 * base1[valid] + base2[valid]
        dinuc_count += np.bincount(dinuc_idx, minlength=16)
    return dinuc_count


def calculate_dbs_opportunity(coverage_bed_path, ref_h5, regions):
    """Genome-wide observed DBS78 opportunity (raw-144 resolution — fold
    with combine_raw_dbs_to_dbs78 for the 78 canonical channels), derived
    from the existing per-locus SBS coverage bed.

    NOT called by do_estimate's main path anymore — it now sums
    _dbs_by_duplex_group.txt instead (call.py's _compute_dbs_opportunity),
    which gives true per-duplex-group resolution instead of this
    function's single-base-independence approximation, and enables the
    group-size-stratified DBS burden breakdown (estimate_dbs78_by_group).
    Kept here as a working, self-contained alternative (needs no call.py
    output beyond the coverage bed) in case it's useful for
    cross-checking or a bed file predating _dbs_by_duplex_group.txt.

    For every pair of directly ADJACENT covered positions (start,
    start+1) with a valid A/T/C/G reference base at both, and for each of
    the 9 possible simultaneous-both-positions alt combinations,
    opportunity is approximated as coverage(alt1 at position 1) *
    coverage(alt2 at position 2) — treating the two positions' single-base
    calling processes as independent. This is a simplifying approximation
    (adjacent positions' error/damage processes aren't always
    independent — correlated adjacent damage, e.g. UV-induced CC>TT, is
    part of what a real DBS signature captures), but it's the standard
    way to approximate DBS opportunity from per-base coverage alone, and
    it automatically inherits every mask (SNP/noise/coverage/mapq/trim)
    already baked into the SBS coverage bed's columns.

    For every pair of directly ADJACENT covered positions (start,
    start+1) with a valid A/T/C/G reference base at both, and for each of
    the 9 possible simultaneous-both-positions alt combinations,
    opportunity is approximated as coverage(alt1 at position 1) *
    coverage(alt2 at position 2) — treating the two positions' single-base
    calling processes as independent. This is a simplifying approximation
    (adjacent positions' error/damage processes aren't always
    independent — correlated adjacent damage, e.g. UV-induced CC>TT, is
    part of what a real DBS signature captures), but it's the standard
    way to approximate DBS opportunity from per-base coverage alone, and
    it automatically inherits every mask (SNP/noise/coverage/mapq/trim)
    already baked into the SBS coverage bed's columns.

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

        # 4x4 reference-base combos (N/out-of-range ref bases encode as 4,
        # so they're automatically excluded — ref1_num/ref2_num only ever
        # range 0-3 here) x their 9 non-ref alt combos each: every
        # combination is a vectorized sum over the whole chromosome's
        # adjacent-pair arrays, not a per-position Python loop.
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
    """74 canonical (strand-symmetry-folded) indel classification labels,
    the estimation-time counterpart of build_indel100_labels (misc.py),
    mirroring how SBS folds 192 raw classes to 96 via combine_raw192_to_sbs96.
    (Name kept as "indel76" despite the actual count, same historical
    convention as "indel100"/"indel83" — see build_indel100_labels.)

    A homopolymer run of Cs on the reference strand is the same physical
    locus as a run of Gs read from the complementary strand, so
    delC_hp{n}/delG_hp{n} are one mutational process, not two — same logic
    pairs A with T. Homopolymer deletion and insertion are pooled this way
    (4 bases x N length bins -> 2 pooled bases x N length bins, for each of
    del/ins — deletion has 6 length bins, insertion only 5: no "0"/rep0 bin
    here, see INDEL100_INS_HP_LEN_BINS). STR deletion/insertion (unit size
    has no base identity to pool) and microhomology deletion pass through
    unchanged.
    6+6 (del CG/AT) + 5+5 (ins CG/AT) + 24 (STR del) + 24 (STR ins) + 4 (MH)
    = 74. Deterministic order, must not be built independently elsewhere.
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
    """Fold the 100 raw indel classes (build_indel100_labels) into the 74
    canonical indel76 classes (build_indel76_labels above), by summing the
    strand-symmetric homopolymer base pairs (C+G, A+T) for deletion and
    insertion; STR and microhomology classes pass through unchanged. This
    is the single point where base-pooling happens, deferred to estimation
    time rather than folded at raw accumulation time — mirrors
    combine_raw192_to_sbs96 for SBS.

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
    """(ref_counts/ref_counts.sum()) / (obs_counts/obs_counts.sum()), i.e.
    how much to scale up each class's mutation count to correct for that
    class being under/over-represented in observed coverage relative to
    the reference genome — but any class with zero observed coverage gets
    ratio 0 instead of Inf/NaN. There's no data to correct for an
    unobserved class, and its mutation count is necessarily 0 too (you
    can't observe a mutation with no coverage), so 0 * 0 = 0 correctly
    contributes nothing to a corrected-burden sum — whereas the unguarded
    division produces Inf there, and Inf * 0 is NaN under IEEE754, which
    poisons the whole sum for a single unobserved class. If the whole
    stratum has zero observed coverage, every class gets ratio 0."""
    obs_sum = obs_counts.sum()
    if obs_sum <= 0:
        return np.zeros_like(obs_counts, dtype=float)
    ref_frac = ref_counts / ref_counts.sum()
    obs_frac = obs_counts / obs_sum
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(obs_frac > 0, ref_frac / obs_frac, 0.0)


def estimate_96(trinuc_cov_96_by_rf, trinuc_mut_by_rf, ref_trinuc_64, n):
    """Estimate SBS96 mutation rates with per-class correction ratios.

    trinuc_cov_96_by_rf : (96, n_groups) — per-SBS96-class L-weighted
        coverage, already combined from the 192 raw (un-folded) classes via
        combine_raw192_to_sbs96.
    trinuc_mut_by_rf : (96, n_groups) — observed mutation counts per SBS96 class.
    ref_trinuc_64    : (64,) — trinucleotide counts in the reference genome,
        raw (un-folded) resolution; folded to the 32 canonical contexts here
        via the i/i+32 reverse-complement pairing.

    Returns burden curves both cumulative ("at least this group size",
    nmin>=1..5, the original behavior) and exact ("at exactly this group
    size", nmin==1..4 plus a pooled nmin>=5 bin) -- see the _exact-suffixed
    return values.
    """
    print("........Estimating mutation rate for each trinucleotide context.......")
    ref_trinuc_32 = ref_trinuc_64[:32] + ref_trinuc_64[32:64]

    ref_trinuc_96 = np.repeat(
        ref_trinuc_32, 3
    )  # (96,) — same ref count for 3 alts of each trinuc

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
    )  # (96,) — alt-base-specific
    mutnum = trinuc_mut.sum()
    cov = trinuc_cov.sum() / 3  # divide by 3 to get per-locus-equivalent coverage
    burden_uncorrected[4] = mutnum / cov
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = poisson_confint(mutnum, cov)
    trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov)
    # Alt-base-specific correction: each of the 96 classes gets its own ratio
    correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov)
    mutnum_corrected = correction_ratio * trinuc_mut
    burden_corrected[4] = mutnum_corrected.sum() / cov
    burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
        mutnum_corrected.sum(), cov
    )
    hap_trinuc[:, 4] = trinuc_rate * ref_trinuc_96
    covs[4] = cov

    for nn in range(4, 0, -1):
        trinuc_mut = trinuc_mut + trinuc_mut_by_rf[:, nmin == nn].sum(axis=1)
        trinuc_cov = trinuc_cov + trinuc_cov_96_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = trinuc_mut.sum()
        cov = trinuc_cov.sum() / 3
        burden_uncorrected[nn - 1] = mutnum / cov
        burden_uncorrected_lb[nn - 1], burden_uncorrected_ub[nn - 1] = poisson_confint(
            mutnum, cov
        )
        trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov)
        covs[nn - 1] = cov
        correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov)
        mutnum_corrected = trinuc_mut * correction_ratio
        burden_corrected[nn - 1] = mutnum_corrected.sum() / cov
        burden_corrected_lb[nn - 1], burden_corrected_ub[nn - 1] = poisson_confint(
            mutnum_corrected.sum(), cov
        )
        hap_trinuc[:, nn - 1] = trinuc_rate * ref_trinuc_96

    # Exact (non-cumulative) group-size bins: nmin==1,2,3,4, and nmin>=5 for
    # the last bin (5+ groups are pooled, same as the cumulative curve's
    # bottom bin) -- "burden AT this exact group size" rather than the
    # cumulative curve's "burden at LEAST this group size" above. Each bin
    # is independent (no running total carried across bins), unlike the
    # cumulative loop.
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
        ) = poisson_confint(mutnum_e, cov_e)
        correction_ratio_e = _safe_correction_ratio(ref_trinuc_96, trinuc_cov_e)
        mutnum_corrected_e = (trinuc_mut_e * correction_ratio_e).sum()
        burden_corrected_exact[nn - 1] = (
            mutnum_corrected_e / cov_e if cov_e > 0 else 0.0
        )
        (
            burden_corrected_exact_lb[nn - 1],
            burden_corrected_exact_ub[nn - 1],
        ) = poisson_confint(mutnum_corrected_e, cov_e)

    # A channel with zero observed coverage (trinuc_cov, after subtracting
    # each mutant locus from its own class's coverage above) already gets
    # correction_ratio 0 (_safe_correction_ratio) -- but its raw
    # (uncorrected) mutnum can still be nonzero, e.g. a class whose only
    # covered locus is itself the mutant one, subtracted away to 0. Force
    # it to 0 too, consistent with correction_ratio, so
    # "mutation_number_uncorrected" never shows a mutation count with no
    # observed genome behind it.
    # 0 (not 0.0): trinuc_mut is int-dtype (mutation counts); np.where with
    # a float fill would upcast the whole array, breaking downstream
    # consumers (e.g. Summarize.py's int() parse of "Uncorrected mutation
    # number") that expect an integer-formatted total.
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


# Number of microhomology-length sub-channels folded under each deletion-
# length group ("2","3","4","5+") in the 83-channel ID83 scheme: a length-2
# deletion can only have mh_len==1 (mh_len<del_len by definition of "not a
# repeat-unit deletion"), length-3 gives {1,2}, length-4 gives {1,2,3}, and
# the "5+" group (true del_len>=5) gives {1,2,3,4,5+} — 1+2+3+5=11 channels.
INDEL83_MH_GROUP_SIZES = [1, 2, 3, 5]


def build_indel83_labels():
    """83-channel ID83-style indel classification — the estimation-time
    counterpart of funcs/misc.py's classify_indel_channel (which resolves a
    single observed call to one of these exact label strings), and the
    successor to build_indel76_labels above now that classification is
    sequence-based rather than repeat-annotation-based. Matches the
    standard COSMIC ID83 scheme:
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
    """Map build_indel83_labels' label strings (our internal ID83 syntax,
    e.g. "Cdelhp1", "2delstr3", "5+delMH2") onto SigProfilerPlotting's
    COSMIC-standard ID83 label strings (e.g. "1:Del:C:0", "2:Del:R:2",
    "5:Del:M:2") — confirmed identical grouping/order against
    site-packages/sigProfilerPlotting/reference_formats/ID83.txt, so this
    mirrors build_indel83_labels' own loop structure (same base/unit-size
    groups, same length/repeat-count bin order) rather than hardcoding an
    83-entry lookup table that could silently drift out of sync with it.

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
    """Fold the 74-canonical-class indel vector (combine_indel100_to_indel76)
    into the 83-channel ID83 scheme, by replicating each of the 4
    unstratified delMH_len{2,3,4,5+} values across its 1/2/3/5
    microhomology-length sub-channels. The coverage side (_indel_by_
    duplex_group.txt, still written by the old hp.h5-annotation-based
    call.py) has no finer-than-deletion-length resolution for
    microhomology — only the observed-mutation side gets true per-mh-
    length resolution now, via classify_indel_channel — so every
    sub-channel in a group shares its parent's coverage value verbatim.

    Del-hp (12), STR del (24) and STR ins (24) are identical between the
    two schemes and pass through unchanged. Ins-hp is not: indel76 has no
    "0"/rep0 bin at all (INDEL100_INS_HP_LEN_BINS excludes it — that
    opportunity is always thrown away and replaced anyway, see
    override_inshp0_with_next_base_opportunity below), so ID83's two hp0
    slots (Cinshp0 index 12, Tinshp0 index 18) are left at 0 here; the
    caller must run override_inshp0_with_next_base_opportunity afterward
    to fill them. The remaining 10 ins-hp values slot into ID83's other 10
    ins-hp positions around those two gaps.

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
    """Fill the ID83 Cinshp0/Tinshp0 entries (indices 12 and 18) — left at
    0 by expand_indel76_to_indel83, since indel76/indel100 carry no hp0
    bin at all — with a computation based purely on the reference's own
    next-base identity: funcs/misc.py's 1bpins{base} columns, indices
    96-99 of the raw 100-class scheme (build_indel100_labels), rather than
    hp.h5 repeat annotation. rep0 for inserting base N requires only that
    the following reference base differs from N; it has nothing to do with
    what's actually at the insertion site itself (unlike every other
    ID83 channel, which does depend on the local repeat context).

    Column 96+b (b=A/T/C/G, base2num order) is already exactly "the
    opportunity to insert base b with rep0" at both resolutions: unweighted
    position counts for ref_indel100 (indel100_reference_bucket_indices),
    and L-table power using base b itself (not the position's own
    reference base) for indel100_by_rf_np (funcs/call.py's per-family
    accumulation, cols 6-9 of cov_mat_indel). So no further transformation
    is needed here beyond the usual C+G / A+T pooling every other ID83
    homopolymer channel already gets.

    values_83  : (83,) or (83, n_groups) ID83-resolution array (ref
        composition or coverage), indices 12/18 still 0 from
        expand_indel76_to_indel83 — not modified in place.
    values_100 : (100,) or (100, n_groups) raw array at the same
        resolution as indel100_by_rf_np/ref_indel100.
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
    """Per-channel correction ratio (ref_frac/obs_frac) for the first 72
    ID83 channels, but computed per deletion-length GROUP — not per
    individual microhomology-length sub-channel — for the last 11: the 4
    delMH groups (channels 72-82) each share one ratio across their
    sub-channels. This mirrors the coverage side's actual resolution
    (expand_indel76_to_indel83 replicates one value per group already) and
    is also just correct statistics — computing a separate ratio per
    individual mh sub-channel would split already-thin counts even
    further for no benefit, since the coverage data can't distinguish them
    anyway.

    Channels (individual or grouped) with zero observed fraction get
    ratio 0 instead of Inf/NaN — see _safe_correction_ratio's docstring
    for why that's the correct fallback rather than letting 0/0 or
    nonzero/0 propagate into the corrected-burden sum.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(obs_frac > 0, ref_frac / obs_frac, 0.0)
    for sl in _indel83_mh_group_slices():
        obs_group = obs_frac[sl].sum()
        ratio[sl] = ref_frac[sl].sum() / obs_group if obs_group > 0 else 0.0
    return ratio


def estimate_indel83(indel83_cov_by_rf, indel83_mut_by_rf, ref_indel83, n):
    """ID83-resolution counterpart of estimate_indel76 above: identical
    read-depth-stratified (nmin) burden/correction machinery, except the
    correction ratio for the 11 microhomology channels (72-82) is grouped
    per deletion length rather than computed per individual channel — see
    _grouped_indel83_correction_ratio.
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
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = poisson_confint(mutnum, cov)
    indel_rate = _safe_rate(indel_mut, indel_cov)
    with np.errstate(divide="ignore", invalid="ignore"):
        obs_frac = indel_cov / indel_cov.sum()
    correction_ratio = _grouped_indel83_correction_ratio(ref_frac, obs_frac)
    mutnum_corrected = correction_ratio * indel_mut
    burden_corrected[4] = mutnum_corrected.sum() / cov if cov > 0 else float("nan")
    burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
        mutnum_corrected.sum(), cov
    )
    hap_indel83[:, 4] = indel_rate * ref_indel83
    covs[4] = cov

    for nn in range(4, 0, -1):
        indel_mut = indel_mut + indel83_mut_by_rf[:, nmin == nn].sum(axis=1)
        indel_cov = indel_cov + indel83_cov_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = indel_mut.sum()
        cov = indel_cov.sum()
        burden_uncorrected[nn - 1] = mutnum / cov if cov > 0 else float("nan")
        burden_uncorrected_lb[nn - 1], burden_uncorrected_ub[nn - 1] = poisson_confint(
            mutnum, cov
        )
        indel_rate = _safe_rate(indel_mut, indel_cov)
        covs[nn - 1] = cov
        with np.errstate(divide="ignore", invalid="ignore"):
            obs_frac = indel_cov / indel_cov.sum()
        correction_ratio = _grouped_indel83_correction_ratio(ref_frac, obs_frac)
        mutnum_corrected = indel_mut * correction_ratio
        burden_corrected[nn - 1] = (
            mutnum_corrected.sum() / cov if cov > 0 else float("nan")
        )
        burden_corrected_lb[nn - 1], burden_corrected_ub[nn - 1] = poisson_confint(
            mutnum_corrected.sum(), cov
        )
        hap_indel83[:, nn - 1] = indel_rate * ref_indel83

    # Exact (non-cumulative) group-size bins -- see estimate_96's matching
    # block for the cumulative-vs-exact distinction; identical structure,
    # just using this function's own grouped microhomology correction
    # ratio instead of the per-channel one.
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
        ) = poisson_confint(mutnum_e, cov_e)
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
        ) = poisson_confint(mutnum_corrected_e, cov_e)

    # See estimate_96's matching comment: a channel with zero observed
    # coverage (indel_cov, after subtracting each mutant locus from its own
    # class's coverage upstream) already gets correction_ratio 0, but its
    # raw mutnum can still be nonzero -- force it to 0 too.
    # 0, not 0.0 -- see estimate_96's matching comment (preserve int dtype).
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
    """Genome-wide reference-composition opportunity for each of the 78
    canonical DBS78 channels.

    A flat per-ref broadcast (repeating one ref-level count across every
    alt-channel of that ref) is correct for the 6 non-self-complementary
    refs: every one of their 9 canonical alt channels receives
    contributions from exactly 2 raw (ref,alt) combinations (one from the
    ref's own reading, one from its distinct RC-partner ref's reading) —
    uniform, so a flat broadcast matches. It is NOT correct for the 4
    self-complementary refs (AT, TA, CG, GC — see build_dbs78_labels): a
    single locus's own 9 possible alts do not map uniformly onto that
    ref's 6 canonical channels there. 3 alts are self-inverse
    (RC(alt)==alt) and fold 1:1 onto their own channel, while the other 6
    pair up under reverse-complementing into 3 merged channels, each
    receiving 2 of the ref's own 9 alts. This weighting must match
    dbs_opp_78/dbs_mut_78 (which inherit this same 2x/1x split via
    combine_raw_dbs_to_dbs78), or correction_ratio picks up a spurious
    ~2x split confined to exactly these 4 ref groups.

    Modeled by crediting each raw ref's full genome count to every one of
    its OWN 9 possible raw alts (dinuc_count_16[ref] each — the correct
    null model for "no alt-specific bias": absent any information to the
    contrary, every possible mutation at a given ref-dinucleotide type is
    equally likely), then folding with the same combine_raw_dbs_to_dbs78
    used for observed mutations and opportunity — so ref/obs/mut all go
    through identical folding logic and are guaranteed self-consistent by
    construction, rather than two independent approximations happening to
    line up (or not). No normalizing division is applied here (e.g. by the
    9 raw alts) — every one of the 144 raw entries gets the same
    unscaled ref count, so any such constant would be a uniform multiplier
    across all 78 folded channels and cancel out identically in
    correction_ratio's ref_frac = ref_counts/ref_counts.sum() anyway;
    leaving it out just keeps this column's numbers whole and directly
    interpretable as genome dinucleotide counts, rather than introducing a
    fraction that carries no additional information.

    dinuc_count_16: (16,) raw reference dinucleotide counts, as returned
        by calculate_ref_dbs (un-folded).
    Returns (78,) canonical-channel reference opportunity.
    """
    raw144 = np.zeros(144)
    for label, idx in DBS_RAW144_LABEL2NUM.items():
        raw_ref = label.split(">")[0]
        raw144[idx] = dinuc_count_16[_dinuc_num(raw_ref)]
    return combine_raw_dbs_to_dbs78(raw144, DBS_RAW144_LABEL2NUM)


def estimate_dbs78(dbs_mut_78, dbs_opp_78, ref_dbs_78):
    """Estimate DBS78 mutation rate and per-class correction ratios: one
    whole-sample figure, no group-size breakdown (see estimate_dbs78_by_group
    below for that).

    dbs_mut_78 : (78,) observed DBS78 mutation counts (dbs_raw_event_to_label
        applied to every PASS record in {sample}_dbs.vcf).
    dbs_opp_78 : (78,) DBS78 opportunity, folded with combine_raw_dbs_to_dbs78
        from _dbs_by_duplex_group.txt summed across all duplex groups.
    ref_dbs_78 : (78,) reference-genome dinucleotide composition, folded
        from the 16 raw ref groups via _dbs78_ref_weighted.
    """
    print("........Estimating mutation rate for each DBS78 context.......")
    # Every raw (un-folded) reference dinucleotide has exactly 9 possible
    # alt combinations (build_dbs_raw144_labels), regardless of whether the
    # ref is self-complementary -- folding 144->78 only regroups entries,
    # it never changes their sum (combine_raw_dbs_to_dbs78 is a pure
    # reindexing sum) -- so dividing the folded opportunity total by 9
    # recovers the same per-locus-pair-equivalent coverage that dividing
    # the raw144 total by 9 would, without needing per-channel handling of
    # the varying (6 or 9) group widths.
    cov = dbs_opp_78.sum() / 9
    mutnum = dbs_mut_78.sum()
    burden_uncorrected = mutnum / cov if cov > 0 else float("nan")
    burden_uncorrected_lb, burden_uncorrected_ub = poisson_confint(mutnum, cov)

    correction_ratio = _safe_correction_ratio(ref_dbs_78, dbs_opp_78)
    mutnum_corrected = correction_ratio * dbs_mut_78
    burden_corrected = mutnum_corrected.sum() / cov if cov > 0 else float("nan")
    burden_corrected_lb, burden_corrected_ub = poisson_confint(
        mutnum_corrected.sum(), cov
    )

    # See estimate_96's matching comment: a channel with zero observed
    # opportunity already gets correction_ratio 0, but its raw mutnum can
    # still be nonzero -- force it to 0 too.
    # 0, not 0.0 -- see estimate_96's matching comment (preserve int dtype).
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
    """Group-size-stratified DBS78 burden — the DBS counterpart of
    estimate_96/estimate_indel83's cumulative ("at least this group size")
    and exact ("at exactly this group size") curves, using
    _dbs_by_duplex_group.txt (call.py's _compute_dbs_opportunity) for
    per-duplex-group DBS opportunity resolution.

    dbs_mut_by_rf    : (78, n_groups) observed DBS78 mutation counts per
        duplex group (each PASS _dbs.vcf record's F1R2/F2R1 info fields
        resolve which group column it belongs to).
    dbs_opp144_by_rf : (144, n_groups) raw-144 DBS opportunity per duplex
        group, straight from _dbs_by_duplex_group.txt — folded to 78
        per-stratum below (each stratum sums its own opportunity slice
        before folding, mirroring how ref_trinuc_96/covs are handled per
        stratum in estimate_96).
    ref_dbs_78       : (78,) reference-genome dinucleotide composition.
    n                : duplex-group column labels ("F1R2+F2R1"), same order
        as both matrices' columns (and as trinuc_by_rf.columns, by
        construction — see Caller.py's shared non_zero_keys writer order).
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
        u_lb, u_ub = poisson_confint(mutnum, cov)
        correction_ratio = _safe_correction_ratio(ref_dbs_78, opp78)
        mutnum_corrected = (correction_ratio * mut78).sum()
        corrected = mutnum_corrected / cov if cov > 0 else 0.0
        c_lb, c_ub = poisson_confint(mutnum_corrected, cov)
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
    ## Calculate burden when take 10 as mininum
    trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = muts_by_rf[:, nmin >= 5].sum(axis=1).sum()
    cov = trinuc_cov.sum() / 3
    burden_indel[4] = mutnum / cov
    burden_indel_lb[4], burden_indel_ub[4] = poisson_confint(mutnum, cov)
    # trinuc_rate[:,9] = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
    for nn in range(4, 0, -1):
        trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= nn].sum(axis=1)
        mutnum = muts_by_rf[:, nmin >= nn].sum(axis=1).sum()
        cov = trinuc_cov.sum() / 3
        burden_indel[nn - 1] = mutnum / cov
        burden_indel_lb[nn - 1], burden_indel_ub[nn - 1] = poisson_confint(mutnum, cov)

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
    """Write the combined burden-by-group-size outputs shared across SBS,
    INDEL, and DBS: one wide-format table with parallel min_group_size
    (cumulative "at least group size N", the original per-type behavior)
    and group_size (exact "at group size N", non-cumulative) column
    blocks sharing one read-number (1-5) axis, plus one 2-panel PDF (one
    panel per breakdown) each showing corrected+uncorrected burden with
    95% CI. Every *_uncorrected/*_corrected/*_cov argument is a (5,) array
    indexed by read-number-1 (bins 1,2,3,4,5+ of min(F1R2,F2R1) reads),
    already computed by estimate_96/estimate_indel83/estimate_dbs78_by_group.

    mutation_type: lowercase file-name infix, e.g. "sbs"/"indel"/"dbs" --
    output files are {outdir}/{sample}_{mutation_type}_burden_by_group_size.{txt,pdf}.
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


def do_estimate(args):
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
        # Raw-144 DBS opportunity per duplex group (call.py's
        # _compute_dbs_opportunity), columns in the same order as
        # trinuc_by_rf.columns by construction (Caller.py writes both from
        # the same non_zero_keys list) -- used below for genome-wide DBS
        # opportunity and for the group-size-stratified DBS burden
        # breakdown.
        dbs144_by_rf = pd.read_csv(
            dbs_dir + "/" + sample + "_dbs_by_duplex_group.txt", sep="\t", index_col=0
        )
        dbs144_by_rf_np = dbs144_by_rf.to_numpy().astype(float)
        ###Estimate SNV burden
        vcf = VCF(sbs_dir + "/" + sample + "_sbs.vcf", "r")
        mut_by_rf = dict()
        # trinuc_by_rf is the raw (un-folded) 192-row file written by
        # Caller.py: 64 trinuc contexts x 3 non-ref alt bases, labeled
        # "{trinuc}>{alt}". Reverse-complement pairs are combined into the
        # 96 canonical classes later, via combine_raw192_to_sbs96.
        label2num_192 = {label: i for i, label in enumerate(trinuc_by_rf.index)}
        trinuc_by_rf_np = trinuc_by_rf.to_numpy().astype(float)
        trinuc_mut_np = np.zeros([96, len(trinuc_by_rf.columns)], dtype=int)
        duplex_no_dict = dict()
        revcomp = _REVCOMP
        for nn, duplex_no in enumerate(trinuc_by_rf.columns):
            duplex_no_dict[duplex_no] = nn

        # Dictionary to track unique mutations and their counts
        unique_mutations = {}

        if args.dilute:
            vcf_out = VCF(
                sbs_dir + "/" + sample + "_sbs_flt.vcf", "w", header=vcf.header
            )
        stats = parse_stats_file(f"{prefix}/{sample}_stats.txt")
        # "Unmasked Coverage" in stats.txt is the opportunity coverage
        # (L-weighted coverage summed over the 3 possible alt bases per
        # locus); divide by 3 to get the base (per-locus) coverage, as is
        # done for trinuc_cov elsewhere in this file.
        unmasked_cov = int(float(stats["Unmasked Coverage"])) / 3
        unmasked_indel_cov = int(float(stats["Unmasked Indel Coverage"]))
        ###Read DBS events up front: classify each PASS event into its
        ### DBS78 channel (used later, in the DBS burden section), and
        ### record its two constituent genomic positions so the SBS96
        ### counting below can skip them. Every PASS DBS event's two
        ### positions are, by construction (call.py's _detect_dbs_pairs),
        ### ALSO individually PASS-called SNVs in _sbs.vcf — counting both
        ### would double-count the same physical mutation as 1 DBS event
        ### AND 2 separate SBS mutations.
        vcf_dbs = VCF(dbs_dir + "/" + sample + "_dbs.vcf", "r")
        dbs_mut_78 = np.zeros(78, dtype=int)
        # Per-duplex-group counterpart of dbs_mut_78, same duplex_no_dict
        # column indexing as trinuc_mut_np/indel83_mut_np above -- lets the
        # DBS burden section below stratify by group size the same way SBS
        # and indel already do, using the per-group opportunity resolution
        # _dbs_by_duplex_group.txt provides.
        dbs_mut_by_rf = np.zeros([78, len(trinuc_by_rf.columns)], dtype=int)
        dbs_snv_positions = set()
        for rec in vcf_dbs.fetch():
            if "PASS" not in rec.filter:
                continue
            dbs_snv_positions.add((rec.chrom, rec.pos))
            dbs_snv_positions.add((rec.chrom, rec.pos + 1))
            label = dbs_raw_event_to_label(rec.ref, rec.alts[0])
            dbs_mut_78[DBS78_LABEL2NUM[label]] += 1
            duplex_no = f"{rec.info['F1R2']}+{rec.info['F2R1']}"
            dbs_mut_by_rf[DBS78_LABEL2NUM[label], duplex_no_dict[duplex_no]] += 1
        vcf_dbs.close()

        ###Calculate masked mutation rate
        # Unmasked burden = PASS (from _sbs.vcf) + snp/noise-masked
        # candidates that cleared LR and got real depth extracted (from
        # _sbs_fail.vcf, filter=="masked" with SNPM/NOISEM set and DP>0 --
        # see call.py's emission logic and Caller.py's deferred_depth_keys).
        # Any other "masked"/rescue-label record in the fail VCF (below LR
        # threshold, or a non-snp/noise mask without --rescue) has no real
        # depth and must not count.
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
            if (rec.chrom, rec.pos) in dbs_snv_positions:
                continue
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
            if (rec.chrom, rec.pos) in dbs_snv_positions:
                continue
            unmasked_mut_count += 1
        vcf_fail.close()
        unmasked_sbs_burden = (
            float(unmasked_mut_count) / float(unmasked_cov)
            if unmasked_cov > 0
            else float("nan")
        )
        unmasked_sbs_burden_lb, unmasked_sbs_burden_ub = poisson_confint(
            unmasked_mut_count, unmasked_cov
        )
        snv_progress = {"start": time.time(), "last": time.time()}
        for i, rec in enumerate(vcf.fetch()):
            log_progress(
                snv_progress,
                "Processing SNV mutations",
                i,
                extra=f"{rec.chrom}:{rec.pos}",
            )
            count_flag = 0
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
                count_flag = 1
            unique_mutations[mutation_key][0] += 1
            """
            if count_flag == 0:
                continue
            """
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
            # duplex_no = str(min(F1R2, F2R1)) + "+" + str(max(F1R2, F2R1))
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
            # Skip the SBS96 numerator for a DBS-forming position (already
            # counted once in dbs_mut_78 above) — but still subtract it from
            # the raw (un-folded) (trinuc, alt) coverage cell below: the
            # locus genuinely isn't non-mutant coverage either way, it's
            # just categorized as a DBS event rather than an SBS one.
            if (chrom, pos) not in dbs_snv_positions:
                trinuc_mut_np[
                    TRINUCSBS2NUM_96[trinucSbs], duplex_no_dict[duplex_no]
                ] += 1
            # Subtract from the raw (un-folded) (trinuc, alt) coverage cell —
            # using the mutation's actual, as-reported orientation, not the
            # pyrimidine-folded one — so this mutant locus isn't also
            # counted as effective non-mutant coverage.
            raw_idx_192 = label2num_192[f"{raw_trinuc}>{raw_alt}"]
            trinuc_by_rf_np[raw_idx_192, duplex_no_dict[duplex_no]] -= 1

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
            trinuc_cov_96_by_rf, trinuc_mut_np, ref_trinuc, trinuc_by_rf.columns
        )
        # Total reference bases actually considered (non-N, non-noise-
        # masked) across args.regions -- ref_trinuc.sum() (all 64 raw,
        # un-folded classes) is exactly this count, since every considered
        # position is credited to exactly one of those 64 classes by its
        # own literal ref base. Reused as the shared denominator for a
        # genuine mutations-per-base rate across SBS/DBS/indel below,
        # unlike e.g. ref_indel83.sum() which double/quadruple-counts
        # positions across its many opportunity columns and so isn't a
        # literal base count.
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
        # SigProfilerPlotting renders one page per DataFrame column (COSMIC-
        # standard SBS96 spectrum plot); mutation_number_genome is dropped
        # here (genome-wide extrapolated panel not needed) but stays in the
        # _sbs_96_corrected.txt table above. Output filename is controlled
        # by the library: SBS_96_plots_{sample}.pdf. process_input() only
        # trusts an incoming DataFrame's index if its name is already one
        # of MutationType/index/Mutation Types/classification — otherwise
        # it discards the real index and tries to use the first data
        # column's own values as the index instead, which collides on
        # duplicate values (e.g. multiple zero-mutation channels).
        # rename_axis pins the index name so our real SBS96 labels survive.
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
            f.write(f"Mutation number per genome\t{corrected_trinuc_num.sum()}\n")
            # Effective coverage at min_group_size=1 (cov_by_minread[0], the
            # same cumulative "at least 1 read on each strand" figure as
            # _sbs_burden_by_group_size.txt's min_group_size=1 row) rather
            # than the reference genome's raw trinucleotide total -- that
            # total is still available below as Reference base number.
            f.write(f"Genome coverage\t{cov_by_minread[0]}\n")
            f.write(f"Unmasked burden\t{unmasked_sbs_burden}\n")
            f.write(f"Unmasked burden 95% lower\t{unmasked_sbs_burden_lb}\n")
            f.write(f"Unmasked burden 95% upper\t{unmasked_sbs_burden_ub}\n")
            f.write(f"Reference base number\t{reference_base_number}\n")

        vcf.close()
        ###Calculate indel burdens
        print("......Estimating Indel Burden.......")
        vcf_indel = VCF(f"{indel_dir}/{sample}_indel.vcf", "r")
        # indel100_by_rf is the raw 100-class coverage/opportunity file
        # written by Caller.py (build_indel100_labels, misc.py), still
        # based on hp.h5 repeat annotation — folded down to the 74-class
        # resolution (combine_indel100_to_indel76) and then expanded to the
        # 83-channel ID83 layout (expand_indel76_to_indel83) as an
        # approximation, since call.py doesn't yet compute per-locus
        # coverage at true microhomology-length resolution. Observed
        # mutations, by contrast, ARE classified at full ID83 resolution
        # below, directly from reference sequence (classify_indel_channel)
        # — so the numerator is exact while the denominator is a coarser
        # stand-in until call.py's coverage loop is rewritten to match.
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
        hp_h5 = h5py.File(args.reference + ".hp.h5", "r")
        str_h5 = h5py.File(args.reference + ".str.h5", "r")

        # Unmasked burden numerator: PASS (below, from vcf_indel) +
        # noise-masked candidates that cleared LR and got real depth
        # extracted (from _indel_fail.vcf, filter=="masked" with NOISEM set
        # and DP>0 -- see the matching SNV block above for the full
        # eligibility explanation; indels have no SNPM tag, only NOISEM).
        vcf_indel_fail = VCF(f"{indel_dir}/{sample}_indel_fail.vcf", "r")
        unmasked_indel_count = 0
        for rec in vcf_indel_fail.fetch():
            if not (
                "masked" in rec.filter
                and rec.info.get("NOISEM")
                and rec.samples["TUMOR"]["DP"] > 0
            ):
                continue
            unmasked_indel_count += 1
        vcf_indel_fail.close()

        indel_count = 0
        indel_progress = {"start": time.time(), "last": time.time()}
        for i, rec in enumerate(vcf_indel.fetch()):
            log_progress(
                indel_progress,
                "Processing indel mutations",
                i,
                extra=f"{rec.chrom}:{rec.pos}",
            )
            unmasked_indel_count += 1
            indel_count += 1
            mutation_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            TAC = rec.samples["TUMOR"]["AC"]
            TDP = rec.samples["TUMOR"]["DP"]

            # Classify by actual sequence content (classify_indel_channel,
            # funcs/misc.py) as the ground truth for category (repeat-unit
            # vs microhomology), with repeat-count read straight from
            # hp.h5 (homopolymer events) or str.h5 (STR events, wherever
            # its unit_len agrees with what was actually observed) — see
            # classify_indel_channel's `anno` docstring for why that's
            # safe, and why it isn't trusted for the category call itself.
            # rec.pos is pysam's 1-based
            # VCF POS (the anchor base), which numerically equals the
            # 0-based index of the first affected base in ref.h5 — see
            # funcs/misc.py's classify_indel_record docstring for why. That
            # anchor is always left-aligned (funcs/indels.py's
            # left_align_indel, applied at call time before this VCF was
            # ever written), which is also why only ref_after is needed
            # here — see the same docstring.
            F1R2 = rec.info["F1R2"]
            F2R1 = rec.info["F2R1"]
            duplex_no = str(F1R2) + "+" + str(F2R1)
            indel_len = len(rec.alts[0]) - len(rec.ref)
            chrom_len = ref_h5[rec.chrom].shape[0]
            if indel_len < 0:
                del_len = -indel_len
                indel_seq = _decode_ref_seq(
                    ref_h5[rec.chrom][rec.pos : rec.pos + del_len]
                )
                after_start = rec.pos + del_len
            else:
                indel_seq = rec.alts[0][1:].upper()
                after_start = rec.pos
            after_end = min(chrom_len, after_start + INDEL_CONTEXT_WINDOW)
            ref_after = _decode_ref_seq(ref_h5[rec.chrom][after_start:after_end])
            # Annotation lookup at rec.pos: for a deletion that's the first
            # deleted base (start of the deleted unit's own tract, since
            # left-aligned); for an insertion that's the next reference
            # base after the insertion point ("1bp next to the insertion
            # locus") — both are exactly where the pre-existing repeat
            # tract's own unit_len/repeat_count apply, matching the
            # opportunity side's own convention (call.py). Read raw from
            # hp.h5 and str.h5 independently (not load_repeat_context's
            # STR-priority merge) since classify_indel_channel picks
            # between them itself, keyed on the observed event's own unit
            # length rather than on which annotation wins positionally.
            anno = (
                int(hp_h5[rec.chrom][0, rec.pos]),
                int(str_h5[rec.chrom][0, rec.pos]),
                int(str_h5[rec.chrom][1, rec.pos]),
            )
            label = classify_indel_channel(indel_seq, ref_after, indel_len, anno)
            duplex_idx = duplex_no_dict[duplex_no]
            indel83_mut_np[label2num_83[label], duplex_idx] += 1
            if not unique_mutations.get(mutation_key):
                # 4th element: this event's own coverage.bed.gz column
                # index (0-13) -- lets duplex_allele_counts.txt read this
                # mutation's actual effective coverage at ITS alt/category,
                # not a sum across all 14 unrelated indel categories.
                unique_mutations[mutation_key] = [
                    0,
                    TAC,
                    TDP,
                    indel_coverage_category_index(indel_seq, ref_after, indel_len),
                ]
            unique_mutations[mutation_key][0] += 1
        unmasked_indel_burden = (
            unmasked_indel_count / unmasked_indel_cov
            if unmasked_indel_cov > 0
            else float("nan")
        )
        unmasked_indel_burden_lb, unmasked_indel_burden_ub = poisson_confint(
            unmasked_indel_count, unmasked_indel_cov
        )

        print("......Estimating indel83 profile and corrected burden........")
        indel76_cov_by_rf = combine_indel100_to_indel76(
            indel100_by_rf_np, label2num_100
        )
        ref_indel100 = calculate_ref_indel100(args)
        ref_indel76 = combine_indel100_to_indel76(ref_indel100, label2num_100)
        indel83_cov_by_rf = expand_indel76_to_indel83(indel76_cov_by_rf)
        ref_indel83 = expand_indel76_to_indel83(ref_indel76)
        # Cinshp0/Tinshp0 (rep0 insertion) depend only on the reference's
        # own next-base identity, not repeat annotation — filled in here
        # from the raw 100-class 1bpins{base} totals, since
        # expand_indel76_to_indel83 leaves them at 0 (indel76/indel100
        # carry no hp0 bin at all to expand from).
        indel83_cov_by_rf = override_inshp0_with_next_base_opportunity(
            indel83_cov_by_rf, indel100_by_rf_np
        )
        ref_indel83 = override_inshp0_with_next_base_opportunity(
            ref_indel83, ref_indel100
        )
        # Subtract each observed mutant locus from its own channel's
        # coverage so it isn't also counted as effective non-mutant
        # coverage — mirrors the same step for SNVs above, just applied
        # after folding/expansion since the mutation side and the coverage
        # side are no longer classified at the same raw resolution.
        indel83_cov_by_rf = indel83_cov_by_rf - indel83_mut_np
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
            indel83_cov_by_rf, indel83_mut_np, ref_indel83, indel100_by_rf.columns
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
        # SigProfilerPlotting's ID83 convention uses different label syntax
        # ("1:Del:C:0" vs our "Cdelhp1") for the same 83 channels in the
        # same order — INDEL83_TO_SIGPROFILER_LABELS translates the index
        # before handing off. mutation_number_genome is dropped from the
        # plot (genome-wide extrapolated panel not needed) but stays in the
        # _indel_83_corrected.txt table above. Output filename is
        # controlled by the library: ID_83_plots_{sample}.pdf.
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
        # Every burden below must be reported per surveyed (duplex-covered)
        # reference base, matching the SBS/DBS convention -- there, dividing
        # the raw coverage total by its known fixed per-locus multiplier (3
        # alt bases for SBS, 9 alt-pair combos for DBS, baked into
        # estimate_96/estimate_dbs78(_by_group)'s own `cov`) already yields a
        # literal per-locus count, so mutations/cov there IS mutations-per-
        # base. Indel coverage has no such fixed multiplier:
        # INDEL_COVERAGE_CATEGORY_LABELS sums coverage across overlapping,
        # locus-varying opportunity columns (one locus can count toward
        # several ID83 channels at once), so estimate_indel83's own burden/
        # cov arrays -- unlike estimate_96/estimate_dbs78_by_group's -- are
        # mutations/coverage per opportunity-COLUMN, not per locus.
        # genome_indel_cov/reference_base_number is that same column-per-
        # locus multiplicity computed from the reference genome itself (a
        # pure sequence-composition ratio, independent of sequencing depth)
        # -- assuming covered loci have the same average multiplicity as the
        # region overall, multiplying by it converts the opportunity-column
        # rate into the same per-covered-base rate SBS/DBS already report
        # directly. Computed here (rather than only at the _indel_burden.txt
        # write below) so it can also rescale the group-size arrays going
        # into write_burden_by_group_size -- those come straight from
        # estimate_indel83 and would otherwise carry the same un-rescaled
        # per-opportunity-column units into _indel_burden_by_group_size.txt.
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
        # Field order/naming mirrors _sbs_burden.txt exactly (uncorrected
        # block, corrected block, per-genome number, opportunity/coverage,
        # ref base number) — the first 4 lines are read by Summarize.py by
        # fixed line number, and still correspond to the uncorrected
        # burden/mutation count, so no Summarize.py index changes are
        # needed.
        with open(indel_dir + "/" + sample + "_indel_burden.txt", "w") as f:
            # Uncorrected burden must share the exact same (ID83-resolution)
            # coverage denominator as Corrected burden below -- both are
            # "mutations at this same population / this same opportunity",
            # differing only in whether the numerator got per-channel
            # correction-ratio reweighting. indel_burden_uncorrected[0]
            # (from estimate_indel83, same array Corrected reads index [0]
            # of) is that consistent denominator: it's expanded across the
            # 83 overlapping microhomology/repeat channels the same way
            # indel83_cov_by_rf is, unlike a flat indel_count / indel_cov
            # ratio (which is not expanded the same way, and reads too
            # high).
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
            f.write(f"Mutation number per genome\t{corrected_indel_num.sum()}\n")
            # Effective coverage at min_group_size=1, rescaled to the same
            # per-reference-locus units as the burden values above (see
            # indel_locus_multiplier comment above) -- matches
            # _indel_burden_by_group_size.txt's min_group_size=1
            # Coverage_base -- rather than the reference genome's raw
            # opportunity-column total (Reference base number below).
            f.write(
                f"Genome coverage\t"
                f"{indel_cov_by_minread[0] / indel_locus_multiplier}\n"
            )
            f.write(
                f"Unmasked burden\t{unmasked_indel_burden * indel_locus_multiplier}\n"
            )
            f.write(
                f"Unmasked burden 95% lower\t"
                f"{unmasked_indel_burden_lb * indel_locus_multiplier}\n"
            )
            f.write(
                f"Unmasked burden 95% upper\t"
                f"{unmasked_indel_burden_ub * indel_locus_multiplier}\n"
            )
            f.write(f"Reference base number\t{reference_base_number}\n")

        ###Calculate DBS burden
        print("......Estimating DBS78 Burden.......")
        # Genome-wide DBS opportunity, summed straight from
        # _dbs_by_duplex_group.txt (call.py's _compute_dbs_opportunity) --
        # coverage(alt1)*coverage(alt2) math, attributed per-duplex-group at
        # call time.
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
        # Our DBS78 labels ("AC>CA") already exactly match
        # SigProfilerPlotting's DBS78 reference labels, so no relabeling is
        # needed (unlike ID83). Output filename is controlled by the
        # library: DBS_78_plots_{sample}.pdf.
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

        # Field order/naming mirrors _sbs_burden.txt (uncorrected block,
        # corrected block, per-genome number, opportunity/coverage, ref base
        # number, per-base rate); the 3 "Unmasked burden" lines are omitted
        # since DBS calls have no masked tier to compute an unmasked variant
        # from (_detect_dbs_pairs only ever emits filter=="PASS").
        with open(dbs_dir + "/" + sample + "_dbs_burden.txt", "w") as f:
            f.write(f"Uncorrected burden\t{dbs_burden_uncorrected}\n")
            f.write(f"Uncorrected burden 95% lower\t{dbs_burden_uncorrected_lb}\n")
            f.write(f"Uncorrected burden 95% upper\t{dbs_burden_uncorrected_ub}\n")
            f.write(f"Uncorrected mutation number\t{int(dbs_mut.sum())}\n")
            f.write(f"Corrected burden\t{dbs_burden_corrected}\n")
            f.write(f"Corrected burden 95% lower\t{dbs_burden_corrected_lb}\n")
            f.write(f"Corrected burden 95% upper\t{dbs_burden_corrected_ub}\n")
            f.write(f"Corrected mutation number\t{dbs_mut_corrected.sum()}\n")
            # dbs_mut_corrected.sum() is already the genome-wide
            # extrapolated DBS count (correction_ratio == ref/opp per
            # channel here, with no MH-style grouping the way indel83
            # has, so it's algebraically identical to a rate*ref "hap"
            # sum -- no separate hap array needed) -- hence the same value
            # is used for both "Corrected mutation number" and "Mutation
            # number per genome" below.
            f.write(f"Mutation number per genome\t{dbs_mut_corrected.sum()}\n")
            # Effective coverage at min_group_size=1 (dbs_cov_by_minread[0],
            # matching _dbs_burden_by_group_size.txt's min_group_size=1 row)
            # rather than the reference genome's raw dinucleotide total
            # (Reference base number below).
            f.write(f"Genome coverage\t{dbs_cov_by_minread[0]}\n")
            f.write(f"Reference base number\t{reference_base_number}\n")

        # Process unique mutations to extract duplex depths and create allele counts table
        print("......Processing unique mutations for duplex allele counts.......")
        coverage_file = f"{prefix}/{sample}_coverage.bed.gz"
        allele_counts_data = []

        # Open the coverage file using TabixFile
        tbx = TabixFile(coverage_file)

        n_unique_muts = len(unique_mutations)
        allele_counts_progress = {"start": time.time(), "last": time.time()}
        for i, ((chrom, pos, ref, alt), counts) in enumerate(unique_mutations.items()):
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
                # Columns 3-6 are per-alt-base L-weighted coverage (floats);
                # columns 7-20 are power-weighted indel coverage per category
                # (INDEL_COVERAGE_CATEGORY_LABELS, funcs/misc.py), analogous
                # to columns 3-6 but for indel calling. Kept as a float
                # (rounded to 3dp, not int-truncated) since this is a
                # power-weighted effective depth, not a raw read count --
                # truncating would silently zero out (and drop, see the
                # duplex_depth==0 check below) any sub-1.0 but real
                # single-family site, and bias duplex_vaf's denominator low.
                alt_col = {"A": 3, "T": 4, "C": 5, "G": 6}
                for row in tbx.fetch(chrom, pos - 1, pos):
                    parts = row.split("\t")
                    if len(parts) >= 21:
                        if len(ref) == 1 and len(alt) == 1 and alt in alt_col:
                            # SNV: use the coverage column for this mutation's alt base
                            duplex_depth = round(float(parts[alt_col[alt]]), 3)
                        else:
                            # INDEL: use the single coverage column matching
                            # this specific event's own category (element 3
                            # of counts, stashed when this mutation_key was
                            # first seen above) -- not a sum across all 14,
                            # which would mix in unrelated categories' coverage.
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

        # Create DataFrame and write to file
        allele_counts_df = pd.DataFrame(allele_counts_data)
        allele_counts_file = args.prefix + "/" + sample + "_duplex_allele_counts.txt"
        allele_counts_df.to_csv(allele_counts_file, sep="\t", index=False)
        print(f"Duplex allele counts written to: {allele_counts_file}")
        if args.genebed:
            print(f"Calculating per-gene duplex coverage")
            gene_dict = dict()
            cds_len = dict()
            for rec in TabixFile(args.genebed).fetch():
                (chrom, start, end, gene_exon) = rec.split("\t")[0:4]
                start = int(start)
                end = int(end)
                gene = gene_exon.split("_")[0]
                gene_exon_cov = 0
                for loc in TabixFile(f"{prefix}/{sample}_coverage.bed.gz").fetch(
                    chrom, start, end
                ):
                    gene_exon_cov += int(loc.split("\t")[3])
                if gene_dict.get(gene):
                    gene_dict[gene] += gene_exon_cov
                    cds_len[gene] += end - start + 1
                else:
                    gene_dict[gene] = gene_exon_cov
                    cds_len[gene] = end - start + 1
            with open(args.prefix + "/" + sample + "_gene_coverage.txt", "w") as f:
                for gene in gene_dict.keys():
                    f.write(f"{gene}\t{gene_dict[gene]/cds_len[gene]}\n")

    else:
        print("Re-estimating mutational burden from provided bed file")
        ref_trinuc = calculate_ref_trinuc(args)
        prefix = args.prefix
        if len(prefix.split("/")[-1]) == 0:
            sample = prefix.split("/")[-2]
        else:
            sample = prefix.split("/")[-1]
        sbs_dir, indel_dir, dbs_dir = _ensure_type_subdirs(prefix)
        tn_int = h5py.File(args.reference + ".tn.h5", "r")
        # trinuc_sbs_cov_64[t, b]: L-weighted coverage for raw (un-folded)
        # trinuc context t (0-63), alt base b. Reverse-complement pairs are
        # only combined into the 96 canonical classes at the final burden
        # computation below, via combine_raw192_to_sbs96.
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
        indel_cov_total = 0
        revcomp = _REVCOMP
        vcf = VCF(sbs_dir + "/" + sample + "_sbs.vcf", "r")
        vcf_indel = VCF(indel_dir + "/" + sample + "_indel.vcf", "r")
        indel_count = 0
        for interval in TabixFile(args.reestimatebed, parser=pysam.asBed()).fetch():
            for loc in TabixFile(f"{prefix}/{sample}_coverage.bed.gz").fetch(
                interval.contig, interval.start, interval.end
            ):
                # Coverage bed format: chrom start end A_cov T_cov C_cov G_cov
                # del_unit ins_unit del_2 del_3 del_4 del_5plus
                # ins_A ins_T ins_C ins_G ins_len2 ins_len3 ins_len4 ins_len5plus
                parts = loc.split("\t")
                a_cov = float(parts[3])
                t_cov = float(parts[4])
                c_cov = float(parts[5])
                g_cov = float(parts[6])
                indel_cov_total += int(sum(float(v) for v in parts[7:21]))
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
                if args.dilute:
                    vcf_out.write(rec)
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

        # Select the 192 raw (un-folded) cells, then combine reverse-
        # complement pairs into the 96 canonical classes here, at final
        # estimation time.
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
        burden_uncorrected_lb, burden_uncorrected_ub = poisson_confint(mutnum, cov)
        trinuc_rate = _safe_rate(trinuc_mut, trinuc_cov_96)
        # Per-SBS96-class correction ratio — no np.repeat needed
        correction_ratio = _safe_correction_ratio(ref_trinuc_96, trinuc_cov_96)
        mutnum_corrected = correction_ratio * trinuc_mut
        burden_corrected = mutnum_corrected.sum() / cov
        burden_corrected_lb, burden_corrected_ub = poisson_confint(
            mutnum_corrected.sum(), cov
        )
        hap_trinuc = trinuc_rate * ref_trinuc_96
        mut_per_genome = hap_trinuc.sum()

        with open(sbs_dir + "/" + sample + "_sbs_burden_re_estimate.txt", "w") as f:
            f.write(f"Uncorrected burden\t{burden_uncorrected}\n")
            f.write(f"Uncorrected burden 95% lower\t{burden_uncorrected_lb}\n")
            f.write(f"Uncorrected burden 95% upper\t{burden_uncorrected_ub}\n")
            f.write(f"Uncorrected mutation number\t{int(mutnum)}\n")
            f.write(f"Corrected burden\t{burden_corrected}\n")
            f.write(f"Corrected burden 95% lower\t{burden_corrected_lb}\n")
            f.write(f"Corrected burden 95% upper\t{burden_corrected_ub}\n")
            f.write(f"Corrected mutation number\t{mutnum_corrected.sum()}\n")
            f.write(f"mutation number per genome\t{mut_per_genome}\n")
            f.write(f"genome coverage\t{genome_cov}\n")
        # indel_cov_total sums the same 14 overlapping opportunity columns
        # (parts[7:21] above) as INDEL_COVERAGE_CATEGORY_LABELS elsewhere in
        # this file -- mutations/indel_cov_total is therefore per
        # opportunity-COLUMN, not per locus, exactly the same units mismatch
        # _indel_burden.txt's indel_locus_multiplier corrects for. This
        # branch has no access to that multiplier (it's computed in the
        # `if not args.reestimatebed` branch above, using ref_indel83 from
        # the main indel83 pipeline, neither of which runs here), so
        # recompute the reference-genome side of it directly from the index
        # files the same way the main branch does.
        label2num_100_re, _ = build_indel100_labels()
        ref_indel100_re = calculate_ref_indel100(args)
        ref_indel76_re = combine_indel100_to_indel76(ref_indel100_re, label2num_100_re)
        ref_indel83_re = expand_indel76_to_indel83(ref_indel76_re)
        ref_indel83_re = override_inshp0_with_next_base_opportunity(
            ref_indel83_re, ref_indel100_re
        )
        indel_locus_multiplier_re = ref_indel83_re.sum() / ref_trinuc.sum()
        indel_burden_re = (
            indel_count / float(indel_cov_total) * indel_locus_multiplier_re
            if indel_cov_total > 0
            else float("nan")
        )
        indel_burden_re_lb, indel_burden_re_ub = poisson_confint(
            indel_count, indel_cov_total
        )
        indel_burden_re_lb *= indel_locus_multiplier_re
        indel_burden_re_ub *= indel_locus_multiplier_re
        with open(indel_dir + "/" + sample + "_indel_burden_re_estimate.txt", "w") as f:
            f.write(f"Indel burden\t{indel_burden_re}\n")
            f.write(f"Indel burden 95% lower\t{indel_burden_re_lb}\n")
            f.write(f"Indel burden 95% upper\t{indel_burden_re_ub}\n")
            f.write(f"Indel number\t{indel_count}\n")
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
