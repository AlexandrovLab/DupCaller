import math
import os
import re
import time
from multiprocessing import Pool
import numpy as np
import h5py
import pysam
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
import pandas as pd
from scipy.optimize import brentq


from .prob import calculateSSPosterior, calculateDSPosterior, indelErrorProbs


# ======================================================================
# Generic / file / progress utilities
# ======================================================================
def _ensure_type_subdirs(prefix):
    """Per-mutation-type output subfolders under a sample's prefix
    directory — SBS/, INDEL/, DBS/ — created (if missing). Every
    SBS/indel/DBS-specific output (VCFs, corrected tables, burden files,
    plots, by-read-group breakdowns, and each type's own
    _*_by_duplex_group.txt coverage file) lives under the matching
    subfolder; only genuinely whole-sample, type-agnostic files
    (_stats.txt, _coverage.bed.gz, _duplex_allele_counts.txt, learned
    error matrices under ERROR/, etc.) stay at the top level of prefix.
    Shared by Caller.py (writer of the VCFs) and Estimate.py (writer of
    everything else, and reader of the VCFs) so both agree on the exact
    same subfolder layout.

    Returns (sbs_dir, indel_dir, dbs_dir).
    """
    sbs_dir = os.path.join(prefix, "SBS")
    indel_dir = os.path.join(prefix, "INDEL")
    dbs_dir = os.path.join(prefix, "DBS")
    for d in (sbs_dir, indel_dir, dbs_dir):
        os.makedirs(d, exist_ok=True)
    return sbs_dir, indel_dir, dbs_dir


def log_progress(state, label, index, total=None, extra="", interval=5.0):
    """Print a throttled "clock" progress line — elapsed time, item
    count, and percent-complete when the total is known — while iterating
    a long-running loop (mutation filtering, VCF streaming, per-chromosome
    reference composition, etc). Printed at most once every `interval`
    seconds (always on the final item when `total` is known, so loops
    that finish quickly still get a 100% line), so it doesn't spam logs
    for large inputs.

    state: a dict the caller keeps around across the whole loop —
        {"start": <loop start time>, "last": <time of last print>}.
    total: item count if known upfront (e.g. len(some_list)), enabling a
        "i/total (pct%)" display. Pass None for loops with no cheap total
        (e.g. streaming a VCF record-by-record) to fall back to a plain
        running count instead; pass 0 to no-op (nothing to iterate).
    extra: optional short string identifying the current item (e.g.
        "chr1:12345" or a chromosome name), appended to the line.
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


# Labels for the 14 power-weighted indel coverage categories, in the same
# order as the trailing 14 columns of the coverage bed file (after the 4
# SBS columns) and the "Effective/Unmasked Indel Coverage (<label>)" lines
# in _stats.txt. Shared between Caller.py (writer) and Estimate.py/
# Summarize.py (readers) so the label strings can't drift out of sync.
#
# Deletion Length 2-5+ and Insertion Length 2-5+ are always computed at a
# fixed homopolymer-length-1 context (microhomology-style for deletion;
# "novel repeat unit" / rep0 for insertion, mirroring
# classify_indel_record's is_hp/indel_len>=2 branch), not the position's
# actual repeat context — Insertion Length columns are further restricted
# to is_hp positions, since only there can no-repeat-yet even apply.
#
# Insertion A/T/C/G are the per-inserted-base rep0 opportunity for a 1bp
# insertion: each uses the L-table power for inserting *that* base
# specifically (fixed hps=1 context) at every position, except positions
# where the immediately following reference base already equals that
# base — there, the value is 0, since inserting that base would extend an
# existing run rather than start a new one. This has nothing to do with
# the position's own reference base (see funcs/misc.py's
# indel100_reference_bucket_indices / classify_indel_channel).
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
    # Columns 0/1 ("...of Repeat Unit") are STR-only now -- opportunity to
    # delete/insert exactly one copy of THIS position's own annotated STR
    # unit, 0 wherever no real STR is annotated (unit_len<2). The
    # 1bp-homopolymer-slip opportunity, which exists everywhere regardless
    # of STR status (genotypeDSIndel always scores a +-1bp event via the
    # HP table, never the STR one -- see funcs/prob.py's indelErrorProbs),
    # lives entirely here instead: unconditional L_indel_1bp coverage
    # using the position's true self-derived homopolymer run (hp_raw_np).
    # Splitting it out this way -- rather than merging it into columns
    # 0/1 via is_hp, as those used to do -- keeps every column
    # non-overlapping, so summing across all 16 (e.g. indel_cov_total)
    # gives a real total instead of double-counting the 1bp opportunity
    # at every non-STR position. Matches
    # indel100_reference_bucket_indices' own convention (funcs/misc.py):
    # unconditional HP-bucket credit, real_str-gated STR-bucket credit,
    # no other overlap. See call.py's cov_mat_indel[:, 14]/[:, 15].
    "Deletion of Homopolymer",
    "Insertion of Homopolymer",
]


def parse_stats_file(path):
    """Parse a DupCaller _stats.txt file into a {label: value_string} dict,
    keyed by the text before the first tab on each line. Label-based lookup
    avoids the fragility of reading by fixed line number, which breaks
    silently if a line is ever added/removed/reordered upstream."""
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
# in a way that ndim/shape checks alone can't detect — e.g. version 1's
# hp.h5 was a single 3-row (unit_len, is_start, repeat_count) file built
# from --repeatBed/--strbed for both homopolymers and STRs; version 2 (this
# one) splits that into hp.h5 (2 rows: hp, cut — self-derived directly from
# the reference sequence, homopolymers only) and a new str.h5 (3 rows:
# unit_len, hp, cut — still built from --repeatBed/--strbed, restricted to
# unit_len>=2). See Index.py's do_index and funcs/misc.py's
# load_repeat_context. `DupCaller.py index` (Index.py) stamps this onto
# each h5 file's root attrs; check_h5_usable below verifies it.
INDEX_SCHEMA_VERSION = 2


def check_h5_usable(h5_path, expected_ndim=None):
    """Try to open and read from an h5 file; return (ok, error_message).

    Checks the index_schema_version file attribute first: its absence means
    the file predates version stamping entirely (an outdated index built by
    an old `DupCaller.py index`), and a mismatched value means it was built
    by a version of DupCaller with a different on-disk layout — neither of
    which the ndim/shape checks below can catch on their own.
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
    for indices 0-63. Index i (i<32) and i+32 are always exact reverse
    complements of each other, by construction of these two loops (verified
    by direct execution, not incidental). Returns (trinuc2num, num2trinuc).
    Shared by call.py, Caller.py, and Estimate.py so this ordering can't
    drift out of sync between the writer and readers of trinuc-indexed data.
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
    formatted "{trinuc}>{alt}". No reverse-complement folding — a context's
    own reference base is skipped entirely rather than included as a
    self-mutation class. Returns (label2num, labels); this order is the
    contract between Caller.py (writer of _trinuc_by_duplex_group.txt) and
    Estimate.py (reader) and must not be built independently in either.
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
    canonical class with its reverse-complement partner from the other half
    of the 64-context space. This is the single point where reverse-
    complement pairs get folded together, deferred to estimation time
    rather than done at raw accumulation time.

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
    """np.where(cov > 0, mut / cov, 0), without the RuntimeWarning:
    np.where evaluates BOTH branches in full before selecting, so
    `mut / cov` still runs — and warns — at every position where cov==0,
    even though the result there is thrown away in favor of the literal
    0. The final values are identical either way; this just silences the
    warning for a divide-by-zero that was always going to be discarded."""
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
    """The 9 possible alt dinucleotides for a reference dinucleotide: each
    of the 2 positions independently becomes one of its 3 non-reference
    bases (a DBS event changes both positions simultaneously — if only
    one changed it would be a plain SBS, not a DBS)."""
    return [
        b1 + b2
        for b1 in _DBS_BASES
        if b1 != ref[0]
        for b2 in _DBS_BASES
        if b2 != ref[1]
    ]


def build_dbs_raw144_labels():
    """The 144 raw (un-folded) "REF>ALT" DBS labels — all 16 raw
    reference dinucleotides x their 9 alts each — the same role for DBS
    that build_trinuc192_labels plays for SBS: a label space wide enough
    to hold any observed/raw event before folding down to the 78 canonical
    DBS78 channels (Estimate.py's combine_raw_dbs_to_dbs78). Shared by
    call.py/Caller.py (writer of _dbs_by_duplex_group.txt) and Estimate.py
    (reader, and genome-wide reference-composition folding) so this order
    can't drift out of sync between them — same contract as
    build_trinuc192_labels/build_indel100_labels.
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
    """[4,4,4,4] int lookup grid mapping (ref1_num, ref2_num, alt1_num,
    alt2_num) -- all in "ATCG"-order 0-3 encoding, the same convention
    call.py's ref_int/alt_int arrays already use (funcs/call.py:64,459) --
    directly to a row index into the 144-raw-class DBS opportunity space
    (build_dbs_raw144_labels), for cheap vectorized accumulation in
    call.py's per-family hot loop without building label strings per
    position pair. Entries where alt1==ref1 or alt2==ref2 (not a valid DBS
    event -- both positions must actually change) are left at -1, a
    sentinel the caller filters out before np.bincount.
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


# [4,4,4,4] (ref1_num, ref2_num, alt1_num, alt2_num) -> row index into
# the 144-raw-class DBS opportunity space -- built once at import time
# (cheap: 144 valid entries out of 256) rather than per-call, and reused
# by every duplex-group DBS-opportunity accumulation below.
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
):
    """Detect DBS (dinucleotide substitution) events within one duplex
    read group's SNV candidate list: two positions changing
    SIMULTANEOUSLY in the SAME physical molecule, not two independent
    SNVs that happen to sit next to each other in the genome. Because
    mut_positions/muts_ind/ref_int/alt_int all come from THIS one read
    group's own consensus calling, a genomically-adjacent pair found here
    is, by construction, supported by the same underlying reads — unlike
    scanning the already-written SNV VCF after the fact, which has no way
    to tell same-molecule adjacency from coincidence.

    Only pairs where BOTH constituent SNVs individually reach full PASS
    status are called as a DBS (recomputed here from pass_bool/
    LR_pass_bool/flt_rs — the same three inputs the existing per-position
    loop above already uses to decide flt — rather than threading a new
    list through that loop, so this can be added without touching any of
    its existing lines).

    Returns a list of DBS mut dicts (ref/alt are 2-character
    dinucleotides), in the same shape as the SNV mut dict built above,
    minus the SNV-only INFO fields (LR/LM/TC/BC/TN/HP/STR) that don't
    apply to a 2-base event.
    """

    def _flt(i):
        if pass_bool[i] and LR_pass_bool[i]:
            return flt_rs
        return "masked"

    dbs_muts = []
    for k in range(len(mut_positions) - 1):
        if mut_positions[k + 1] != mut_positions[k] + 1:
            continue
        i0, i1 = muts_ind[k], muts_ind[k + 1]
        if _flt(i0) != "PASS" or _flt(i1) != "PASS":
            continue
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
                },
                "formats": ["AC", "RC", "DP"],
                # No independent raw-BAM depth re-verification for DBS
                # yet (mirroring extractDepthSnv/extractDepthBatchSnv for
                # SNVs) — both constituent bases already passed full
                # duplex-consensus LR calling, which is a strong
                # signal on its own, but there's currently no equivalent
                # of the SNV/indel maxAF / normal-VAF post-processing
                # check for DBS. Placeholder AC=1/RC=0/DP=1 (tumor),
                # 0/0/0 (normal) so downstream VCF writing has values to
                # write; a real depth-extraction pass would replace this.
                "samples": [[1, 0, 1], [0, 0, 0]],
            }
        )
    return dbs_muts


def _compute_dbs_opportunity(cov_mat, ref_int, pass_bool):
    """Per-duplex-group DBS opportunity contribution from one family's
    window: for every pair of directly ADJACENT in-window positions
    (i, i+1) that both individually pass calling (pass_bool[i] and
    pass_bool[i+1]) and have a valid A/T/C/G reference base at both,
    opportunity is coverage(alt1 at i) * coverage(alt2 at i+1) summed over
    all 9 non-ref (alt1, alt2) combinations -- the same product-of-
    coverages approximation Estimate.py's calculate_dbs_opportunity uses
    genome-wide from the coverage bed, just computed here per-family (so
    it can be attributed to that family's own duplex_no) from the same
    cov_mat/ref_int/pass_bool the trinuc accumulator right after this
    function's call site already uses -- window-local arrays, aligned the
    same way.

    cov_mat  : (window_len, 4) L-weighted per-alt-base coverage.
    ref_int  : (window_len,) reference base, 0-3=A/T/C/G, 4=N/invalid.
    pass_bool: (window_len,) bool, calling-eligible positions.
    Returns a (144,) raw-DBS-class contribution vector (funcs/misc.py's
    build_dbs_raw144_labels order), zero wherever no adjacent pass_bool
    pair with valid reference bases exists.
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
# Deletion length/count bins start at 1 (you need an existing run of at
# least 1 to delete from); insertion STR/count bins start at 0 (you can
# insert into a context with no existing repeat at all). A repeat_count of
# 0 or 1 is "not really a repeat" (see funcs/call.py's classification
# code): for deletion both fold into bucket "1"; for STR insertion both
# are counted under *both* bucket "0" and bucket "1".
INDEL_DEL_HP_LEN_BINS = ["1", "2", "3", "4", "5", "6+"]
# Homopolymer-insertion length bins at ID83 (output) resolution — includes
# "0" (rep0: inserting into a context with no existing repeat), since
# Cinshp0/Tinshp0 are real COSMIC ID83 channels. Used only by
# build_indel83_labels; NOT by the raw100/indel76 pipeline below (see
# INDEL100_INS_HP_LEN_BINS) since the hp.h5-annotation-based rep0
# opportunity computed at that resolution is always discarded and
# replaced by the exact next-base computation in
# Estimate.py's override_inshp0_with_next_base_opportunity.
INDEL_INS_HP_LEN_BINS = ["0", "1", "2", "3", "4", "5+"]
# Same, but for the raw100/indel76 pipeline, which never computes a real
# rep0 opportunity for homopolymer insertion at all (see above) — so
# there's no "0" bin to carry through here.
INDEL100_INS_HP_LEN_BINS = ["1", "2", "3", "4", "5+"]
INDEL_STR_UNIT_BINS = ["2", "3", "4", "5+"]
INDEL_DEL_STR_COUNT_BINS = ["1", "2", "3", "4", "5", "6+"]
INDEL_INS_STR_COUNT_BINS = ["0", "1", "2", "3", "4", "5+"]
INDEL_MH_DEL_LEN_BINS = ["2", "3", "4", "5+"]


def build_indel100_labels():
    """100 raw indel classification labels (fine resolution; the "100"
    name predates the 4 1bpins{base} columns at the end, mirroring how
    "indel76"/"indel83" elsewhere are also historical rather than renamed
    on every change): base-specific homopolymer deletion (4 bases x 6
    length bins = 24), base-specific homopolymer insertion (4 bases x 5
    length bins = 20 — no "0"/rep0 bin here, see INDEL100_INS_HP_LEN_BINS;
    that opportunity is always thrown away and replaced anyway, computed
    instead from the 1bpins{base} columns below), unit-size/repeat-count
    specific STR deletion/insertion (4 unit-size buckets x 6 count bins =
    24 each), microhomology deletion length (4 bins), and 4 additional
    fixed-context "rep0 opportunity for inserting exactly this base"
    columns (see indel100_reference_bucket_indices and funcs/call.py's
    matching per-family accumulation) — used to compute the rep0 ("no
    repeat existed yet") bin of the ID83 Cinshp/Tinshp channels correctly
    (Estimate.py's override_inshp0_with_next_base_opportunity). Each
    1bpins{base} column represents the opportunity/power to insert *that
    specific base* with no pre-existing repeat: every position counts,
    EXCEPT where the immediately following reference base already equals
    that base (there it's 0/excluded, since inserting it would extend an
    existing run rather than start a new one) — independent of the
    position's own reference base or any repeat annotation.
    24+20+24+24+4+4 = 100. Deterministic order, shared by the writer
    (Caller.py) and any future reader — must not be built independently
    elsewhere.
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
    """Genome-wide "opportunity" bucket contributions for the indel100
    scheme (build_indel100_labels above), given whole-chromosome repeat-
    context arrays read directly from the .hp.h5/.str.h5/.ref.h5 index
    files: hp_run_arr = hp.h5 row 0 (self-derived homopolymer run length,
    always valid, unconditional on any STR annotation), str_unit_len_arr/
    str_repeat_count_arr = str.h5 rows 0/1 (BED-derived STR annotation,
    unit_len==0 wherever no STR is annotated), ref_base_arr = ref.h5 (all
    base2num-encoded: A=0,T=1,C=2,G=3). Every position contributes weight
    1 to each column it's eligible for.

    hp and STR opportunity are computed independently here, deliberately
    NOT merged with STR taking priority (unlike load_repeat_context, used
    elsewhere for indel *calling*'s classification, where a position can
    only be one category at a time): a position inside an annotated STR
    interval that also starts an embedded run of identical bases (e.g.
    the "AA" in a (AAT)n repeat) contributes to BOTH the homopolymer and
    the STR opportunity columns.

    Real STR opportunity (an annotated tract's own unit_len/repeat_count)
    is read straight off str.h5, cut-gated, no re-scanning. Positions
    str.h5 doesn't annotate additionally get a flat, uniform rep1
    (deletion) / rep0+rep1 (insertion) credit for every unit size — every
    such position is a candidate for an arbitrary, non-repeating U-bp
    indel, but discovering whether a *longer* run coincidentally exists
    there would require rescanning the reference, which this deliberately
    doesn't do. That flat credit is a straight per-position count here
    (this function has no read-depth/family concept), and call.py's
    matching per-family credit is the analogous fixed-context L-table
    lookup (see cov_mat_indel's del_2..del_5plus / "Insertion Length
    2..5+" columns) — both scan-free, so this side and the coverage side
    agree by construction rather than by two independent approximations
    happening to line up.

    hp_cut_arr, str_cut_arr (optional): hp.h5 row 1 / str.h5 row 2
    (start-of-run booleans). A real, annotated repeat tract is one
    molecular locus, not one-opportunity-per-base-it-spans — a slippage
    event there is always reported at a single canonical (left-aligned)
    anchor, no matter how long the tract is. Crediting every base of a
    length-L tract as an independent opportunity (the old, unconditional
    behavior) makes a bucket's aggregate opportunity scale with L, and
    since coverage/mappability degrades with L too (see funcs/call.py's
    nm_mask), that L-weighting amplifies exactly the loci most likely to
    be under-covered — inflating correction_ratio well beyond what the
    coverage gap alone would cause. Passing these arrays restricts
    homopolymer/STR opportunity to cut==1 positions only, one credit per
    tract regardless of its length. If omitted, every position is
    credited as before (backward compatible with older callers).

    next_ref_base_arr (optional): ref_base_arr shifted by one position —
    next_ref_base_arr[i] is the reference base immediately following
    position i (base2num-encoded; use an out-of-range sentinel, e.g. -1,
    for positions with no valid next base, such as a chromosome's last
    base). Feeds the 4 "1bpins{base}" columns (96-99, base2num order
    A/T/C/G): column 96+b counts every position where the *following*
    base differs from b — the rep0 opportunity to insert base b there —
    independent of unit_len/repeat_count or the position's own base. If
    omitted, those 4 columns get no contribution (backward compatible
    with older callers).

    This mirrors — but does not share code with — the L-table-weighted
    per-position classification in funcs/call.py's per-window indel
    coverage loop (search that file for "Fine-grained 100-class raw indel
    classification"); the two independently derive the same anchor
    convention (a position P is "opportunity for a 1bp/unit event whose
    VCF-reported anchor would be P-1"), which is also why
    Estimate.py.classify_indel_record below looks up context at rec.pos
    directly (VCF's 1-based anchor position numerically equals P in
    0-based terms). Keep bucket boundaries here in sync with both if
    build_indel100_labels ever changes.

    Returns the (100,) accumulated opportunity-count array directly: each
    contribution is folded straight into the running (100,) total as it's
    computed -- a scalar increment for columns shared by every masked
    position, or a small per-call np.bincount for columns that vary
    per-position. Every contribution has weight exactly 1, so this avoids
    keeping up to ~20 separate whole-chromosome-length index/weight array
    pairs alive simultaneously for one final concatenate, which for a real
    chromosome would be the dominant memory cost of computing
    reference-genome indel composition (see calculate_ref_indel100).
    """
    n = hp_run_arr.shape[0]
    ref_base_safe = np.where((ref_base_arr < 0) | (ref_base_arr > 3), 0, ref_base_arr)
    real_str = str_unit_len_arr >= 2
    base_cgta = _INDEL_BASE2NUM_TO_CGTA[ref_base_safe]
    hp_cut = (
        np.asarray(hp_cut_arr, dtype=bool)
        if hp_cut_arr is not None
        else np.ones(n, dtype=bool)
    )
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

    # 1-4: homopolymer deletion opportunity (cols 0-23) — one credit per
    # tract (hp_cut==1), not one per base spanned. hp_run==1 ("delrep1",
    # an isolated, non-repeated base) is unaffected either way: a
    # length-1 run's single position is always its own cut by
    # definition, so gating never excludes it.
    del_bucket = np.clip(hp_run_arr, 1, 6) - 1
    add_varying(base_cgta * 6 + del_bucket, hp_cut)

    # 5-8: homopolymer insertion opportunity (cols 24-43) — length bins
    # 1..5+ only, one credit per tract (hp_cut==1). repeat_count 0 (rep0:
    # no existing repeat at all) is deliberately NOT credited here: the
    # hp.h5-annotation-based opportunity for that bin is always thrown
    # away and replaced by the exact next-base computation in the
    # 1bpins{base} columns below (see Estimate.py's
    # override_inshp0_with_next_base_opportunity), so there's no point
    # computing it here at all.
    ins_bucket = np.clip(hp_run_arr, 1, 5) - 1
    add_varying(24 + base_cgta * 5 + ins_bucket, (hp_run_arr >= 1) & hp_cut)

    # 9: STR deletion opportunity (cols 44-67) — real part: gated by
    # real_str (str_unit_len_arr>=2) AND str_cut==1, one credit per
    # annotated tract, not one per base it spans.
    unit_bucket = np.clip(str_unit_len_arr, 2, 5) - 2
    count_bucket_del = np.clip(str_repeat_count_arr, 1, 6) - 1
    add_varying(44 + unit_bucket * 6 + count_bucket_del, real_str & str_cut)

    # Flat rep1 opportunity: every position is a candidate for an
    # arbitrary, non-repeating U-bp deletion (repeat_count=1) for each
    # unit size 2/3/4/5 — a fixed, uniform credit, NOT an attempt to
    # discover a true higher count that might coincidentally exist there
    # (that would need rescanning the reference, which is exactly what
    # call.py's matching per-family credit deliberately doesn't do either
    # — see its own comment). Kept deterministic and scan-free so this
    # side and the coverage side (call.py's cov_mat_indel/indel100
    # accumulation) agree by construction rather than by two
    # independently-approximated scans happening to line up.
    #
    # Deliberately NOT gated by ~real_str (unlike an earlier version of
    # this function): a hypothetical, non-matching U-bp deletion is a
    # real, independent opportunity even at a position that also has a
    # real, different-unit STR annotation -- e.g. inside ATGATGATG (a
    # real unit=3 tract), an arbitrary 2bp deletion that doesn't match
    # the ATG unit is still a genuine "delSTR2_rep1"-type event. real_str
    # positions already get their own unit's real credit above; this
    # flat credit is for every *other* hypothetical unit size,
    # independently, at every position (including real_str ones).
    for u_idx in range(4):
        out[44 + u_idx * 6] += n

    # 10: STR insertion opportunity (cols 68-91) — real part: one credit
    # per annotated tract (real_str & str_cut), with the same rep-count
    # 0/1 ambiguity hedge as the homopolymer deletion case above (STR
    # insertion still has a real "0" bin, unlike homopolymer insertion).
    real_str_cut = real_str & str_cut
    count_bucket_ins = np.clip(str_repeat_count_arr, 0, 5)
    add_varying(68 + unit_bucket * 6 + count_bucket_ins, real_str_cut)
    is_rep1 = str_repeat_count_arr == 1
    add_varying(68 + unit_bucket * 6 + 0, real_str_cut & is_rep1)

    # Flat rep0/rep1 opportunity ("what if a U-bp unit were inserted
    # here"): same uniform, scan-free credit as the deletion case above,
    # into both the "0" (novel, no repeat yet) and "1" (first extension)
    # buckets for every unit size, at every position -- deliberately not
    # gated by ~real_str, same reasoning as the deletion case above.
    for u_idx in range(4):
        base_col = 68 + u_idx * 6
        out[base_col] += n
        out[base_col + 1] += n

    # 11: microhomology deletion opportunity (cols 92-95) — every position,
    # fixed context, independent of actual repeat annotation.
    for k in range(4):
        out[92 + k] += n

    # 12: 1bp-insertion "next base" opportunity (cols 96-99) — a position
    # is opportunity for "insert base N here, rep0" purely if its own
    # following base isn't N; independent of unit_len/repeat_count/the
    # position's own base entirely, per classify_indel_channel's rep0
    # definition (funcs/misc.py).
    if next_ref_base_arr is not None:
        next_base_valid = (next_ref_base_arr >= 0) & (next_ref_base_arr <= 3)
        for b in range(4):
            add_uniform(96 + b, next_base_valid & (next_ref_base_arr != b))

    return out


def load_repeat_context(hp_dataset, str_dataset, start=None, end=None):
    """Merge the self-derived homopolymer index (hp.h5: 2 rows, hp/cut —
    see Index.py's do_index) and the BED-derived STR repeat index (str.h5:
    3 rows, unit_len/hp/cut, unit_len>=2 only) for one chromosome (or a
    start:end slice of it) into the single 3-row (unit_len, cut,
    repeat_count) layout every consumer already expects from the old
    combined hp.h5 (unit_len=row0, is_start/cut=row1, repeat_count=row2) —
    so nothing downstream of this needs to know the index is now split
    across two files.

    Wherever a position falls inside an annotated STR interval (str
    unit_len >= 2), the STR annotation wins; everywhere else this falls
    back to the self-derived homopolymer run (unit_len implicitly 1) --
    safe because an annotated STR interval always describes an actual
    repeat, with no separate "background" source to conflict with.

    hp_dataset, str_dataset: h5py Datasets for one chromosome from hp.h5
        and str.h5 respectively (e.g. hp_h5[chrom], str_h5[chrom]).
    start, end: optional slice bounds; defaults to the whole chromosome.

    Returns a (3, L) int16 array: row0=unit_len, row1=cut, row2=repeat_count.
    Values never exceed 127 (Index.py clips them before writing), so int16
    is already generous — it's only used instead of int8/uint8 to leave
    headroom for downstream arithmetic (e.g. subtracting a small constant)
    without any risk of silent unsigned wraparound, while avoiding the 8x
    widening a plain `int` (numpy int64) would cost on the on-disk uint8
    data -- several GB of unnecessary temporary allocation for a whole
    chromosome (the common case for Estimate.py's genome-composition
    calls, as opposed to call.py's windowed calls).
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


def _prefix_match_len(a, b, max_len):
    """Length of the longest common prefix of `a` and `b`, capped at
    max_len. Matching is inherently monotonic (if the first n characters
    agree, so do the first n-1), so a simple greedy walk from the start
    finds the maximum directly — no need to test every n independently."""
    n = 0
    while n < max_len and a[n] == b[n]:
        n += 1
    return n


def classify_indel_record(indel_seq, ref_after, indel_len):
    """Classify a single OBSERVED indel call by actual sequence content,
    rather than by pre-computed repeat annotation (unit_len/repeat_count
    from .hp.h5) as the previous version of this function did.

    Only ever looks AFTER the event, never before: every indel reaching
    this function has already been left-aligned (funcs/indels.py's
    left_align_indel, applied to every raw CIGAR-derived indel before it
    ever enters a duplex family's consensus set — see genotypeDSIndel/
    profileTriNucMismatches). A left-aligned event's anchor is, by
    construction, the leftmost position that still represents the same
    net change — if a matching copy of the unit existed immediately
    before it, left-alignment would already have shifted the anchor to
    absorb it. So checking "before" would always find zero overlap; it's
    not computed at all rather than computed and asserted to be 0.

    indel_seq: the deleted or inserted bases themselves (uppercase, length
        abs(indel_len)).
    ref_after: reference bases immediately following the indel (only the
        first min(len(ref_after), abs(indel_len)) matter).
    indel_len: signed length of the event (+insertion, -deletion); its
        magnitude must equal len(indel_seq).

    Deletion (indel_len < 0): microhomology length is the longest prefix
        of ref_after matching a prefix of indel_seq (e.g. deleting "ATCG"
        out of "[ATCG]ATCGAAAA" has mh_len=4, since deleted=="ATCG" ==
        the next 4 bases). If mh_len equals the deletion length (the
        entire deleted run duplicates the following sequence), it's a
        repeat-unit ("STR") deletion rather than a microhomology one — a
        partial match is what makes it MH.

    Insertion (indel_len > 0): only an exact, full-length match is checked
        (no partial-overlap microhomology concept for insertions) — if
        indel_seq equals the equal-length window of ref_after, it's
        extending an existing repeat; otherwise it's a novel/rep0
        insertion.

    Returns a dict describing the classification:
        {"category": "str_deletion", "unit_len": L}
        {"category": "mh_deletion", "del_len": L, "mh_len": n}
        {"category": "repeat_insertion", "ins_len": L}
        {"category": "novel_insertion", "ins_len": L}
    """
    assert len(indel_seq) == abs(indel_len)

    if indel_len < 0:
        del_len = -indel_len
        if del_len == 1:
            # A single deleted base is always a homopolymer-run deletion,
            # never microhomology (ID83 has no 1delMH* bucket — MH channels
            # start at del_len==2). This holds regardless of whether the
            # deleted base happens to match a flanking base: a match means
            # a real run of length >=2 (mh_len==del_len below would also
            # catch it), and a non-match just means the run's length is 1
            # (an isolated, non-repeated base) — classify_indel_channel's
            # repeat-count scan already expects and handles that case
            # ("repeat count is 1 ... plus however many more flank it").
            return {"category": "str_deletion", "unit_len": 1}
        mh_len = _prefix_match_len(ref_after, indel_seq, min(len(ref_after), del_len))
        if mh_len == del_len:
            return {"category": "str_deletion", "unit_len": del_len}
        return {"category": "mh_deletion", "del_len": del_len, "mh_len": mh_len}

    ins_len = indel_len
    matches_after = len(ref_after) >= ins_len and ref_after[:ins_len] == indel_seq
    if matches_after:
        return {"category": "repeat_insertion", "ins_len": ins_len}
    return {"category": "novel_insertion", "ins_len": ins_len}


def classify_indel_channel(indel_seq, ref_after, indel_len, anno=None):
    """Resolve a single observed indel call to its exact ID83-style channel
    label (see Estimate.py's build_indel83_labels for the full ordered
    list), building on classify_indel_record's MH-vs-repeat-unit
    determination with an additional repeat-count scan: how many times the
    deleted/inserted unit already repeats in the reference, counting
    outward from the event. That count is what picks the hp1-6+ /
    STR-rep1-6+ (deletion) or hp0-5+ / STR-rep0-5+ (insertion) sub-bucket
    — classify_indel_record alone only determines unit size and
    MH-vs-repeat-unit status, not how long the run actually is. Only
    scans AFTER the event, for the same left-alignment reason
    classify_indel_record only looks after it — see that function's
    docstring.

    For a "str_deletion" (classify_indel_record's mh_len==del_len case),
    the deleted sequence's own length IS the repeat unit — one base for
    homopolymers (unit_len==1, pooled C+G -> "C" / A+T -> "T"), or the
    deletion's own length for STR unit sizes 2/3/4/5+. Repeat count is 1
    (the deleted unit itself) plus however many more of that same unit
    immediately follow it (never "before": that side is always empty for
    a left-aligned event, see classify_indel_record).

    For insertions, "repeat_insertion" (matched after) and
    "novel_insertion" (matched neither) are really the same computation at
    two different repeat counts, not two different code paths: a
    non-repeat-insertion is exactly the repeat_count==0 case, so both
    categories are handled identically here — a 1bp insertion matching
    nothing after it naturally resolves to "{pool}inshp0", never a
    separate "novel" bucket.

    For "mh_deletion", no repeat-count scan applies — the microhomology
    channels are keyed only by (deletion length, microhomology length),
    with the microhomology length floored at 1 and capped at 5 ("5+") per
    the standard ID83 convention (a deletion with zero true overlap still
    falls in the "...delMH1" bucket, since ID83 has no "0" bucket).

    anno: optional (hp_run, str_unit_len, str_repeat_count) read directly
        from hp.h5 (self-derived homopolymer run length) and str.h5
        (unit_len/repeat_count) at this event's anchor position — the same
        raw values the opportunity side (call.py's cov_mat_indel/indel100
        accumulation, funcs/misc.py's indel100_reference_bucket_indices)
        already reads. Purely mechanical, no live re-scanning of the
        reference: a unit_len==1 event (homopolymer) always takes its
        repeat_count from hp_run — hp.h5 is self-derived and always valid,
        so this never needs a fallback. A unit_len>=2 event (STR) takes it
        from str_repeat_count only when str_unit_len matches the event's
        own unit length; otherwise (str.h5 has no matching annotation
        here) repeat_count defaults to 2 — classify_indel_record has
        already confirmed at least one matching copy follows the event
        (that's what makes it "str_deletion"/"repeat_insertion" rather
        than "mh_deletion"/"novel_insertion" in the first place), so 2 is
        the minimum count consistent with what was actually observed,
        without hypothesizing a longer run str.h5 doesn't know about. If
        anno is omitted entirely, the same defaults apply throughout (1
        for homopolymers, 2 for STR).
    """
    result = classify_indel_record(indel_seq, ref_after, indel_len)
    cat = result["category"]

    if cat in ("str_deletion", "repeat_insertion", "novel_insertion"):
        is_ins = cat != "str_deletion"
        unit = indel_seq
        unit_len = len(unit)
        if cat == "novel_insertion":
            # No matching copy exists at all (that's the definition of
            # "novel") -- repeat_count is 0 by construction, not a lookup.
            repeat_count = 0
        elif unit_len == 1:
            repeat_count = anno[0] if anno is not None else 1
        elif anno is not None and anno[1] == unit_len:
            repeat_count = anno[2]
        else:
            repeat_count = 2
        pool = "C" if unit[0] in ("C", "G") else "T"

        if is_ins:
            count_label = "5+" if repeat_count >= 5 else str(repeat_count)
            if unit_len == 1:
                return f"{pool}inshp{count_label}"
            unit_label = "5+" if unit_len >= 5 else str(unit_len)
            return f"{unit_label}insstr{count_label}"
        else:
            count_label = "6+" if repeat_count >= 6 else str(repeat_count)
            if unit_len == 1:
                return f"{pool}delhp{count_label}"
            unit_label = "5+" if unit_len >= 5 else str(unit_len)
            return f"{unit_label}delstr{count_label}"

    # mh_deletion
    del_len = result["del_len"]
    mh_len = result["mh_len"]
    del_group = "5+" if del_len >= 5 else str(del_len)
    max_mh = {"2": 1, "3": 2, "4": 3, "5+": 5}[del_group]
    reported_mh = min(max(mh_len, 1), max_mh)
    mh_label = "5+" if (del_group == "5+" and reported_mh == 5) else str(reported_mh)
    return f"{del_group}delMH{mh_label}"


def indel_coverage_category_index(indel_seq, ref_after, indel_len):
    """Map a single observed indel event to its one matching column index
    (0-15) into coverage.bed.gz's 16 power-weighted indel opportunity
    columns, i.e. INDEL_COVERAGE_CATEGORY_LABELS above -- the exact
    channel this specific mutation's own detection power was screened/
    thresholded against (call.py's cov_mat_indel columns, same order),
    not an average across all 16 unrelated categories.

    Reuses classify_indel_record's repeat-vs-mismatch determination for
    the 1bp cases (the only place this coarser scheme distinguishes
    repeat-unit from novel/mismatched events -- deletions and insertions
    of length >=2 are bucketed by length alone here, regardless of
    true-STR vs microhomology, since that finer distinction only exists
    in the 83-channel scheme classify_indel_channel builds on top of
    this).

    Every +-1bp repeat-matching event lands on the homopolymer columns
    (14/15), never 0/1 (STR-only, one copy of a real annotated repeat
    unit -- a multi-bp concept a 1bp event never satisfies): matches
    genotypeDSIndel, which always scores a +-1bp event via the HP table
    regardless of STR membership (see funcs/prob.py's indelErrorProbs).

    indel_seq/ref_after/indel_len: same arguments as classify_indel_record.
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
        # -R restricts calling to the bed file's regions, so windows with
        # no overlap contribute nothing and can skip the bam.count() call
        # entirely -- otherwise splitBamRegions ends up counting reads
        # across the whole genome even when -R covers only a small
        # fraction of it.
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
    """Filter out any region tuple spanning zero real bases: a
    (contig, start, end) tuple with start>=end, or the open-ended
    (contig, start) "rest of this contig" tail once start already sits
    at (or past) that contig's actual length. Both can arise from
    splitBamRegions's cut sites landing exactly at a contig boundary
    (see its own comments) -- a benign edge case at high thread counts
    or on small/sparse inputs, but an empty region reaching callBam would
    otherwise, for the 3-tuple case, write a zero-width coverage.bed.gz
    interval that breaks tabix's sort-order requirement, and for the
    2-tuple case fetch an inverted/out-of-range range.
    """
    filtered = []
    for region in regionSequence:
        if len(region) == 3 and region[1] >= region[2]:
            continue
        if len(region) == 2 and region[1] >= bamObject.get_reference_length(region[0]):
            continue
        filtered.append(region)
    return filtered


def splitBamRegions(bams, num, contigs, step, ref=None, regionFile=None):
    bamObject = getAlignmentObject(bams[0], "rb", ref)
    contigs_set = set(contigs)
    contigs_sorted = [_ for _ in bamObject.references if _ in contigs]

    total_read_num = 0
    # for bam in bams:
    # for stat in BAM(bam,"rb").get_index_statistics():
    # if stat.contig in contigs:
    # total_read_num += stat.total
    # chunkSize = math.ceil(total_read_num / num)
    # tidList = [bamObject.get_tid(c) for c in contigs_sorted]

    cutSite = [(0, 0)]
    break_flag = 0
    current_read_count = 0
    contig_lens = np.array([bamObject.get_reference_length(c) for c in contigs_sorted])
    window_nums = np.ceil(contig_lens / step).astype(int)
    window_nums_cumulative = np.cumsum(window_nums)
    total_reads_by_windows = np.zeros(window_nums_cumulative[-1])
    # One pool, sized to the requested thread count and reused across every
    # contig/bam below, instead of spinning up a fresh (unbounded, defaulting
    # to os.cpu_count()) Pool() per num-window chunk -- for a chromosome with
    # hundreds of windows that was hundreds of pool start/teardown cycles for
    # what is, per chunk, only `num` worth of work.
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
    cut_pos = (
        cut_inds - np.concatenate([[0], window_nums_cumulative])[cut_contigs]
    ) * step
    # window_nums's last window per contig is a partial window whenever
    # contig_len isn't an exact multiple of step, so window_index * step
    # can overshoot that contig's real length. That only actually happens
    # when cut_inds saturates at window_nums_cumulative[-1] (i.e. a
    # requested thread count -p asks for more chunks than there's read
    # data to fill -- the trailing np.arange(1, num) * chunkSize targets
    # exceed the total read count and searchsorted pins them at the very
    # last window index). Left unclamped this produces a cut site past
    # the chromosome end, which downstream becomes a region with
    # start > end (or start beyond the contig) and crashes when that
    # range is fetched. Clamping to the contig's real length instead
    # collapses those trailing, over-saturated cut sites onto the
    # legitimate end-of-contig boundary, where the dedup below then
    # merges them into one.
    cut_pos = np.minimum(cut_pos, contig_lens[cut_contigs])
    for nn in range(cut_inds.size):
        candidate = (int(cut_contigs[nn]), int(cut_pos[nn]))
        # At high thread counts (small chunkSize) a single window whose
        # read count exceeds chunkSize can have more than one of the
        # np.arange(1, num) * chunkSize thresholds land inside it, so
        # searchsorted returns the same index for consecutive nn --
        # duplicating the previous cut site here. cutSite entries become
        # (contig, start, end) region boundaries downstream (see the
        # cutSites[1:] consuming loops in Caller.py/Learn.py), so a
        # duplicate produces a zero-width start==end region, which in
        # turn writes a zero-width coverage.bed.gz interval and breaks
        # its sort order for tabix. Skipping the duplicate just means
        # this run ends up with one fewer real parallel chunk than -p
        # requested for this contig boundary; callers already tolerate
        # fewer chunks than threads.
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


def get_duplex_barcode(rec, nanoseq_bam=False):
    """Return the "bc1-bc2" duplex barcode string for an aligned read.

    Normally this is just the DB tag written by DupCaller trim. With
    nanoseq_bam, the bam instead carries NanoSeq-style per-mate RB/MB tags
    (own-read barcode / mate barcode); the DB-equivalent string is
    reconstructed as {MB}-{RB} for read 1 and {RB}-{MB} for read 2, so it
    can be split the same way DB is everywhere else.
    """
    if not nanoseq_bam:
        return rec.get_tag("DB")
    mb_tag = rec.get_tag("MB")
    rb_tag = rec.get_tag("RB")
    if rec.is_read2:
        return f"{rb_tag}-{mb_tag}"
    return f"{mb_tag}-{rb_tag}"


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
        # -R restricts calling to the bed file's regions, so only fetch
        # reads inside the bed intervals that overlap this worker's chunk
        # region, instead of every read across the whole chunk. A read can
        # overlap more than one bed interval within the same chunk (e.g.
        # adjacent/nearby targets), so dedup on (query_name, is_read1,
        # reference_start) within the chunk; cross-chunk dedup still relies
        # on the same reference_start >= region[1] rule used above, which
        # is unaffected by which bed interval a read was fetched through.
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
        ### Mask overlap of reverse read
        if processed_flag:
            mate_cigar = seq.get_tag("MC")
            cigar_m = re.findall(r"(\d+)M", mate_cigar)
            cigar_d = re.findall(r"(\d+)D", mate_cigar)
            mate_reference_length = sum([int(_) for _ in cigar_m]) + sum(
                [int(_) for _ in cigar_d]
            )
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
    elif chrom == regions_end_chrom and pos > regions_end_pos:
        return locus_bed_next
    else:
        return locus_bed


def _rescue_reason_label(
    ncov_flag, nm_flag, trim_flag, extra_flag=None, extra_label=None
):
    """Pick a single, specific rescue-reason filter label for a
    --rescue-eligible candidate that reached the masked-rescue fallback
    (snp_mask/noise_mask no longer block calling at all -- see the
    unmasked_pass_bool/unmasked_antimask switch below -- and include_mask
    is never rescuable, so by construction whatever's left blocking this
    candidate is one or more of: n_cov_mask, nm_mask, trim, and -- indel
    side only -- indel_mask). Args are booleans or boolean arrays covering
    the candidate's own position(s); a plain scalar bool works too since
    np.any() accepts one. extra_flag/extra_label covers indel_mask, which
    has no SNV analogue. Checked in a fixed priority order so a position
    blocked by more than one mask still gets exactly one label rather than
    silently picking whichever happened to be checked last.
    """
    if extra_flag is not None and np.any(extra_flag):
        return extra_label
    if np.any(ncov_flag):
        return "ncov_rescued"
    if np.any(nm_flag):
        return "nm_rescued"
    if np.any(trim_flag):
        return "trim_rescued"
    return "masked_rescued"


def _compute_read_label(rec, params):
    """Duplex-family label for rec: barcode pair plus signed template_length."""
    bc1, bc2 = get_duplex_barcode(rec, params.get("nanoSeqBam")).split("-")
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
    disagree on next_reference_start, point every minority read's mate at
    the majority (largest) mate position via rugged_reads_index, and
    pre-open its rugged_reads_pool bucket.

    Pool/index keys are (chrom, position), not bare position: a redirect
    queued here isn't guaranteed to ever be drained (e.g. the downstream
    mate never turns up in this process's region -- filtered upstream, or
    mapped outside the assigned span). An orphaned bare-position key would
    silently collide with an unrelated position on a *later* chromosome
    (every chromosome restarts its own coordinate numbering near 0),
    merging a stale read from one chromosome into another chromosome's
    batch -- which corrupts reference_mat_chrom for that flush and writes
    a non-contiguous chromosome block into the coverage bed file (fatal
    for tabix indexing). Qualifying by chromosome makes a stale entry
    inert instead of a silent cross-chromosome misattribution.
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


# ======================================================================
# Error-matrix normalization and loading
# ======================================================================
def _normalize_indel_hp_mat(mat):
    """Normalize a raw (10, 12) hp.txt count matrix (rows hp run length
    1-10+, columns ref_allele*3+(idLen+1) for idLen in {-1,0,1}) into
    per-context probabilities: within each base's own 3-column group
    (ref/del/ins counts for that base), each row is divided by its own
    sum so the three probabilities for a given (hp_len, base) add to 1.
    Then enforce monotonic non-decrease across increasing hp_len within
    each group (row n >= row n-1 elementwise) --
    a regularization against noisy per-row estimates at the sparser,
    longer-homopolymer rows, done independently per base group since
    they're otherwise unrelated contexts.
    """
    mat_new = np.zeros_like(mat, dtype=float)
    for g in range(4):
        block = mat[:, g * 3 : g * 3 + 3]
        row_sum = block.sum(axis=1, keepdims=True)
        row_nonzero = (row_sum != 0).flatten()
        block_new = np.zeros_like(block, dtype=float)
        block_new[row_nonzero, :] = block[row_nonzero, :] / row_sum[row_nonzero, :]
        for nn in range(1, block_new.shape[0]):
            current_row = block_new[nn, :]
            smaller_entries = current_row <= block_new[nn - 1, :]
            current_row[smaller_entries] = block_new[nn - 1, :][smaller_entries]
            block_new[nn, :] = current_row
        mat_new[:, g * 3 : g * 3 + 3] = block_new
    return mat_new


def _normalize_indel_str_mat(mat):
    """Normalize a raw (5, 11) str.txt count matrix (rows STR-length bin
    0="0-1"/not a real repeat through 4="40+", columns idLen+5 for idLen
    in -5..5) into per-context probabilities: each row divided by its own
    sum across all 11 columns. Real-STR rows (1-4) that end up entirely
    zero (no observations at that length bin) fall back to the pooled
    distribution across whichever real-STR rows do have data. Row 0 is a
    different population (not a real repeat at all) and is never pooled
    into or from.
    """
    row_sum = mat.sum(axis=1, keepdims=True)
    row_nonzero = (row_sum != 0).flatten()
    mat_new = np.zeros_like(mat, dtype=float)
    mat_new[row_nonzero, :] = mat[row_nonzero, :] / row_sum[row_nonzero, :]
    str_rows_total = mat[1:5, :].sum(axis=0)
    if str_rows_total.sum() > 0:
        pooled = str_rows_total / str_rows_total.sum()
        for r in range(1, 5):
            if row_sum[r, 0] == 0:
                mat_new[r, :] = pooled
    return mat_new


def regularizeErrorMat(mat, pseudocount):
    """Add `pseudocount` to every entry of mat, so no cell is ever exactly
    zero and a real observation can never be assigned zero likelihood
    downstream. NaN entries (from a 0/0 row-sum division upstream -- a
    trinuc/indel context with zero observed counts) are treated as 0
    before the pseudocount is added. Replaces the old scheme of patching
    only invalid (zero/NaN) cells with the matrix's own smallest valid
    value (or a last-resort `minerr` fallback if the whole matrix was
    invalid) -- a uniform additive pseudocount makes that special-casing,
    and the fallback, unnecessary."""
    mat = np.nan_to_num(mat, nan=0.0)
    return mat + pseudocount


def load_error_matrices(params):
    """Load and regularize the SBS (amperr/dmgerr) and indel
    (amperri/dmgerri) error-rate matrices from the files named in params,
    populating params in place with ampmat/ampmat_rev/dmgmat_top/
    dmgmat_rev_top/dmgmat_bot/dmgmat_rev_bot/ampmat_hp/dmgmat_hp/
    ampmat_str/dmgmat_str/trinuc_convert/
    trinuc2num_dict/num2trinuc_list -- exactly the fields callBam needs to
    build its L/L_indel_1bp/L_indel_len power tables. Factored out (rather
    than left inlined in callBam) so Caller.py's post-hoc per-channel FDR
    threshold refinement -- re-simulating a single context's detection
    power at a new, stricter LR threshold -- can load the identical
    matrices without re-running callBam/rescanning the bam.

    isLearn is read from params (defaulting False); in learn mode the amp/
    dmg matrices are left as zeros (they're not used for calling there).
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
    ### Load amp error matrix
    if isLearn:
        ampmat = np.zeros([64, 4])
    else:
        # amperr_file now points at a Learn.py/Caller.py-auto-learn-
        # produced SRD (single-read-damage) rate matrix
        # (funcs/learn.py's estimate_sbs_srd_rates) -- already a
        # per-row probability distribution (p/p_b1/p_b2/p_b3 summing to
        # 1, EM-fit directly against the per-read base-quality
        # distribution) rather than raw mismatch counts, so no row-sum
        # normalization here any more: the old in-situ "divide raw
        # amp.tn.txt counts by their row sum" step now happens once at
        # learn time instead of on every call-time load.
        ampmat = (
            pd.read_csv(params["amperr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )
    # ampmat_avg_error = (1 - ampmat.max(axis=1,keepdims=True))/3
    ampmat_min_error = ampmat.min(axis=1, keepdims=True)
    ampmat = np.concatenate([ampmat, ampmat_min_error], axis=1)
    ampmat = regularizeErrorMat(ampmat, 1e-6)
    params["ampmat"] = ampmat

    ampmat_rev = np.zeros([64, 4])
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            ampmat_rev[trinuc2num[trinuc], nn] = ampmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    # ampmat_rev_avg_error = (1 - ampmat_rev.max(axis=1,keepdims=True))/3
    ampmat_rev_min_error = ampmat_rev.min(axis=1, keepdims=True)
    ampmat_rev = np.concatenate([ampmat_rev, ampmat_rev_min_error], axis=1)
    params["ampmat_rev"] = ampmat_rev

    if isLearn:
        params["ampmat_hp"] = np.zeros([10, 12])
        params["ampmat_str"] = np.zeros([5, 11])
    else:
        # amperri_file's prefix, split into the two files funcs/learn.py
        # writes (see Caller.py/Learn.py's write-out).
        amperri_hp_file = params["amperri_file"].replace(".amp.id.txt", ".amp.hp.txt")
        amperri_str_file = params["amperri_file"].replace(".amp.id.txt", ".amp.str.txt")
        ampmat_hp = pd.read_csv(amperri_hp_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        ampmat_str = pd.read_csv(amperri_str_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        ampmat_hp = _normalize_indel_hp_mat(ampmat_hp)
        ampmat_hp = regularizeErrorMat(ampmat_hp, 1e-6)
        params["ampmat_hp"] = ampmat_hp

        ampmat_str = _normalize_indel_str_mat(ampmat_str)
        ampmat_str = regularizeErrorMat(ampmat_str, 1e-6)
        params["ampmat_str"] = ampmat_str

    # params["ampmat_indel_mean"] = np.mean(ampmat_indel,axis=1)
    # params["ampmat_indel_rev_mean"] = np.mean(ampmat_indel_rev,axis=1)

    ### Load damage matrix
    if isLearn:
        dmgmat = np.zeros([64, 4])
    else:
        dmgmat = (
            pd.read_csv(params["dmgerr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )
        # dmgmat += 1

    # dmgmat += 0.5
    # Same zero-row guard as ampmat above: a trinuc context with zero
    # observed counts must stay exact zero here, not become NaN via 0/0,
    # so regularizeErrorMat below can patch it correctly.
    dmgmat_row_sum = dmgmat.sum(axis=1, keepdims=True)
    dmgmat_row_nonzero = (dmgmat_row_sum != 0).flatten()
    dmgmat_normalized = np.zeros_like(dmgmat)
    dmgmat_normalized[dmgmat_row_nonzero, :] = (
        dmgmat[dmgmat_row_nonzero, :] / dmgmat_row_sum[dmgmat_row_nonzero, :]
    )
    dmgmat = dmgmat_normalized
    # dmgmat_ref_error = 1 -  dmgmat.max(axis=1, keepdims=True)
    dmgmat_ref_error = dmgmat.min(axis=1, keepdims=True)
    dmgmat = np.concatenate([dmgmat, dmgmat_ref_error], axis=1)
    dmgmat = regularizeErrorMat(dmgmat, 1e-8)

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
    # dmgmat_rev_ref_error = 1 -  dmgmat_rev.max(axis=1, keepdims=True)
    dmgmat_rev_ref_error = dmgmat_rev.min(axis=1, keepdims=True)
    dmgmat_rev = np.concatenate([dmgmat_rev, dmgmat_rev_ref_error], axis=1)
    params["dmgmat_rev_top"] = dmgmat_rev

    dmgmat_b = np.vstack((dmgmat[32:64, [1, 0, 3, 2]], dmgmat[:32, [1, 0, 3, 2]]))
    # dmgmat_b_ref_error = 1 - dmgmat_b.max(axis=1, keepdims=True)
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
        params["dmgmat_str"] = np.zeros([5, 11])
    else:
        dmgerri_hp_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.hp.txt")
        dmgerri_str_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.str.txt")
        dmgmat_hp = pd.read_csv(dmgerri_hp_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        dmgmat_str = pd.read_csv(dmgerri_str_file, sep="\t", index_col=0).to_numpy(
            dtype=float
        )
        dmgmat_hp = _normalize_indel_hp_mat(dmgmat_hp)
        dmgmat_hp = regularizeErrorMat(dmgmat_hp, 1e-8)
        params["dmgmat_hp"] = dmgmat_hp

        dmgmat_str = _normalize_indel_str_mat(dmgmat_str)
        dmgmat_str = regularizeErrorMat(dmgmat_str, 1e-8)
        params["dmgmat_str"] = dmgmat_str


# ======================================================================
# FDR / detection-power calibration
# ======================================================================
def threshold_rng(base_seed, threshold):
    """Deterministic, order-independent RNG for one simulate_power_grid
    Monte Carlo draw at a specific LR `threshold`.

    Deriving a fresh seed from (base_seed, threshold) means the same
    threshold always gets the same seed -- and thus the same simulated
    grid -- no matter which process/thread computes it or in what order.
    Results are therefore reproducible across both repeated runs (given
    the same base_seed) and different -p values: a shared per-worker RNG
    would instead make results depend on which worker happens to process
    which threshold, which itself depends on how many -p workers exist.

    Seeded from base_seed plus the threshold's first 5 significant digits
    (threshold values here are LR thresholds, typically single- to
    low-double-digit floats) -- np.random.default_rng mixes a sequence of
    ints via SeedSequence, so combining the two this way (rather than
    concatenating/adding them by hand) avoids accidental seed collisions
    between different (base_seed, threshold) pairs.

    threshold == +-inf (_refine_channel's deliberate "this channel has no
    mutation-rate evidence, reject every LR" sentinel for a zero-coverage
    channel) can't be scaled into a finite digit string -- but which seed
    gets used doesn't matter in that case anyway: simulate_power_grid's
    `(LL_B1 - LL_B2) >= threshold` comparison is already deterministically
    False (or True, for -inf) for every possible draw once threshold is
    infinite, regardless of which quality values got sampled. A fixed
    sentinel digit string is therefore exact, not an approximation.
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

    Shared by callBam's L (SBS) and L_indel_1bp/L_indel_len (indel)
    power-table construction below, and reused as-is by Caller.py to
    re-simulate one channel's grid at a new, FDR-controlled threshold
    without rescanning the bam -- both need bit-identical math, just at a
    different threshold.
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


_REFINE_WORKER_CTX = {}


def init_refine_worker(
    depth_by_trinuc,
    depth_by_hpstr,
    all_quals,
    ampmat,
    ampmat_rev,
    dmgmat_top,
    dmgmat_rev_top,
    dmgmat_bot,
    dmgmat_rev_bot,
    trinuc_convert,
    ampmat_hp,
    dmgmat_hp,
    ampmat_str,
    dmgmat_str,
    seed,
):
    """Pool initializer for the per-channel FDR-threshold refinement pool
    (refine_channel_task below): stashes the read-only context every
    channel's Eeff simulation needs as worker-global state once per
    process, instead of re-pickling it into each of the ~200 per-channel
    tasks Caller.py's do_call dispatches. trinuc2num/num2trinuc are
    rebuilt locally (build_trinuc64_order takes no args and is
    deterministic) rather than also being passed through initargs.

    `seed` is the run's base Monte Carlo seed (params["seed"], resolved
    once in Caller.py's do_call), stashed here rather than a live RNG
    instance -- _channel_eeff_at_threshold derives a fresh, threshold-
    specific RNG per call (see threshold_rng) instead of mutating one
    shared generator, so results don't depend on which worker or how many
    -p workers end up evaluating which threshold.
    """
    global _REFINE_WORKER_CTX
    trinuc2num, num2trinuc = build_trinuc64_order()
    _REFINE_WORKER_CTX = dict(
        depth_by_trinuc=depth_by_trinuc,
        depth_by_hpstr=depth_by_hpstr,
        all_quals=all_quals,
        ampmat=ampmat,
        ampmat_rev=ampmat_rev,
        dmgmat_top=dmgmat_top,
        dmgmat_rev_top=dmgmat_rev_top,
        dmgmat_bot=dmgmat_bot,
        dmgmat_rev_bot=dmgmat_rev_bot,
        trinuc_convert=trinuc_convert,
        ampmat_hp=ampmat_hp,
        dmgmat_hp=dmgmat_hp,
        ampmat_str=ampmat_str,
        dmgmat_str=dmgmat_str,
        trinuc2num=trinuc2num,
        num2trinuc=num2trinuc,
        base2num={"A": 0, "T": 1, "C": 2, "G": 3},
        seed=seed,
        n1_mask=np.minimum(*np.indices((10, 10))) >= 1,
    )


def _channel_eeff_at_threshold(kind, ctx_key, threshold):
    """Detection-power-weighted opportunity (Eeff) for one channel at one
    LR threshold. Picklable top-level dispatch on channel kind, so one
    per-channel worker task (refine_channel_task) can call this instead of
    needing a different closure per kind -- closures aren't picklable
    across a Pool.

    rng is derived fresh from (base seed, threshold) on every call via
    threshold_rng rather than reused from worker state, so the same
    threshold always simulates the same grid regardless of which worker
    or how many -p workers are running.
    """
    c = _REFINE_WORKER_CTX
    rng = threshold_rng(c["seed"], threshold)
    if kind == "sbs96":
        t_fwd, b_fwd, t_rc = ctx_key
        ref_base_idx = c["base2num"][c["num2trinuc"][t_fwd][1]]
        tc = int(c["trinuc_convert"][t_fwd, b_fwd])
        probs = (
            c["ampmat"][tc, ref_base_idx],
            c["ampmat_rev"][tc, ref_base_idx],
            c["dmgmat_top"][tc, ref_base_idx],
            c["dmgmat_rev_top"][tc, ref_base_idx],
            c["dmgmat_bot"][tc, ref_base_idx],
            c["dmgmat_rev_bot"][tc, ref_base_idx],
        )
        grid = simulate_power_grid(*probs, threshold, c["all_quals"], rng)
        n1_mask = c["n1_mask"]
        eeff = float(np.sum(c["depth_by_trinuc"][:, :, t_fwd] * grid * n1_mask))
        eeff += float(np.sum(c["depth_by_trinuc"][:, :, t_rc] * grid * n1_mask))
        return eeff
    if kind == "hp":
        # HP channels only ever exist for id_len in {-1, 1} now -- there's
        # no hp.txt column for multi-bp lengths any more (those are always
        # STR-context, see indelErrorProbs). inserted_base is passed equal
        # to `base` (the position's own reference base) since this is
        # enumerating "the true opportunity for a real same-base
        # homopolymer-extending event at a position with this base" --
        # exactly what depth_by_hpstr's base axis already represents, not
        # an arbitrary hypothetical mismatched insertion (that background
        # rate has its own str.txt row-0 context, not tied to hp_len/base
        # at all). Restricted to just this channel's own pool's 2 bases
        # (base2num order A,T,C,G -- "T" pool is A/T, "C" pool is C/G,
        # matching classify_indel_channel's own base pooling and
        # Caller.py's raw_lr_hp accumulation) rather than all 4, now that
        # each pool is its own channel with its own threshold.
        hp_len, id_len, pool = ctx_key
        total = 0.0
        for base in (0, 1) if pool == "T" else (2, 3):
            probs = indelErrorProbs(
                hp_len,
                0,
                id_len,
                base,
                base,
                c["ampmat_hp"],
                c["dmgmat_hp"],
                c["ampmat_str"],
                c["dmgmat_str"],
            )
            grid = simulate_power_grid(*probs, threshold, c["all_quals"], rng)
            total += float(
                np.sum(c["depth_by_hpstr"][:, :, base * 10 + hp_len - 1] * grid)
            )
        return total
    if kind == "str":
        # str_bin in {1,2,3,4} for real STR-length contexts (>=2bp
        # events): ref_allele/inserted_base are irrelevant there (no base
        # identity in that branch of indelErrorProbs). str_bin==0 with
        # id_len==1 is the mismatched-insertion background channel (the
        # only way to reach str.txt row 0 for a +-1bp event) -- forced by
        # passing two different dummy base values (0 != 1) so
        # indelErrorProbs takes the mismatch branch regardless; harmless
        # for the multi-bp case since that branch ignores both args.
        str_bin, id_len = ctx_key
        probs = indelErrorProbs(
            1,
            str_bin,
            id_len,
            0,
            1,
            c["ampmat_hp"],
            c["dmgmat_hp"],
            c["ampmat_str"],
            c["dmgmat_str"],
        )
        grid = simulate_power_grid(*probs, threshold, c["all_quals"], rng)
        if str_bin == 0:
            # No per-context depth bucket exists for "not a real repeat"
            # (depth_by_hpstr's 39+strs_for_row buckets only ever populate
            # for strs_for_row>=1 -- 39+0 would collide with the last real
            # HP bucket, base G/hp_len=10, not a dedicated background
            # slot). A mismatched-insertion opportunity isn't tied to any
            # hp/str context anyway (row 0's Pamp/Pdmg have no context
            # dependence), so the correct denominator is genuinely
            # "covered at this (n_top, n_bot) depth, anywhere" --
            # depth_by_trinuc summed across all 64 contexts is exactly
            # that count, already computed for the SBS side.
            depth = c["depth_by_trinuc"].sum(axis=2)
        else:
            depth = c["depth_by_hpstr"][:, :, 39 + str_bin]
        return float(np.sum(depth * grid))
    raise ValueError(f"unknown channel kind {kind!r}")


def _refine_channel(mu0, threshold0, fdr_thr):
    """Non-iterative FDR-controlled threshold for one channel.

    mu0 (this channel's round-1 MLE mixture weight, from the direct
    brentq solve in refine_channel_task) is taken as the channel's
    mutation rate outright -- no re-simulation/re-MLE step-up loop. The LR
    threshold whose local fdr 1/(1+LR*mu0) equals fdr_thr is solved for
    directly (LR = (1-fdr_thr)/(fdr_thr*mu0)); mu0 itself is also what
    every PASS call's own local FDR is computed against downstream
    (Caller.py's per-call FDR stamping).

    Clamped to never go below threshold0: round 1's bam scan only ever
    records a candidate whose LR clears threshold0 in the first place
    (funcs/call.py's LR_pass_bool gate, evaluated long before any
    per-channel refinement exists) -- a candidate below threshold0 isn't
    filtered out, it's never written into mutsAll/indelsAll at all, so
    round 2's post-hoc re-filter (no bam rescan) could never recover it
    even if mu0 alone would justify a looser threshold here.

    mu0 == 0 (this channel had zero effective coverage in round 1 -- see
    refine_channel_task's Eeff0 == 0 short-circuit) makes the local-fdr
    formula's implied LR threshold a division by zero (mu0/(1-mu0) == 0).
    There is no mutation-rate evidence for this channel at all, so instead
    of solving for a finite cutoff, report nothing as PASS in round 2 for
    it: a threshold of +inf fails every finite LR in Caller.py's
    LR < threshold re-filter.
    """
    if mu0 == 0:
        return float("inf")
    lr_threshold = (1.0 - fdr_thr) / (fdr_thr * mu0 / (1 - mu0))
    return max(threshold0, float(np.log10(lr_threshold)))


def refine_channel_task(job):
    """One independent unit of work for the FDR-threshold pool: computes
    this channel's Eeff at its default threshold, then solves directly for
    this channel's MLE mixture weight (mu0) and the LR threshold at which
    its local fdr (using mu0) equals fdr_thr, clamped to never go below
    threshold0 -- see _refine_channel. Channels (each SBS96 class / each
    HP-length x indel-length / STR-bin x indel-length combo, ~200 total)
    are fully independent -- own raw_lr list, own Eeff formula -- so
    do_call dispatches these across a Pool (one task per channel) instead
    of running them one at a time in the main process.

    mu0 solves g(mu) = sum(raw_lr/(1-mu+mu*raw_lr)) - Eeff0 + pseudocount/mu
    == 0 directly via brentq. The pseudocount/mu term sends g(0) to
    literally +inf (np.divide(pseudocount, 0.0) rather than plain float
    division, which would raise ZeroDivisionError instead); g(1) works out
    to n - Eeff0 + pseudocount (n = len(raw_lr_list), since each
    raw_lr/(1-1+1*raw_lr) term is exactly 1), so a bracketing sign change
    is only guaranteed when Eeff0 > n + pseudocount, not for every
    Eeff0 > 0. In practice this always holds -- candidate mutations are
    rare relative to a channel's total effective coverage (n/Eeff0 << 1)
    -- except when Eeff0 itself is essentially zero (rounds down to
    around 1) while the channel still has n >= 1 candidates, the one
    realistic way to violate it (without this term, g has a trivial root
    sitting at mu=0 itself).

    Eeff0 == 0 (zero effective coverage for this channel in round 1) is
    the one case brentq can't solve: every term of g is non-negative
    then, so g never crosses zero anywhere in (0, 1) and brentq raises
    for failing to bracket a root. Short-circuit to mu0 = 0 directly --
    there's no coverage to estimate a mutation rate from anyway.
    """
    name, kind, ctx_key, raw_lr_list, threshold0, fdr_thr, pseudocount = job
    Eeff0 = _channel_eeff_at_threshold(kind, ctx_key, threshold0)
    raw_lr = np.asarray(raw_lr_list, dtype=float)

    if Eeff0 == 0:
        mu0 = 0.0
    else:

        def g(mu):
            with np.errstate(divide="ignore"):
                return (
                    np.sum(raw_lr / (1.0 - mu + mu * raw_lr))
                    - Eeff0
                    + np.divide(pseudocount, mu)
                )

        mu0 = brentq(g, 0.0, 1.0)
    new_threshold = _refine_channel(mu0, threshold0, fdr_thr)
    return name, kind, ctx_key, Eeff0, mu0, new_threshold


# ======================================================================
# Depth-matrix accumulation
# ======================================================================
def _accumulate_depth_matrix(depth_mat, n_top, n_bot, category, valid, n_cat):
    """
    Add one count per valid position to depth_mat[n_top, n_bot, category],
    via a single bincount instead of a Python loop.

    depth_mat: (10, 10, n_cat) running total, updated in place.
    n_top, n_bot: (window_len,) read-family top/bottom-strand counts,
        already capped at 9 (matching the L-table convention used
        elsewhere in this module).
    category: (window_len,) context bucket index in [0, n_cat) -- e.g.
        trinuc (0-63), merged hp/str bucket (0-22), or dinuc (0-15).
        Entries outside [0, n_cat) are fine as long as they're excluded
        by `valid` (this function never reads them).
    valid: (window_len,) bool mask of positions to actually count.

    This only ever touches the small, fixed-size aggregate total -- no
    per-locus record is kept, so there is nothing to carry across window
    boundaries and nothing to merge across worker processes' regions.
    """
    if not np.any(valid):
        return
    idx = (
        n_top[valid].astype(np.int64) * 10 * n_cat
        + n_bot[valid].astype(np.int64) * n_cat
        + category[valid].astype(np.int64)
    )
    depth_mat += np.bincount(idx, minlength=10 * 10 * n_cat).reshape(10, 10, n_cat)
