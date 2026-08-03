import math
import os
import time
from multiprocessing import Pool
import numpy as np
import h5py
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED


def _ensure_type_subdirs(prefix):
    """Per-mutation-type output subfolders under a sample's prefix
    directory — SBS/, INDEL/, DBS/ — created (if missing). Every
    SBS/indel/DBS-specific output (VCFs, corrected tables, burden files,
    plots, by-read-group breakdowns) lives under the matching subfolder;
    shared/whole-sample files (_stats.txt, _coverage.bed.gz,
    _duplex_allele_counts.txt, _*_by_duplex_group.txt, etc.) stay at the
    top level of prefix. Shared by Caller.py (writer of the VCFs) and
    Estimate.py (writer of everything else, and reader of the VCFs) so
    both agree on the exact same subfolder layout.

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
]


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

    Returns the (100,) accumulated opportunity-count array directly.
    (Previously returned (col_idx, weight) arrays for the caller to run
    through np.bincount(..., minlength=100) — but every contribution here
    has weight exactly 1, and up to ~20 separate whole-chromosome-length
    index/weight array pairs had to stay alive simultaneously (one per
    call below) until one final concatenate at the end. For a real
    chromosome that pattern was the dominant memory cost of computing
    reference-genome indel composition — see calculate_ref_indel100. Each
    contribution is now folded straight into the running (100,) total
    instead: a scalar increment for columns shared by every masked
    position, or a small per-call np.bincount for columns that vary
    per-position — numerically identical since every weight was always 1,
    just without the concatenated-array detour.)
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

    # Flat rep1 opportunity at non-STR-annotated (~real_str) positions:
    # every such position is a candidate for an arbitrary, non-repeating
    # U-bp deletion (repeat_count=1) for each unit size 2/3/4/5 — a fixed,
    # uniform credit, NOT an attempt to discover a true higher count that
    # might coincidentally exist there (that would need rescanning the
    # reference, which is exactly what call.py's matching per-family
    # credit deliberately doesn't do either — see its own comment). Kept
    # deterministic and scan-free so this side and the coverage side
    # (call.py's cov_mat_indel/indel100 accumulation) agree by
    # construction rather than by two independently-approximated scans
    # happening to line up.
    not_real_str = ~real_str
    for u_idx in range(4):
        add_uniform(44 + u_idx * 6, not_real_str)

    # 10: STR insertion opportunity (cols 68-91) — real part: one credit
    # per annotated tract (real_str & str_cut), with the same rep-count
    # 0/1 ambiguity hedge as the homopolymer deletion case above (STR
    # insertion still has a real "0" bin, unlike homopolymer insertion).
    real_str_cut = real_str & str_cut
    count_bucket_ins = np.clip(str_repeat_count_arr, 0, 5)
    add_varying(68 + unit_bucket * 6 + count_bucket_ins, real_str_cut)
    is_rep1 = str_repeat_count_arr == 1
    add_varying(68 + unit_bucket * 6 + 0, real_str_cut & is_rep1)

    # Flat rep0/rep1 opportunity at non-STR-annotated positions ("what if
    # a U-bp unit were inserted here"): same uniform, scan-free credit as
    # the deletion case above, into both the "0" (novel, no repeat yet)
    # and "1" (first extension) buckets for every unit size.
    for u_idx in range(4):
        base_col = 68 + u_idx * 6
        add_uniform(base_col, not_real_str)
        add_uniform(base_col + 1, not_real_str)

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
    back to the self-derived homopolymer run (unit_len implicitly 1) —
    mirrors the old apply_repeat_bed's last-write-wins priority, where
    --strbed/--repeatBed entries always described actual annotated repeats
    and there was no separate "background" source to conflict with.

    hp_dataset, str_dataset: h5py Datasets for one chromosome from hp.h5
        and str.h5 respectively (e.g. hp_h5[chrom], str_h5[chrom]).
    start, end: optional slice bounds; defaults to the whole chromosome.

    Returns a (3, L) int16 array: row0=unit_len, row1=cut, row2=repeat_count.
    Values never exceed 127 (Index.py clips them before writing), so int16
    is already generous — it's only used instead of int8/uint8 to leave
    headroom for downstream arithmetic (e.g. subtracting a small constant)
    without any risk of silent unsigned wraparound. This used to be a
    plain `int` (numpy int64) via `.astype(int)` on every row, an 8x
    widening of the on-disk uint8 data for no reason — for a whole
    chromosome (the common case for Estimate.py's genome-composition
    calls, as opposed to call.py's windowed calls) that was several GB of
    unnecessary temporary allocation per call.
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


import pysam


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


def bamReadCount(bam, chrom, start, end, ref):  # , regionFile):
    # if not regionFile:
    count = getAlignmentObject(bam, "rb", ref).count(chrom, start, end)
    """
    else:
        regionFile = BED(regionFile,parser = pysam.asBed())
        count = 0
        if chrom in regionFile.contigs:
            for interval in regionFile.fetch(chrom,start,end):
                count += getAlignmentObject(bam, "rb", ref).count(interval.contig, interval.start, interval.end)
        else:
            count += 0
    """
    return count


def splitBamRegions(bams, num, contigs, step, ref):  # , regionFile):
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
                        # regionFile,
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
    for nn in range(cut_inds.size):
        cutSite.append((cut_contigs[nn], cut_pos[nn]))
    return cutSite, chunkSize, contigs_sorted


"""
Detect if the input is a bam or cram and read accordingly. Commit by Luka Culibrk (github:@lculibrk)
"""


## Allow reading either bam or cram as bamObject
def getAlignmentObject(bam, mode, refpath=None):
    if bam.endswith(".bam"):
        bamObject = BAM(bam, mode)
    elif bam.endswith(".cram"):
        bamObject = BAM(bam, mode, reference_filename=refpath)
    else:
        raise NameError(f"{bam} should have .bam or .cram as file extension")
    return bamObject
