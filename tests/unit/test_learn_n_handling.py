"""How profileTriNucMismatches treats N/off-read bases within a read
family, per the current (intentional) design:

1. There is NO EXPLICIT family-wide antimask veto for a position just
   because one member read has an N there (that broader, direct veto --
   `F1R2_antimask[(F1R2_seq_mat == 4).any(axis=0)] = False` -- has been
   deliberately left disabled since commit ba964c4, which scoped the
   N-handling fix to a narrower per-read exclusion instead of
   re-enabling the family-wide one).
2. BUT an N can still indirectly mask a position for the WHOLE group:
   F1R2_count_mat only counts real ATCG bases (`for nn in range(4):
   F1R2_count_mat[nn] = (F1R2_seq_mat == nn).sum(...)`), never N. So
   `F1R2_antimask[F1R2_count_sum < min_depth] = False` silently drops a
   position for every read in the group whenever too few reads have a
   real base there -- if the family is small enough (at or near
   min_depth), a single read's N is enough to trip this.
3. When the group's depth stays comfortably above min_depth despite one
   read's N, that read's OTHER positions still contribute normally; only
   the N position itself is excluded on that one read
   (`valid = (qual_masked > 0) & (seq_masked < 4)`), not vetoed for the
   whole group. This holds regardless of how many N's/mismatches a read
   carries -- there is no per-read mismatch-count gate any more.
"""
import numpy as np
import pysam

from DupCaller_sub.funcs import learn as learnmod
from DupCaller_sub.funcs.misc import build_trinuc64_order

READ_LEN = 10
QUAL = 37

trinuc2num, num2trinuc = build_trinuc64_order()
base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
ROW = trinuc2num["ACA"]
REF_BASE = base2num["ACA"[1]]  # 'C'

_bam = pysam.AlignmentFile(
    "/dev/null", "wb", header={"SQ": [{"SN": "chr1", "LN": 1000}]}
)
_QUALS = pysam.qualitystring_to_array("".join(chr(QUAL + 33) for _ in range(READ_LEN)))


def _make_read(bases, is_read1):
    """bases: string of length READ_LEN over 'ATCGN'."""
    rec = pysam.AlignedSegment(_bam.header)
    rec.query_name = "r"
    rec.flag = (
        1 | 2 | (64 if is_read1 else 128)
    )  # paired, proper_pair, read1/read2, forward
    rec.reference_id = 0
    rec.reference_start = 1000
    rec.cigarstring = f"{READ_LEN}M"
    rec.query_sequence = bases
    rec.query_qualities = _QUALS
    return rec


def _profile(f1r2_seqs, f2r1_seqs):
    n = READ_LEN
    reference_int = np.full(n, REF_BASE, dtype=int)
    trinuc_int = np.full(n, ROW, dtype=int)
    hp_raw = np.zeros((2, n))
    str_raw = np.zeros((3, n))
    antimask = np.ones(n, dtype=bool)
    params = {
        "trinuc2num_dict": trinuc2num,
        "minBq": 18,
        "minAltQual": 90,
        "minRef": 3,
        "minAlt": 3,
    }
    seqs = [_make_read(b, is_read1=True) for b in f1r2_seqs] + [
        _make_read(b, is_read1=False) for b in f2r1_seqs
    ]
    mismatch_profile, _, _, _, _, _, sbs_hist = learnmod.profileTriNucMismatches(
        seqs, 1000, reference_int, trinuc_int, hp_raw, str_raw, antimask, params
    )
    return mismatch_profile, sbs_hist


REF_READ = "C" * READ_LEN  # every position matches reference ("ACA"'s ref base)


def test_single_n_excludes_only_that_read_at_that_position():
    """Family well above min_depth (4 clean reads + 1 with an N), so the
    N can't drop group depth below the floor -- isolates the pure
    per-read exclusion (point 3 above) from the depth-interaction
    (point 2), which the next test covers separately."""
    one_n_at_5 = "CCCCCNCCCC"  # exactly one N, at index 5
    f1r2 = [
        one_n_at_5,
        REF_READ,
        REF_READ,
        REF_READ,
        REF_READ,
    ]  # depth=4 ATCG even with the N
    f2r1 = [REF_READ, REF_READ, REF_READ]

    mismatch_profile, sbs_hist = _profile(f1r2, f2r1)

    ref_count = mismatch_profile[ROW, REF_BASE]
    alt_count = mismatch_profile[ROW].sum() - ref_count
    n_counted = sbs_hist[ROW].sum()

    # 8 reads x 10 positions - 1 (the single N position, excluded only on
    # the one read that has it) = 79.
    assert n_counted == 79
    assert ref_count == 79
    assert alt_count == 0


def test_n_can_indirectly_mask_the_whole_group_at_low_depth():
    """Family sitting exactly at min_depth=3: F1R2_count_mat never counts
    N, so one read's N drops the group's real-base depth at that
    position to 2 (< min_depth), tripping `F1R2_antimask[F1R2_count_sum <
    min_depth] = False` for ALL THREE F1R2 reads at that position -- not
    just the one with the N. This is a real, currently-live side effect,
    distinct from the (deliberately disabled) direct per-position N
    veto."""
    one_n_at_5 = "CCCCCNCCCC"
    f1r2 = [one_n_at_5, REF_READ, REF_READ]  # exactly 3 = min_depth
    f2r1 = [REF_READ, REF_READ, REF_READ]

    mismatch_profile, sbs_hist = _profile(f1r2, f2r1)

    ref_count = mismatch_profile[ROW, REF_BASE]
    n_counted = sbs_hist[ROW].sum()

    # F1R2: position 5 is masked for the whole group (depth drops to 2
    # there) -> all 3 F1R2 reads lose that one position -> 3*9=27.
    # F2R1: untouched (its own depth stays 3 throughout) -> 3*10=30.
    # Total 57 -- NOT 59 (what you'd get if only the N-bearing read's own
    # position were excluded, per the previous test's family size).
    assert n_counted == 57
    assert ref_count == 57


def test_real_alt_base_still_credited_alongside_reference_opportunity():
    one_alt = "CCCCCACCCC"  # one real A/C-context mismatch at index 5, no N's
    f1r2 = [one_alt, REF_READ, REF_READ, REF_READ, REF_READ]
    f2r1 = [REF_READ, REF_READ, REF_READ]

    mismatch_profile, sbs_hist = _profile(f1r2, f2r1)

    ref_count = mismatch_profile[ROW, REF_BASE]
    alt_count = mismatch_profile[ROW, base2num["A"]]
    n_counted = sbs_hist[ROW].sum()

    # No N's anywhere, so every one of the 80 positions is a real base:
    # the mismatching read still contributes its 9 matching positions as
    # reference opportunity plus 1 as the real alt observation.
    assert n_counted == 80
    assert ref_count == 79
    assert alt_count == 1
