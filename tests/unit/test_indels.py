"""Regression tests for funcs/indels.py.

Two bugs fixed together:

1. findIndels mis-walked the N cigar op (reference skip) as if it were S
   (soft clip): both advanced readPos, but N consumes reference bases,
   not query bases. Any indel downstream of an N in the same read was
   reported at the wrong genomic position.

2. getIndelArr's deletion branch scored an unaligned base immediately
   after a deletion candidate's anchor as -1 (prob.py's multiallele-
   conflict sentinel) unconditionally -- including when the read simply
   ends there (a soft clip), which has nothing to do with the candidate.
   Since -1 from a single read drops the whole candidate via
   mask_multiallele, an unrelated trailing soft clip on one read could
   silently kill a real deletion call for the whole family. The
   insertion branch already special-cased this; the deletion branch now
   does too, scoring 0 (uninformative) when the unaligned base is
   outside the read's aligned span.
"""
import pysam

from DupCaller_sub.funcs.indels import findIndels, getIndelArr

CHROM = "chr1"
CONTIG_LEN = 1000
HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}]}


def _read(cigar, seq=None, ref_start=100, quals=None):
    rec = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    length = sum(n for op, n in cigar if op in (0, 1, 4))
    rec.query_name = "r1"
    rec.query_sequence = seq if seq is not None else "A" * length
    rec.query_qualities = (
        quals if quals is not None else pysam.qualitystring_to_array("I" * length)
    )
    rec.reference_id = 0
    rec.reference_start = ref_start
    rec.mapping_quality = 60
    rec.cigar = cigar
    return rec


# ---------------------------------------------------------------- findIndels


def test_deletion_and_insertion_reported_at_correct_anchor():
    # 10M, 3bp deletion, 5M -> anchor is the last matched base before the
    # event (VCF POS convention): reference_start + 10 - 1.
    read = _read([(0, 10), (2, 3), (0, 5)])
    assert findIndels(read) == ["109:-3"]


def test_n_op_advances_reference_not_query():
    # 5M, 100bp reference skip (e.g. a spliced alignment), 5M, 1bp
    # insertion, 4M. The insertion's anchor must account for the 100
    # skipped reference bases, not be computed as if N had consumed
    # query bases like a soft clip would.
    read = _read([(0, 5), (3, 100), (0, 5), (1, 1), (0, 4)], seq="A" * 15)
    # ref consumed before the insertion: 5 (M) + 100 (N) + 5 (M) = 110
    assert findIndels(read) == ["209:1:A"]


def test_soft_clip_still_advances_query_only():
    # 5S, 5M, 1bp insertion, 4M -- the leading soft clip shifts the read
    # index but not the reference index.
    read = _read([(4, 5), (0, 5), (1, 1), (0, 4)], seq="A" * 15)
    assert findIndels(read) == ["104:1:A"]


# -------------------------------------------------------------- getIndelArr


def test_deletion_soft_clipped_right_after_anchor_is_uninformative():
    """A read that ends (soft clip) immediately after a deletion
    candidate's anchor must not be flagged as a competing allele (-1) --
    that sentinel would drop the whole candidate for every read via
    prob.py's mask_multiallele, even though this read simply doesn't
    extend far enough to say anything about the candidate."""
    read = _read([(0, 5), (4, 10)], seq="A" * 15, ref_start=100)
    seqArr, _ = getIndelArr(read, ["104:-2"], min_bq=0)
    assert seqArr[0] == 0


def test_deletion_real_competing_insertion_after_anchor_is_conflict():
    """Sanity check the fix isn't over-broad: a genuine competing indel
    (an insertion right where the candidate deletion would start) is
    still a real multiallele conflict and must stay -1."""
    read = _read([(0, 5), (1, 3), (0, 7)], seq="A" * 15, ref_start=100)
    seqArr, _ = getIndelArr(read, ["104:-2"], min_bq=0)
    assert seqArr[0] == -1


def test_deletion_matching_read_scores_alt():
    read = _read([(0, 5), (2, 2), (0, 8)], seq="A" * 13, ref_start=100)
    seqArr, qualArr = getIndelArr(read, ["104:-2"], min_bq=0)
    assert seqArr[0] == 1
    assert qualArr[0] > 0


def test_deletion_ref_read_scores_ref():
    read = _read([(0, 15)], seq="A" * 15, ref_start=100)
    seqArr, _ = getIndelArr(read, ["104:-2"], min_bq=0)
    assert seqArr[0] == 0


def test_insertion_soft_clip_window_checked_against_inserted_seq():
    """Pre-existing behavior this change mirrors: an insertion candidate
    whose inserted bases fall in a trailing soft clip is scored by
    whether the clipped bases actually match the candidate sequence, not
    unconditionally treated as a conflict."""
    read_match = _read([(0, 5), (4, 2)], seq="AAAAA" + "CC", ref_start=100)
    seqArr, _ = getIndelArr(read_match, ["104:2:CC"], min_bq=0)
    assert seqArr[0] == 1

    read_mismatch = _read([(0, 5), (4, 2)], seq="AAAAA" + "GG", ref_start=100)
    seqArr, _ = getIndelArr(read_mismatch, ["104:2:CC"], min_bq=0)
    assert seqArr[0] == 0
