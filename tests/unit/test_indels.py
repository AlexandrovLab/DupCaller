"""Regression tests for funcs/indels.py.

Two bugs fixed together:

1. findIndels mis-walked the N cigar op (reference skip) as if it were S
   (soft clip): both advanced readPos, but N consumes reference bases,
   not query bases. Any indel downstream of an N in the same read was
   reported at the wrong genomic position.

2. getIndelArr guessed an allele for reads that can't show one: soft
   clips matching an inserted sequence counted as ALT, reads ending right
   after the anchor or partway through a deletion counted as REF, and a
   read not overlapping the candidate at all counted as a conflict
   (dropping the candidate). It now returns ALT only for the read's own
   anchored CIGAR indel, REF only for a read aligned through the whole
   locus, CONFLICT for a different indel inside the locus, and
   UNINFORMATIVE otherwise.
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
import numpy as np

from DupCaller_sub.funcs.indels import (
    INDEL_ALT,
    INDEL_CONFLICT,
    INDEL_REF,
    INDEL_UNINFORMATIVE,
)

# Non-repetitive around 104-107, so the deletion "104:-2" (removes GC at
# 105-106) has one placement and needs REF reads through 107.
REF = "ACGTTGCAAGCTTACGGATCCATGCAGTCA"
REF_INT = np.array(["ATCG".index(b) for b in REF])


def _arr(read, indel):
    return getIndelArr(read, [indel], 0, REF_INT, 100)


def _ref_read(n, start=100, soft=0):
    cigar = [(0, n)] + ([(4, soft)] if soft else [])
    return _read(
        cigar, seq=REF[start - 100 : start - 100 + n] + "T" * soft, ref_start=start
    )


def test_deletion_alt_and_ref():
    alt = _read([(0, 5), (2, 2), (0, 8)], seq=REF[:5] + REF[7:15])
    seqArr, qualArr = _arr(alt, "104:-2")
    assert seqArr[0] == INDEL_ALT and qualArr[0] > 0
    seqArr, qualArr = _arr(_ref_read(15), "104:-2")
    assert seqArr[0] == INDEL_REF and qualArr[0] > 0


def test_read_ending_at_anchor_is_uninformative():
    assert _arr(_ref_read(5, soft=10), "104:-2")[0][0] == INDEL_UNINFORMATIVE
    assert _arr(_ref_read(5), "104:-2")[0][0] == INDEL_UNINFORMATIVE


def test_deletion_not_fully_spanned_is_uninformative():
    # Must align 104..107 (deleted bases plus the next base) to be REF.
    assert _arr(_ref_read(6), "104:-2")[0][0] == INDEL_UNINFORMATIVE
    assert _arr(_ref_read(7, soft=3), "104:-2")[0][0] == INDEL_UNINFORMATIVE
    assert _arr(_ref_read(8), "104:-2")[0][0] == INDEL_REF


def test_non_overlapping_read_is_uninformative_not_conflict():
    assert _arr(_ref_read(10, start=110), "104:-2")[0][0] == INDEL_UNINFORMATIVE


def test_other_indel_in_locus_is_conflict():
    ins = _read([(0, 5), (1, 3), (0, 7)], seq=REF[:5] + "GGG" + REF[5:12])
    assert _arr(ins, "104:-2")[0][0] == INDEL_CONFLICT
    longer_del = _read([(0, 5), (2, 3), (0, 7)], seq=REF[:5] + REF[8:15])
    assert _arr(longer_del, "104:-2")[0][0] == INDEL_CONFLICT


def test_indel_outside_locus_is_not_conflict():
    far = _read([(0, 12), (2, 2), (0, 5)], seq=REF[:12] + REF[14:19])
    assert _arr(far, "104:-2")[0][0] == INDEL_REF


def test_insertion_soft_clips_are_uninformative_either_way():
    match = _read([(0, 5), (4, 2)], seq=REF[:5] + "CC")
    mismatch = _read([(0, 5), (4, 2)], seq=REF[:5] + "GG")
    assert _arr(match, "104:2:CC")[0][0] == INDEL_UNINFORMATIVE
    assert _arr(mismatch, "104:2:CC")[0][0] == INDEL_UNINFORMATIVE


def test_insertion_at_read_end_needs_anchor():
    trailing = _read([(0, 5), (1, 2)], seq=REF[:5] + "CC")
    assert _arr(trailing, "104:2:CC")[0][0] == INDEL_UNINFORMATIVE
    anchored = _read([(0, 5), (1, 2), (0, 8)], seq=REF[:5] + "CC" + REF[5:13])
    assert _arr(anchored, "104:2:CC")[0][0] == INDEL_ALT
    assert _arr(_ref_read(15), "104:2:CC")[0][0] == INDEL_REF


def test_repeat_ref_needs_whole_run_and_left_aligned_alt():
    # 100 CAAAAT...: +A anywhere in the A run left-aligns to anchor 100.
    ref = "CAAAATGCGT"
    ref_int = np.array(["ATCG".index(b) for b in ref])
    inside = _read([(0, 4)], seq=ref[:4])
    assert (
        getIndelArr(inside, ["100:1:A"], 0, ref_int, 100)[0][0] == INDEL_UNINFORMATIVE
    )
    spans = _read([(0, 6)], seq=ref[:6])
    assert getIndelArr(spans, ["100:1:A"], 0, ref_int, 100)[0][0] == INDEL_REF
    shifted = _read([(0, 3), (1, 1), (0, 5)], seq=ref[:3] + "A" + ref[3:8])
    assert getIndelArr(shifted, ["100:1:A"], 0, ref_int, 100)[0][0] == INDEL_ALT


def test_reference_position_zero_is_aligned():
    # Reference position 0 is a real aligned position, not "unaligned".
    ref_int = np.array(["ATCG".index(b) for b in "ACGTACGTAC"])
    read = _read([(0, 10)], seq="ACGTACGTAC", ref_start=0)
    assert getIndelArr(read, ["0:-1"], 0, ref_int, 0)[0][0] == INDEL_REF
