"""A base is usable iff BQ >= minBq, everywhere: the read matrices
(genotyping and learning zero BQ < minBq), the indel median-BQ gate, and the
pileup depth paths (pysam min_base_quality / samtools -Q keep >= minBq)."""

import numpy as np
import pysam

from DupCaller_sub.funcs.indels import INDEL_LOW_BQ, INDEL_REF, getIndelArr

CHROM = "chr1"
HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
REF = "ACGTTGCAAGCTTACGGATCCATGCAGTCA"
REF_INT = np.array(["ATCG".index(b) for b in REF])


def _ref_read(bq):
    rec = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    rec.query_name = "r1"
    rec.query_sequence = REF[:15]
    rec.query_qualities = pysam.qualitystring_to_array(chr(33 + bq) * 15)
    rec.reference_id = 0
    rec.reference_start = 100
    rec.mapping_quality = 60
    rec.cigar = [(0, 15)]
    return rec


def test_indel_median_bq_equal_to_min_bq_is_usable():
    seq_arr, _ = getIndelArr(_ref_read(18), ["104:-2"], 18, REF_INT, 100)
    assert seq_arr[0] == INDEL_REF


def test_indel_median_bq_below_min_bq_is_low_bq():
    seq_arr, _ = getIndelArr(_ref_read(17), ["104:-2"], 18, REF_INT, 100)
    assert seq_arr[0] == INDEL_LOW_BQ


def test_indel_read_without_anchor_is_low_bq_even_at_min_bq_zero():
    # The anchor (104) is not in a read starting at 105.
    rec = _ref_read(30)
    rec.reference_start = 105
    rec.query_sequence = REF[5:20]
    rec.query_qualities = pysam.qualitystring_to_array("?" * 15)
    seq_arr, _ = getIndelArr(rec, ["104:-2"], 0, REF_INT, 100)
    assert seq_arr[0] != INDEL_REF
