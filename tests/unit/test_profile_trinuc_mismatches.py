"""profileTriNucMismatches with a read whose only "mismatch" is an N
(seq_mat sentinel 4) carrying a non-zero base quality above minBq.

seq_mat encodes A/T/C/G as 0-3 and 4 for anything else. The per-read SBS
loop treats a 4 as a mismatch; the BQ-histogram accumulation must keep it out
of the 4*96 ATCG range so the `.reshape([4, 96, NUM_BQ])` succeeds.
Quality-0 Ns never reach this code (`qual_masked > 0` filter).
"""
import numpy as np
import pysam

from DupCaller_sub.funcs.learn import profileTriNucMismatches

CHROM = "chr1"
CONTIG_LEN = 1000
POS = 500  # single-position window under test


def _read(header, name, is_read1, base, qual):
    rec = pysam.AlignedSegment(header)
    rec.query_name = name
    rec.query_sequence = base
    rec.query_qualities = [qual]
    rec.reference_id = 0
    rec.reference_start = POS
    rec.mapping_quality = 60
    rec.cigar = [(0, 1)]
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False  # both mates forward -> read1s are F1R2, read2s F2R1
    return rec


def _params():
    return {
        "trinuc2num_dict": {},
        "minBq": 10,
        "minRef": 1,
        "minAlt": 1,
    }


def test_nonzero_quality_n_does_not_crash_bq_hist():
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}],
        }
    )

    # F1R2 (read1, forward): two clean ref calls + one read whose only
    # "mismatch" is an N with a non-zero (> minBq) quality.
    f1r2 = [
        _read(header, "clean1", True, "A", 30),
        _read(header, "clean2", True, "A", 30),
        _read(header, "n_read", True, "N", 25),
    ]
    # F2R1 (read2, forward): plain clean family, just to clear the
    # per-strand min-group-size floor.
    f2r1 = [
        _read(header, "mate1", False, "A", 30),
        _read(header, "mate2", False, "A", 30),
        _read(header, "mate3", False, "A", 30),
    ]

    result = profileTriNucMismatches(
        seqs=f1r2 + f2r1,
        reference_start=POS,
        reference_int=np.array([0]),  # reference base at POS is 'A'
        trinuc_int=np.array([0]),
        hp_raw=np.zeros([2, 1]),
        str_raw=np.zeros([3, 1]),
        antimask=np.array([True]),
        params=_params(),
    )

    sbs_alt_bq_hist = result[-1]
    assert sbs_alt_bq_hist.shape == (64, 4, 94)
