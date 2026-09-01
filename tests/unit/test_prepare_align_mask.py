"""Regression test for funcs/depth.py's prepareAlignMask: the per-read NM
tag must be indel-corrected (indel bases backed out, one mismatch-
equivalent credited per indel event -- same correction as
funcs/call.py's read-blacklist NM_no_id) before it's averaged into the
NM/AS-XS alignment mask. A prior regression left the CIGAR loop that
computes this correction dead (id_length/id_num were initialized but
never incremented, and the block that referenced them didn't even run --
it raised NameError on any real read), so NM_no_id silently fell back to
the raw NM tag, inflating the apparent mismatch count for any
indel-containing read and risking over-masking real loci.
"""
import pysam

from DupCaller_sub.funcs.depth import prepareAlignMask

CHROM = "chr1"
CONTIG_LEN = 1000
REGION_START = 100
REGION_END = 150
MID = 25  # REGION_START + MID, inside every read's span below

# maxNM chosen so a deletion's raw NM (indel bases counted individually)
# and its indel-corrected NM_no_id (one mismatch-equivalent per indel
# event) land on opposite sides of the align_mask threshold (maxNM/2):
# raw NM=3 -> masked, corrected NM_no_id=1 -> not masked.
DELETION_LEN = 3
MAX_NM = 4
PARAMS = {"maxNM": MAX_NM, "minMeanASXS": 0, "mapq": 0}


def _aligned_segment(header, name, is_read1, cigar, nm, as_tag=100, xs_tag=0):
    seq_len = sum(length for op, length in cigar if op in (0, 1, 4))
    rec = pysam.AlignedSegment(header)
    rec.query_name = name
    rec.query_sequence = "A" * seq_len
    rec.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    rec.reference_id = 0
    rec.reference_start = REGION_START
    rec.mapping_quality = 60
    rec.cigar = cigar
    rec.is_paired = True
    rec.is_proper_pair = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = not is_read1
    rec.mate_is_reverse = is_read1
    rec.next_reference_id = 0
    rec.next_reference_start = REGION_START
    rec.template_length = 50 if is_read1 else -50
    rec.set_tag("NM", nm)
    rec.set_tag("AS", as_tag)
    rec.set_tag("XS", xs_tag)
    return rec


def _build_bam(path, fwd_cigar, fwd_nm, rev_cigar, rev_nm):
    """One read pair (forward + reverse mate, both spanning
    [REGION_START, REGION_START+50)) with the given CIGAR/NM combos."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}],
    }
    unsorted_path = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as bam:
        fwd = _aligned_segment(bam.header, "pair1", True, fwd_cigar, fwd_nm)
        rev = _aligned_segment(bam.header, "pair1", False, rev_cigar, rev_nm)
        bam.write(fwd)
        bam.write(rev)
    pysam.sort("-o", str(path), unsorted_path)
    pysam.index(str(path))


def test_indel_correction_prevents_over_masking(tmp_path):
    """A perfect-match forward mate plus a reverse mate whose only edit
    is one DELETION_LEN-bp deletion (no real mismatches): raw NM tag
    equals DELETION_LEN, but the correct NM_no_id is exactly 1. With
    MAX_NM/2 == 2, the corrected value must NOT mask; the raw value
    would."""
    bam_path = tmp_path / "test.bam"
    _build_bam(
        bam_path,
        fwd_cigar=[(0, 50)],
        fwd_nm=0,
        rev_cigar=[(0, 20), (2, DELETION_LEN), (0, 30)],
        rev_nm=DELETION_LEN,
    )

    align_mask = prepareAlignMask(
        str(bam_path), CHROM, REGION_START, REGION_END, PARAMS
    )

    assert not align_mask[MID], (
        "position was masked -- indel-corrected NM_no_id likely regressed "
        "to the raw (uncorrected) NM tag"
    )


def test_real_mismatches_without_indels_still_mask(tmp_path):
    """Sanity check that the indel correction doesn't neuter ordinary
    mismatch-based masking: a reverse mate with DELETION_LEN real
    mismatches and NO indels has NM_no_id == raw NM == DELETION_LEN
    (id_length/id_num both 0, so the correction is a no-op), which must
    still clear the mask at MAX_NM/2."""
    bam_path = tmp_path / "test_mismatch.bam"
    _build_bam(
        bam_path,
        fwd_cigar=[(0, 50)],
        fwd_nm=0,
        rev_cigar=[(0, 50)],
        rev_nm=DELETION_LEN,
    )

    align_mask = prepareAlignMask(
        str(bam_path), CHROM, REGION_START, REGION_END, PARAMS
    )

    assert align_mask[MID]
