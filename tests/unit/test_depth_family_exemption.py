"""Tumor depth extraction exempts a founding family's reads from minBq and
--mapq; the family identifier is the duplex barcode pair (either orientation)
plus |template_length|, so a same-barcode read from a different molecule
(other template length) is still held to both."""

import pysam

from DupCaller_sub.funcs.depth import (
    extractDepthBatchDbs,
    extractDepthBatchIndel,
    extractDepthBatchSnv,
)

CHROM = "chr1"
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
REF = "ACGTTGCAAGCTTACGGATC"
START = 100
POS = 105  # 1-based; REF[4] == "T"
PARAMS = {"maxDepth": 1000, "mapq": 0, "barcodeTag": "DB"}


def _read(name, db, tlen, bq, alt="G", alt2=None, mapq=60):
    seq = list(REF)
    seq[POS - 1 - START] = alt
    if alt2 is not None:
        seq[POS - START] = alt2
    rec = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    rec.query_name = name
    rec.query_sequence = "".join(seq)
    rec.query_qualities = pysam.qualitystring_to_array(chr(33 + bq) * len(REF))
    rec.reference_id = 0
    rec.reference_start = START
    rec.mapping_quality = mapq
    rec.cigar = [(0, len(REF))]
    rec.template_length = tlen
    rec.set_tag("DB", db)
    return rec


def _bam(tmp_path, reads):
    path = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for r in reads:
            out.write(r)
    pysam.index(path)
    return pysam.AlignmentFile(path, "rb")


def _reads(**kw):
    return [
        _read("same_fam", "AAA-CCC", 200, 5, **kw),
        _read("mate_swapped_bc", "CCC-AAA", -200, 5, **kw),
        _read("same_bc_other_tl", "AAA-CCC", 300, 5, **kw),
        _read("other_bc", "GGG-TTT", 200, 5, **kw),
        _read("high_bq", "GGG-TTT", 250, 30, **kw),
    ]


def test_snv_exemption_requires_barcode_and_template_length(tmp_path):
    bam = _bam(tmp_path, _reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    # same_fam + mate_swapped_bc exempt, high_bq passes on its own.
    assert res[key] == (3, 0, 0, 3)


def test_snv_no_family_match_keeps_minbq(tmp_path):
    bam = _bam(tmp_path, _reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 999)}}
    )
    assert res[key] == (1, 0, 0, 1)


def test_dbs_exemption_requires_barcode_and_template_length(tmp_path):
    bam = _bam(tmp_path, _reads(alt2="A"))
    key = (CHROM, POS, "TG", "GA")
    res = extractDepthBatchDbs(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    assert res[key] == (3, 0, 0, 3)


MAPQ_PARAMS = dict(PARAMS, mapq=30)


def _mapq_reads(**kw):
    return [
        _read("same_fam", "AAA-CCC", 200, 30, mapq=0, **kw),
        _read("mate_swapped_bc", "CCC-AAA", -200, 30, mapq=5, **kw),
        _read("same_bc_other_tl", "AAA-CCC", 300, 30, mapq=0, **kw),
        _read("other_bc", "GGG-TTT", 200, 30, mapq=0, **kw),
        _read("high_mapq", "GGG-TTT", 250, 30, mapq=60, **kw),
    ]


def test_snv_mapq_exemption_requires_barcode_and_template_length(tmp_path):
    bam = _bam(tmp_path, _mapq_reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], MAPQ_PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    # same_fam + mate_swapped_bc exempt from mapq, high_mapq passes on its own.
    assert res[key] == (3, 0, 0, 3)


def test_snv_mapq_without_call_barcodes_is_strict(tmp_path):
    bam = _bam(tmp_path, _mapq_reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(bam, [key], MAPQ_PARAMS, minbq=18)
    assert res[key] == (1, 0, 0, 1)


def test_snv_family_read_needs_no_minbq_or_mapq(tmp_path):
    bam = _bam(
        tmp_path,
        [
            _read("same_fam", "AAA-CCC", 200, 5, mapq=0),
            _read("other", "GGG-TTT", 200, 5, mapq=60),
        ],
    )
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], MAPQ_PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    assert res[key] == (1, 0, 0, 1)


def test_low_mapq_non_family_read_does_not_block_its_mate(tmp_path):
    # Same query name: a low-MAPQ read outside every founding family must not
    # claim the name, so the high-MAPQ mate overlapping the column counts.
    bam = _bam(
        tmp_path,
        [
            _read("pair", "GGG-TTT", 200, 30, mapq=0),
            _read("pair", "GGG-TTT", -200, 30, mapq=60),
        ],
    )
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], MAPQ_PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    assert res[key] == (1, 0, 0, 1)


def test_indel_mapq_exemption(tmp_path):
    def ins(name, db, tlen, mapq):
        rec = _read(name, db, tlen, 30, alt="T", mapq=mapq)
        seq = REF[: POS - START] + "A" + REF[POS - START :]
        rec.query_sequence = seq
        rec.query_qualities = pysam.qualitystring_to_array("?" * len(seq))
        rec.cigar = [(0, POS - START), (1, 1), (0, len(REF) - (POS - START))]
        return rec

    bam = _bam(
        tmp_path,
        [
            ins("same_fam", "AAA-CCC", 200, 0),
            ins("same_bc_other_tl", "AAA-CCC", 300, 0),
            ins("high_mapq", "GGG-TTT", 250, 60),
        ],
    )
    key = (CHROM, POS, "T", "TA")
    res = extractDepthBatchIndel(
        bam, [key], MAPQ_PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    assert res[key] == (2, 0, 0, 2)


def test_dbs_mapq_exemption(tmp_path):
    bam = _bam(tmp_path, _mapq_reads(alt2="A"))
    key = (CHROM, POS, "TG", "GA")
    res = extractDepthBatchDbs(
        bam, [key], MAPQ_PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200)}}
    )
    assert res[key] == (3, 0, 0, 3)
