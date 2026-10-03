"""Tumor depth extraction exempts a founding family's reads from minBq; the
family identifier is the duplex barcode pair (either orientation) plus
|template_length|, so a same-barcode read from a different molecule (other
template length) is still held to minBq."""

import pysam

from DupCaller_sub.funcs.depth import extractDepthBatchDbs, extractDepthBatchSnv

CHROM = "chr1"
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
REF = "ACGTTGCAAGCTTACGGATC"
START = 100
POS = 105  # 1-based; REF[4] == "T"
PARAMS = {"maxDepth": 1000, "mapq": 0, "barcodeTag": "DB"}


def _read(name, db, tlen, bq, alt="G", alt2=None):
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
    rec.mapping_quality = 60
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
