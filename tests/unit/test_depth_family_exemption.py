"""Tumor depth extraction exempts a founding family's reads from minBq and
--mapq; the family identifier is the duplex barcode pair (either orientation)
plus |template_length| plus the fragment's leftmost start, so a same-barcode
read from a different molecule (other template length, or same template
length but another fragment start) is still held to both."""

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


def _read(name, db, tlen, bq, alt="G", alt2=None, mapq=60, start=START, mate=START):
    """Read starting at `start` (covering POS); a -TL read's mate (the
    fragment's leftmost read) starts at `mate`."""
    seq = list(REF[start - START :])
    seq[POS - 1 - start] = alt
    if alt2 is not None:
        seq[POS - start] = alt2
    rec = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    rec.query_name = name
    rec.query_sequence = "".join(seq)
    rec.query_qualities = pysam.qualitystring_to_array(chr(33 + bq) * len(seq))
    rec.reference_id = 0
    rec.reference_start = start
    rec.mapping_quality = mapq
    rec.cigar = [(0, len(seq))]
    rec.template_length = tlen
    rec.next_reference_id = 0
    rec.next_reference_start = mate
    rec.set_tag("DB", db)
    return rec


def _bam(tmp_path, reads):
    path = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
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
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200, START)}}
    )
    # same_fam + mate_swapped_bc exempt, high_bq passes on its own.
    assert res[key] == (3, 0, 0, 3)


def test_snv_no_family_match_keeps_minbq(tmp_path):
    bam = _bam(tmp_path, _reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 999, START)}}
    )
    assert res[key] == (1, 0, 0, 1)


def test_dbs_exemption_requires_barcode_and_template_length(tmp_path):
    bam = _bam(tmp_path, _reads(alt2="A"))
    key = (CHROM, POS, "TG", "GA")
    res = extractDepthBatchDbs(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200, START)}}
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
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
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
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
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
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
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
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
    )
    assert res[key] == (2, 0, 0, 2)


def test_dbs_mapq_exemption(tmp_path):
    bam = _bam(tmp_path, _mapq_reads(alt2="A"))
    key = (CHROM, POS, "TG", "GA")
    res = extractDepthBatchDbs(
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
    )
    assert res[key] == (3, 0, 0, 3)


# Barcode + |TL| collisions: same TAG1/TAG2 and template length, different
# fragment start (a different molecule) -- never exempt.
def _collision_reads(**kw):
    return [
        _read("founding", "AAA-CCC", 200, 5, **kw),
        _read("founding_mate", "CCC-AAA", -200, 5, **kw),
        # +TL read of another molecule starting 2bp later.
        _read("other_start_fwd", "AAA-CCC", 200, 5, start=START + 2, **kw),
        # -TL read whose mate (fragment start) is elsewhere.
        _read("other_start_mate", "CCC-AAA", -200, 5, mate=START - 50, **kw),
    ]


def test_snv_barcode_tl_collision_other_fragment_start_not_exempt(tmp_path):
    bam = _bam(tmp_path, _collision_reads())
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200, START)}}
    )
    assert res[key] == (2, 0, 0, 2)


def test_snv_collision_mapq_exemption_only_founding(tmp_path):
    bam = _bam(tmp_path, _collision_reads(mapq=0))
    key = (CHROM, POS, "T", "G")
    res = extractDepthBatchSnv(
        bam,
        [key],
        dict(MAPQ_PARAMS),
        minbq=0,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
    )
    assert res[key] == (2, 0, 0, 2)


def test_indel_collision_mapq_exemption_only_founding(tmp_path):
    def ins(name, db, tlen, start=START, mate=START):
        rec = _read(name, db, tlen, 30, alt="T", mapq=0, start=start, mate=mate)
        off = POS - start
        seq = REF[start - START : POS - START] + "A" + REF[POS - START :]
        rec.query_sequence = seq
        rec.query_qualities = pysam.qualitystring_to_array("?" * len(seq))
        rec.cigar = [(0, off), (1, 1), (0, len(seq) - off - 1)]
        return rec

    bam = _bam(
        tmp_path,
        [
            ins("founding", "AAA-CCC", 200),
            ins("other_start_fwd", "AAA-CCC", 200, start=START + 2),
            ins("other_start_mate", "CCC-AAA", -200, mate=START - 50),
        ],
    )
    key = (CHROM, POS, "T", "TA")
    res = extractDepthBatchIndel(
        bam,
        [key],
        MAPQ_PARAMS,
        minbq=18,
        call_barcodes={key: {("AAA", "CCC", 200, START)}},
    )
    assert res[key] == (1, 0, 0, 1)


def test_dbs_barcode_tl_collision_not_exempt(tmp_path):
    bam = _bam(tmp_path, _collision_reads(alt2="A"))
    key = (CHROM, POS, "TG", "GA")
    res = extractDepthBatchDbs(
        bam, [key], PARAMS, minbq=18, call_barcodes={key: {("AAA", "CCC", 200, START)}}
    )
    assert res[key] == (2, 0, 0, 2)


def test_family_fragment_starts_and_call_barcode_ids():
    from DupCaller_sub.funcs.call import _collect_call_barcode
    from DupCaller_sub.funcs.depth import family_fragment_starts

    # A -TL family (downstream reads) identifies its molecule by the mates'
    # start, so its reads and their +TL mates share one identifier.
    down = [
        _read("r1", "CCC-AAA", -200, 30, start=START + 2, mate=START),
        _read("r2", "AAA-CCC", -200, 30, start=START + 2, mate=START),
    ]
    assert family_fragment_starts(down) == (START,)
    mut = {"infos": {"TAG1": "AAA", "TAG2": "CCC", "TL": -200, "FS": (START,)}}
    ids = {}
    _collect_call_barcode(ids, "k", mut)
    assert ids == {"k": {("AAA", "CCC", 200, START)}}
