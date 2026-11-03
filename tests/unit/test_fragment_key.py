"""call --fragment-key: _fragment_ends gives a read pair's unclipped 5' fragment ends (read CIGAR + mate MC tag),
and _compute_read_label keys families on them instead of the signed template length, so both mates of a pair and
reads of one molecule clipped differently share a label. Without the flag (or without an MC tag) the label is the
stock barcode + signed template length."""
import pysam

from DupCaller_sub.funcs.misc import _compute_read_label, _fragment_ends

HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]})
FK = {"fragment_key": True}
STOCK = {"fragment_key": False}


def _read(name, read1, forward, start, cigar, mate_start, mate_cigar, tlen, barcode="AAA-CCC"):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    rec.flag = 0x1 | 0x2 | (0x40 if read1 else 0x80) | (0 if forward else 0x10) | (0x20 if forward else 0)
    rec.reference_id = rec.next_reference_id = 0
    rec.reference_start = start
    rec.cigarstring = cigar
    rec.query_sequence = "A" * rec.infer_query_length()
    rec.query_qualities = pysam.qualitystring_to_array("I" * rec.infer_query_length())
    rec.mapping_quality = 60
    rec.next_reference_start = mate_start
    rec.template_length = tlen
    rec.set_tag("DB", barcode)
    if mate_cigar is not None:
        rec.set_tag("MC", mate_cigar)
    return rec


# Fragment 100-250: forward read 105-250 with a 5 bp 5' soft clip, reverse read 105-245 with a 5 bp 5' soft clip.
FWD = _read("p1", True, True, 105, "5S145M", 105, "140M5S", 140)
REV = _read("p1", False, False, 105, "140M5S", 105, "5S145M", -140)


def test_fragment_ends_from_either_mate():
    assert _fragment_ends(FWD) == (100, 250)
    assert _fragment_ends(REV) == (100, 250)


def test_fragment_ends_none_without_mc_tag():
    assert _fragment_ends(_read("p2", True, True, 105, "5S145M", 105, None, 140)) is None


def test_both_mates_share_a_label_with_fragment_key():
    assert _compute_read_label(FWD, FK) == _compute_read_label(REV, FK) == "AAA+CCC+f100-250"


def test_stock_label_uses_signed_template_length():
    assert _compute_read_label(FWD, STOCK) == "AAA+CCC+140"
    assert _compute_read_label(REV, STOCK) == "AAA+CCC+-140"


def test_differently_clipped_reads_of_one_molecule_share_a_label():
    # Same molecule, aligner keeps 3 more bases of the forward read and reports a different template length.
    other = _read("p3", True, True, 102, "2S148M", 105, "140M5S", 143)
    assert _compute_read_label(other, STOCK) != _compute_read_label(FWD, STOCK)
    assert _compute_read_label(other, FK) == _compute_read_label(FWD, FK)


def test_read_without_mc_tag_keeps_template_length_label():
    assert _compute_read_label(_read("p2", True, True, 105, "5S145M", 105, None, 140), FK) == "AAA+CCC+140"
