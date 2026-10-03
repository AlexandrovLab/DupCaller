"""trim fails loudly on bad input instead of writing corrupt FASTQ, detects
gzip per mate, and gives both mates the same name (no /1 /2)."""

import argparse
import gzip

import pytest

from DupCaller_sub.Trim import do_trim


def _fq(path, records, gz=False):
    text = "".join(f"{n}\n{s}\n+\n{q}\n" for n, s, q in records)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


def _run(tmp_path, r1, r2, pattern="NNNXX", pattern2=None):
    args = argparse.Namespace(
        fq=r1, fq2=r2, pattern=pattern, pattern2=pattern2, output=str(tmp_path / "out")
    )
    do_trim(args)
    return (
        (tmp_path / "out_1.fastq").read_text().split("\n"),
        (tmp_path / "out_2.fastq").read_text().split("\n"),
    )


def test_mate_suffixes_are_dropped_and_names_match(tmp_path):
    r1 = _fq(tmp_path / "a_1.fq", [("@frag1/1", "ACGTTGGGG", "IIIIIIIII")])
    r2 = _fq(tmp_path / "a_2.fq", [("@frag1/2 2:N:0", "TTTAAGGGG", "IIIIIIIII")])
    o1, o2 = _run(tmp_path, r1, r2)
    assert o1[0] == o2[0] == "@frag1_ACG+TTT DB:Z:ACG-TTT"
    assert o1[1] == "GGGG" and o2[1] == "GGGG"


def test_mixed_gzip_mates_are_detected_independently(tmp_path):
    r1 = _fq(tmp_path / "b_1.fq.gz", [("@f", "ACGTTGGGG", "IIIIIIIII")], gz=True)
    r2 = _fq(tmp_path / "b_2.fq", [("@f", "TTTAAGGGG", "IIIIIIIII")])
    o1, o2 = _run(tmp_path, r1, r2)
    assert o1[0] == o2[0] == "@f_ACG+TTT DB:Z:ACG-TTT"


def test_read_shorter_than_pattern_fails(tmp_path):
    r1 = _fq(tmp_path / "c_1.fq", [("@f", "AC", "II")])
    r2 = _fq(tmp_path / "c_2.fq", [("@f", "TTTAAGG", "IIIIIII")])
    with pytest.raises(ValueError, match="shorter than"):
        _run(tmp_path, r1, r2)


def test_mismatched_mate_names_fail(tmp_path):
    r1 = _fq(tmp_path / "d_1.fq", [("@f1", "ACGTTGG", "IIIIIII")])
    r2 = _fq(tmp_path / "d_2.fq", [("@f2", "TTTAAGG", "IIIIIII")])
    with pytest.raises(ValueError, match="Mate names differ"):
        _run(tmp_path, r1, r2)


def test_invalid_pattern_fails(tmp_path):
    r1 = _fq(tmp_path / "e_1.fq", [("@f", "ACGTTGG", "IIIIIII")])
    r2 = _fq(tmp_path / "e_2.fq", [("@f", "TTTAAGG", "IIIIIII")])
    with pytest.raises(ValueError, match="Invalid barcode pattern"):
        _run(tmp_path, r1, r2, pattern="NNA")


def test_read2_pattern_clips_read2_separately(tmp_path):
    # -p NNNXX on read 1, -p2 NNNXXXX on read 2: same 3-base barcodes,
    # read 2 loses 7 bases instead of 5.
    r1 = _fq(tmp_path / "f_1.fq", [("@f", "ACGTTGGGG", "IIIIIIIII")])
    r2 = _fq(tmp_path / "f_2.fq", [("@f", "TTTAACCGGGG", "IIIIIIIIIII")])
    o1, o2 = _run(tmp_path, r1, r2, pattern="NNNXX", pattern2="NNNXXXX")
    assert o1[0] == o2[0] == "@f_ACG+TTT DB:Z:ACG-TTT"
    assert o1[1] == "GGGG" and o2[1] == "GGGG"
    assert o1[3] == "IIII" and o2[3] == "IIII"


def test_read2_pattern_defaults_to_read1_pattern(tmp_path):
    r1 = _fq(tmp_path / "g_1.fq", [("@f", "ACGTTGGGG", "IIIIIIIII")])
    r2 = _fq(tmp_path / "g_2.fq", [("@f", "TTTAAGGGG", "IIIIIIIII")])
    assert _run(tmp_path, r1, r2, pattern2=None) == _run(
        tmp_path, r1, r2, pattern2="NNNXX"
    )


def test_read2_pattern_with_other_barcode_length_fails(tmp_path):
    r1 = _fq(tmp_path / "h_1.fq", [("@f", "ACGTTGGGG", "IIIIIIIII")])
    r2 = _fq(tmp_path / "h_2.fq", [("@f", "TTTAAGGGG", "IIIIIIIII")])
    with pytest.raises(ValueError, match="N \\(barcode\\) bases"):
        _run(tmp_path, r1, r2, pattern="NNNXX", pattern2="NNNNX")


def test_read2_shorter_than_its_pattern_fails(tmp_path):
    r1 = _fq(tmp_path / "i_1.fq", [("@f", "ACGTTGGGG", "IIIIIIIII")])
    r2 = _fq(tmp_path / "i_2.fq", [("@f", "TTTAAG", "IIIIII")])
    with pytest.raises(ValueError, match="shorter than"):
        _run(tmp_path, r1, r2, pattern="NNNXX", pattern2="NNNXXXX")


def test_invalid_read2_pattern_fails(tmp_path):
    r1 = _fq(tmp_path / "j_1.fq", [("@f", "ACGTTGG", "IIIIIII")])
    r2 = _fq(tmp_path / "j_2.fq", [("@f", "TTTAAGG", "IIIIIII")])
    with pytest.raises(ValueError, match="Invalid barcode pattern"):
        _run(tmp_path, r1, r2, pattern2="NNA")
