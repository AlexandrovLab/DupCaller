"""trim's B pattern letter: a barcode from --barcode-list, of any length,
matched by walking all listed barcodes at once and dropping each once its
mismatches exceed --max-mismatch; fewest mismatches wins, then shortest."""

import argparse

import pytest

import random

from DupCaller_sub.Trim import (
    AMBIGUOUS,
    NO_MATCH,
    _validate_pattern,
    build_barcode_index,
    do_trim,
    load_barcode_list,
    make_barcode_matcher,
    match_barcode,
    match_barcode_indexed,
)

BCS = ["GGCACCGAAAA", "CTCGGCGATAAA", "CTGAGCTCGTTTT"]  # 11, 12, 13 nt
ME = "AGATGTGTATAAGAGACAG"


def test_exact_match_of_each_length():
    for bc in BCS:
        assert match_barcode(bc + ME + "ACGT", 0, BCS, 1) == (bc, len(bc))


def test_one_mismatch_allowed_by_default_two_not():
    read = "GGCACCGTAAA" + ME  # 1 mismatch vs GGCACCGAAAA
    assert match_barcode(read, 0, BCS, 1) == ("GGCACCGAAAA", 11)
    read2 = "GGCACCGTATA" + ME  # 2 mismatches
    assert match_barcode(read2, 0, BCS, 1) == (NO_MATCH, 0)
    assert match_barcode(read2, 0, BCS, 2) == ("GGCACCGAAAA", 11)


def test_fewest_mismatches_beats_shorter():
    bcs = ["ACGTAC", "ACGTTCGG"]
    # 1 mismatch vs the 6-mer, 0 vs the 8-mer
    assert match_barcode("ACGTTCGGAAAA", 0, bcs, 1) == ("ACGTTCGG", 8)


def test_tie_on_mismatches_picks_shortest():
    bcs = ["ACGTACGG", "ACGTAC"]
    assert match_barcode("ACGTACGGTTTT", 0, bcs, 1) == ("ACGTAC", 6)


def test_same_length_tie_is_ambiguous():
    bcs = ["AAAAAA", "AAAAAT"]
    assert match_barcode("AAAAAGCC", 0, bcs, 1) == (AMBIGUOUS, 0)


def test_barcode_longer_than_remaining_read_does_not_finish():
    assert match_barcode("GGCACCG", 0, BCS, 1) == (NO_MATCH, 0)


def test_match_starts_after_a_fixed_prefix():
    assert match_barcode("NN" + "CTCGGCGATAAA" + ME, 2, BCS, 1) == ("CTCGGCGATAAA", 12)


def test_pattern_validation():
    _validate_pattern("B" + "X" * 19, BCS)
    _validate_pattern("NNXB", BCS)
    with pytest.raises(ValueError, match="needs --barcode-list"):
        _validate_pattern("BXX")
    with pytest.raises(ValueError, match="has no B"):
        _validate_pattern("NNNXX", BCS)
    with pytest.raises(ValueError, match="at most one B"):
        _validate_pattern("BXB", BCS)
    with pytest.raises(ValueError, match="Invalid barcode pattern"):
        _validate_pattern("BQ", BCS)


def test_barcode_list_file_checks(tmp_path):
    p = tmp_path / "bc.txt"
    p.write_text("# comment\nggcaccgaaaa\n\nCTCGGCGATAAA\n")
    assert load_barcode_list(str(p)) == ["GGCACCGAAAA", "CTCGGCGATAAA"]
    p.write_text("ACGT\nACGT\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_barcode_list(str(p))
    p.write_text("ACGN\n")
    with pytest.raises(ValueError, match="may only contain"):
        load_barcode_list(str(p))


def _fq(path, records):
    path.write_text("".join(f"{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in records))
    return str(path)


def test_do_trim_with_barcode_list(tmp_path, capsys):
    lst = tmp_path / "bc.txt"
    lst.write_text("\n".join(BCS) + "\n")
    ins1, ins2 = "ACGTACGTAC", "TTGGCCAATT"
    r1 = _fq(
        tmp_path / "r1.fq",
        [
            ("@p1", "GGCACCGAAAA" + ME + ins1),  # exact
            ("@p2", "CTCGGCGTTAAA" + ME + ins1),  # 1 mismatch, corrected
            ("@p3", "TTTTTTTTTTTTT" + ME + ins1),  # no match -> dropped
        ],
    )
    r2 = _fq(
        tmp_path / "r2.fq",
        [
            ("@p1", "CTGAGCTCGTTTT" + ME + ins2),
            ("@p2", "GGCACCGAAAA" + ME + ins2),
            ("@p3", "GGCACCGAAAA" + ME + ins2),
        ],
    )
    do_trim(
        argparse.Namespace(
            fq=r1,
            fq2=r2,
            pattern="B" + "X" * 19,
            output=str(tmp_path / "out"),
            barcode_list=str(lst),
            max_mismatch=1,
        )
    )
    o1 = (tmp_path / "out_1.fastq").read_text().split("\n")
    o2 = (tmp_path / "out_2.fastq").read_text().split("\n")
    assert (
        o1[0] == o2[0] == "@p1_GGCACCGAAAA+CTGAGCTCGTTTT DB:Z:GGCACCGAAAA-CTGAGCTCGTTTT"
    )
    assert o1[1] == ins1 and o2[1] == ins2
    assert o1[4] == "@p2_CTCGGCGATAAA+GGCACCGAAAA DB:Z:CTCGGCGATAAA-GGCACCGAAAA"
    assert len(o1) == 9  # 2 records + trailing empty
    assert "dropped 1 (read 1 barcode not in list)" in capsys.readouterr().out


def test_negative_max_mismatch_fails(tmp_path):
    lst = tmp_path / "bc.txt"
    lst.write_text("ACGT\n")
    r = _fq(tmp_path / "r.fq", [("@a", "ACGTAAAA")])
    with pytest.raises(ValueError, match="max-mismatch"):
        do_trim(
            argparse.Namespace(
                fq=r,
                fq2=r,
                pattern="B",
                output=str(tmp_path / "o"),
                barcode_list=str(lst),
                max_mismatch=-1,
            )
        )


def _random_reads_near(bcs, rng, n):
    """Reads starting with a listed barcode carrying 0-3 random errors, a
    random stretch, or a truncated barcode -- plus random tails."""
    reads = []
    for _ in range(n):
        kind = rng.random()
        if kind < 0.8:
            b = list(rng.choice(bcs))
            for _ in range(rng.choice([0, 0, 1, 1, 2, 3])):
                p = rng.randrange(len(b))
                b[p] = rng.choice("ACGTN")
            read = "".join(b) + "".join(
                rng.choice("ACGT") for _ in range(rng.randrange(0, 8))
            )
        elif kind < 0.9:
            read = "".join(rng.choice("ACGTN") for _ in range(rng.randrange(0, 16)))
        else:
            read = rng.choice(bcs)[: rng.randrange(1, 10)]
        reads.append(read)
    return reads


@pytest.mark.parametrize("max_mm", [0, 1])
@pytest.mark.parametrize(
    "bcs",
    [
        BCS,
        [
            "AAAAAA",
            "AAAAAT",
            "AAAATT",
            "AAAAAAAA",
            "AAAAAATT",
            "CCCCCC",
        ],  # ties, prefixes
        ["ACGTAC", "ACGTACGG", "ACGTTCGG", "TTTTTTTTTTTT", "ACG"],
    ],
)
def test_indexed_matches_walk(bcs, max_mm):
    """The lookup-table path must give exactly the walk's answer."""
    rng = random.Random(len(bcs) * 10 + max_mm)
    index = build_barcode_index(bcs, max_mm)
    for read in _random_reads_near(bcs, rng, 4000):
        for start in (0, 1):
            assert match_barcode_indexed(read, start, index) == match_barcode(
                read, start, bcs, max_mm
            ), (read, start)


def test_matcher_uses_walk_above_one_mismatch():
    read = "GGCACCGTATA" + ME  # 2 mismatches
    assert make_barcode_matcher(BCS, 1)(read, 0) == (NO_MATCH, 0)
    assert make_barcode_matcher(BCS, 2)(read, 0) == ("GGCACCGAAAA", 11)
