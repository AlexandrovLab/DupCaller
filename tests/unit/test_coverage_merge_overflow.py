"""Regression test for Caller.py's merge_and_combine_coverage_files: an
overflow-extended `next_region` file (bamIterateMultipleRegionWithOverflow,
funcs/misc.py, lets a worker's coverage flush reach arbitrarily far past its
own nominal region_end to correctly attribute a boundary-crossing duplex
family) used to break the plain chunk-order concatenation this function
relied on, since the following worker's own main coverage file could start
at a position *lower* than wherever the overflow reached -- producing an
unsorted final bed file that crashed `tabix -f -p bed` indexing (observed on
16/17 real production PD samples' full-genome runs: "Unsorted positions on
sequence #2: 15864109 followed by 15863937").

The fix: peel any leading rows of the following worker's main file that fall
at or below the overflow reach into the same sum-on-matching-position merge
already used for the small next_region/prev_region overlap, instead of
assuming a plain concatenation is already sorted. Positions appearing in
both the overflow-extended next_region file and the peeled main-file rows
represent genuinely different, non-overlapping families (the complementary
"later worker skips reads whose mate precedes its own region start" rule
only excludes specific reads, not everything at that position) -- so they
must be *summed*, not deduplicated or last-write-wins overwritten.
"""
import os

from Bio import bgzf

from DupCaller_sub.Caller import merge_and_combine_coverage_files


def _write_bed(path, rows):
    with bgzf.open(path, "wt") as f:
        for chrom, start, end, vals in rows:
            f.write(
                "\t".join([chrom, str(start), str(end)] + [str(v) for v in vals]) + "\n"
            )


def _read_bed(path):
    rows = []
    with bgzf.open(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            rows.append(
                (parts[0], int(parts[1]), int(parts[2]), [float(v) for v in parts[3:]])
            )
    return rows


def test_overflow_extended_overlap_is_sorted_and_summed(tmp_path):
    sample = "TEST"
    tmp = str(tmp_path)

    # Worker 0: main file within its own nominal range.
    _write_bed(
        os.path.join(tmp, f"{sample}_0_coverage.bed.gz"),
        [("chr1", 100, 101, [5, 0, 0, 0]), ("chr1", 200, 201, [3, 0, 0, 0])],
    )
    # Worker 0's overflow-extended next_region file reaches past worker 1's
    # own region_start (300) all the way to 350.
    _write_bed(
        os.path.join(tmp, f"{sample}_0_coverage_next_region.tmp.bed.gz"),
        [("chr1", 300, 301, [2, 0, 0, 0]), ("chr1", 350, 351, [1, 0, 0, 0])],
    )
    _write_bed(os.path.join(tmp, f"{sample}_1_coverage_prev_region.tmp.bed.gz"), [])
    # Worker 1's own main file starts at its nominal region_start (300),
    # including a row at the SAME position (300) as the overflow row --
    # from an unrelated family -- to verify summation, plus a row (320)
    # inside the overflow's reach with no matching next_region row (must
    # survive untouched), plus a row (400) safely past the reach.
    _write_bed(
        os.path.join(tmp, f"{sample}_1_coverage.bed.gz"),
        [
            ("chr1", 300, 301, [4, 0, 0, 0]),
            ("chr1", 320, 321, [7, 0, 0, 0]),
            ("chr1", 400, 401, [9, 0, 0, 0]),
        ],
    )

    merge_and_combine_coverage_files(sample, tmp, nprocess=2)

    final = os.path.join(tmp, f"{sample}_coverage.bed.gz")
    assert os.path.exists(final)
    assert os.path.exists(final + ".tbi"), "tabix indexing must not fail"

    rows = _read_bed(final)
    positions = [r[1] for r in rows]
    assert positions == sorted(positions), f"output not sorted: {positions}"

    by_pos = {r[1]: r[3][0] for r in rows}
    assert by_pos[100] == 5.0
    assert by_pos[200] == 3.0
    # Overflow row (2) + worker 1's own unrelated-family row (4) at the same
    # position must be summed, not deduplicated/overwritten.
    assert by_pos[300] == 6.0
    # Inside the overflow's reach but no matching next_region row -- must
    # survive the peel untouched.
    assert by_pos[320] == 7.0
    assert by_pos[350] == 1.0
    assert by_pos[400] == 9.0


def test_no_overflow_matches_prior_plain_concatenation_behavior(tmp_path):
    """The common case (no overflow ever occurred, next_region/prev_region
    empty) must behave exactly as the original plain-concatenation code
    did: each worker's main file untouched and simply stitched together."""
    sample = "BASE"
    tmp = str(tmp_path)

    _write_bed(
        os.path.join(tmp, f"{sample}_0_coverage.bed.gz"),
        [("chr1", 100, 101, [5, 0, 0, 0]), ("chr1", 200, 201, [3, 0, 0, 0])],
    )
    _write_bed(os.path.join(tmp, f"{sample}_0_coverage_next_region.tmp.bed.gz"), [])
    _write_bed(os.path.join(tmp, f"{sample}_1_coverage_prev_region.tmp.bed.gz"), [])
    _write_bed(
        os.path.join(tmp, f"{sample}_1_coverage.bed.gz"),
        [("chr1", 300, 301, [4, 0, 0, 0]), ("chr1", 400, 401, [9, 0, 0, 0])],
    )

    merge_and_combine_coverage_files(sample, tmp, nprocess=2)

    final = os.path.join(tmp, f"{sample}_coverage.bed.gz")
    assert os.path.exists(final) and os.path.exists(final + ".tbi")
    rows = _read_bed(final)
    assert [r[1] for r in rows] == [100, 200, 300, 400]
    assert [r[3][0] for r in rows] == [5.0, 3.0, 4.0, 9.0]


def test_overflow_merge_preserves_next_chromosome_order(tmp_path):
    """A worker spanning chromosomes must not peel low positions on the next one."""
    sample = "MULTI"
    files = {
        "0_coverage.bed.gz": [("chr1", 100, 101, [5])],
        "0_coverage_next_region.tmp.bed.gz": [("chr1", 350, 351, [2])],
        "1_coverage_prev_region.tmp.bed.gz": [],
        "1_coverage.bed.gz": [
            ("chr1", 350, 351, [4]),
            ("chr1", 400, 401, [7]),
            ("chr2", 100, 101, [3]),
            ("chr2", 400, 401, [9]),
        ],
    }
    for suffix, rows in files.items():
        _write_bed(str(tmp_path / f"{sample}_{suffix}"), rows)
    merge_and_combine_coverage_files(sample, str(tmp_path), nprocess=2)
    final = str(tmp_path / f"{sample}_coverage.bed.gz")
    assert os.path.exists(final + ".tbi")
    assert _read_bed(final) == [
        ("chr1", 100, 101, [5]),
        ("chr1", 350, 351, [6]),
        ("chr1", 400, 401, [7]),
        ("chr2", 100, 101, [3]),
        ("chr2", 400, 401, [9]),
    ]


def test_compressed_concatenation_preserves_blocks_and_supports_tabix(tmp_path):
    from DupCaller_sub.Caller import _concatenate_bgzf_files
    import pysam

    first, second, empty = (
        tmp_path / f"{name}.gz" for name in ("first", "second", "empty")
    )
    _write_bed(str(first), [("chr1", i, i + 1, [i % 7]) for i in range(10000)])
    _write_bed(str(second), [("chr2", i, i + 1, [i % 11]) for i in range(10000)])
    _write_bed(str(empty), [])
    output = tmp_path / "merged.bed.gz"
    _concatenate_bgzf_files([str(first), str(empty), str(second)], str(output))
    # All input compressed blocks survive byte-for-byte; only intermediate
    # 28-byte BGZF EOF markers are removed.
    assert output.read_bytes() == first.read_bytes()[:-28] + second.read_bytes()
    pysam.tabix_index(str(output), preset="bed", force=True)
    with pysam.TabixFile(str(output)) as indexed:
        assert len(list(indexed.fetch("chr1"))) == 10000
        assert len(list(indexed.fetch("chr2"))) == 10000
