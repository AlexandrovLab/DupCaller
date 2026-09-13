"""Regression test for a thread-count-dependent miscount in funcs/call.py's
callBam family-batching loop, isolated from the full callBam pipeline (no
reference genome / error-matrix fixtures needed) by exercising the real
batching helpers directly: bamIterateMultipleRegionWithOverflow,
_index_rugged_mates, _drain_rugged_pool, _place_read_in_dict,
_compute_read_label -- and callBam's own "ignore if my mate precedes my
region start" boundary check, mirrored here to match callBam's real loop
exactly.

Background: PCR-duplicate copies of one original molecule can disagree
slightly on where their downstream mate maps (alignment noise) even though
their upstream (5') end -- and hence their family "label" (barcode pair +
signed template_length) -- agrees exactly. _index_rugged_mates reconciles
this: when an upstream batch closes, for every read whose recorded mate
position (next_reference_start) disagrees with the majority, it registers
that read's query_name so its actual downstream mate, when later
encountered, gets rerouted into the majority position's batch instead of
opening its own (spurious) batch. This makes both copies' downstream mates
land in ONE F1R2/F2R1-counted closure, correctly reflecting the true
2-copy family size.

The reroute only completes if the SAME worker's scan later reaches the
exact majority position to trigger the drain. Since `-p`/--threads
controls how many pieces the genome is split into (splitBamRegions),
whether a chunk boundary happens to fall between the minority position (C)
and the majority position (B) changes with thread count. When it does
(without the fix): the earlier worker's own scan ends before reaching B,
so the redirected minority mate's queued entry is silently discarded when
that worker exits; the majority mate, fetched independently by the next
worker (which never saw the upstream anchor), closes as its own singleton.
Net effect: a real 2-copy family ("2+0") is undercounted as a 1-copy
family ("1+0") -- an exact, integer, thread-count-dependent drop, exactly
matching real amp.tn.txt/dmg.tn.txt differences observed between -p 32 and
-p 64 reruns of the same production BAM.

The fix (bamIterateMultipleRegionWithOverflow in funcs/misc.py): the
EARLIER worker, instead of giving up at its own region_end, tracks how far
past its own end any read's *mate* reaches (computed from the MC tag, not
a random-access mate() lookup) and, if anything crosses the boundary, does
a second "overflow" fetch past its own end far enough to reach that mate
-- continuing the SAME in-memory batching state, so the existing
redirect+drain logic completes exactly as it would with no split at all.
To stay exactly-once (no double-count), the LATER worker's own normal scan
skips any read whose mate starts before its own region's start (callBam's
prev_boundary check) -- that read is the earlier worker's overflow pass's
responsibility, not the later worker's.
"""
import pysam
import pytest

from DupCaller_sub.funcs.misc import (
    bamIterateMultipleRegionWithOverflow,
    _compute_read_label,
    _drain_rugged_pool,
    _index_rugged_mates,
    _place_read_in_dict,
)

CHROM = "chr1"
CONTIG_LEN = 10000
U = 1000  # shared upstream (5') start for both copies' read1
B = 1300  # majority downstream mate position (copy1's read2)
C = 1290  # minority/rugged downstream mate position (copy2's read2)
READLEN = 50
MC_TAG = f"{READLEN}M"  # every read in this test uses a plain full-match CIGAR
BARCODE_PARAMS = {"barcodeTag": "DB", "barcodeNormalize": True, "barcodeSep": "-"}


def _read(header, name, is_read1, start, mate_start, tlen, barcode="AAA-TTT"):
    rec = pysam.AlignedSegment(header)
    rec.query_name = name
    rec.query_sequence = "A" * READLEN
    rec.query_qualities = pysam.qualitystring_to_array("I" * READLEN)
    rec.reference_id = 0
    rec.reference_start = start
    rec.mapping_quality = 60
    rec.cigar = [(0, READLEN)]
    rec.is_paired = True
    rec.is_proper_pair = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = (
        not is_read1
    )  # read1 forward (upstream), read2 reverse (downstream)
    rec.mate_is_reverse = is_read1
    rec.next_reference_id = 0
    rec.next_reference_start = mate_start
    rec.template_length = tlen
    rec.set_tag("NM", 0)
    rec.set_tag("AS", 100)
    rec.set_tag("XS", 0)
    rec.set_tag("DB", barcode)
    rec.set_tag("MC", MC_TAG)
    return rec


FILLER_POS = U + 5  # strictly between U and C -- see _build_bam docstring


def _build_bam(path, extra_pairs=()):
    """Two PCR-duplicate copies of one molecule, both read1s cleanly at U
    (same template_length -- same family label). copy1's read2 lands at
    the majority position B; copy2's read2 -- the "rugged" one -- lands a
    few bp off at C, but its read1 still (correctly) claims the majority
    template_length, which is exactly the disagreement
    _index_rugged_mates/next_reference_start detects and reconciles.

    Also includes one unrelated single-end "filler" read at FILLER_POS,
    strictly between U and C. Without it, C -- being the very next
    distinct position after U in a coordinate-sorted scan -- would be the
    very read whose own arrival triggers the U-batch's closure (and thus
    the _index_rugged_mates call that registers it for interception); a
    read can never be intercepted by a registration its own arrival is
    what causes, so it would slip through unintercepted regardless of any
    chunk boundary. In real (non-toy) data there is essentially always
    some other read between U and a rugged mate a few bp downstream; the
    filler reproduces that so this test isolates the boundary-specific
    failure mode instead of a degenerate same-batch artifact.

    extra_pairs: optional list of (name, start, mate_start, tlen, barcode)
    tuples for additional, unrelated read pairs (read1 at `start`, its
    downstream mate at `mate_start`) -- used to test that an unrelated
    family sitting in the overflow zone is neither dropped nor
    double-counted.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CONTIG_LEN}],
    }
    unsorted_path = str(path) + ".unsorted.bam"
    tlen = (B + READLEN) - U
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as bam:
        for name, mate_start in (("copy1", B), ("copy2", C)):
            r1 = _read(bam.header, name, True, U, mate_start, tlen)
            bam.write(r1)
        filler = _read(bam.header, "filler", True, FILLER_POS, FILLER_POS + 200, 250)
        bam.write(filler)
        # Downstream mates written separately so BAM sort order matches a
        # real coordinate-sorted file regardless of write order.
        r1b = _read(bam.header, "copy1", False, B, U, -tlen)
        r2b = _read(bam.header, "copy2", False, C, U, -tlen)
        bam.write(r1b)
        bam.write(r2b)
        for name, start, mate_start, extra_tlen, barcode in extra_pairs:
            r1 = _read(bam.header, name, True, start, mate_start, extra_tlen, barcode)
            bam.write(r1)
            r2 = _read(bam.header, name, False, mate_start, start, -extra_tlen, barcode)
            bam.write(r2)
    pysam.sort("-o", str(path), unsorted_path)
    pysam.index(str(path))


def _run_family_batching(bam_path, regions):
    """Reimplementation of callBam's batching loop (funcs/call.py
    ~1721-1782, ~2163-2176 in the fixed version), using the real
    misc.py/call.py helpers -- bamIterateMultipleRegionWithOverflow plus
    the same prev_boundary skip callBam applies -- without the
    surrounding reference-genome/error-matrix machinery callBam itself
    needs. Returns {label: [F1R2, F2R1]} accumulated over every batch
    closure (mirroring how duplex_read_num_dict[duplex_no] would be fed).
    """
    prev_boundary_chrom = None
    prev_boundary_pos = None
    if len(regions[0]) > 1 and regions[0][1] != 0:
        prev_boundary_chrom = regions[0][0]
        prev_boundary_pos = regions[0][1]

    currentReadDict = {}
    rugged_reads_index = {}
    rugged_reads_pool = {}
    rugged_pool_chrom = None
    currentStart = -1
    closures = []  # one (label, F1R2, F2R1) per batch closure, in order

    def close_batch(rd):
        for label, entry in rd.items():
            closures.append((label, entry["F1R2"], entry["F2R1"]))

    for rec, _region in bamIterateMultipleRegionWithOverflow(bam_path, regions, None):
        if (
            prev_boundary_chrom is not None
            and rec.reference_name == prev_boundary_chrom
            and rec.next_reference_start < prev_boundary_pos
        ):
            continue
        if rec.template_length < 0 and rec.query_name in rugged_reads_index:
            rugged_reads_pool[rugged_reads_index.pop(rec.query_name)].append(rec)
            continue
        start = rec.reference_start
        label = _compute_read_label(rec, BARCODE_PARAMS)
        chrom = rec.reference_name
        if chrom != rugged_pool_chrom:
            rugged_reads_pool.clear()
            rugged_reads_index.clear()
            rugged_pool_chrom = chrom
        if currentStart == -1:
            currentStart = start
            _drain_rugged_pool(
                chrom, currentStart, rugged_reads_pool, BARCODE_PARAMS, currentReadDict
            )
        if start == currentStart:
            _place_read_in_dict(rec, label, currentReadDict)
        else:
            _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
            close_batch(currentReadDict)
            currentReadDict = {}
            _place_read_in_dict(rec, label, currentReadDict)
            currentStart = start
            _drain_rugged_pool(
                chrom, currentStart, rugged_reads_pool, BARCODE_PARAMS, currentReadDict
            )

    _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
    close_batch(currentReadDict)
    return closures


DOWNSTREAM_LABEL = "AAA+TTT+-350"


def test_rugged_mate_merges_into_one_closure_without_boundary(tmp_path):
    """No chunk split: copy1's and copy2's downstream mates (majority B,
    rugged C) merge into ONE closure via the redirect+drain, correctly
    reflecting both copies (F1R2=2) as a single family observation."""
    bam_path = tmp_path / "family.bam"
    _build_bam(bam_path)

    closures = _run_family_batching(str(bam_path), [(CHROM, 0, CONTIG_LEN)])
    downstream_closures = [c for c in closures if c[0] == DOWNSTREAM_LABEL]

    assert downstream_closures == [(DOWNSTREAM_LABEL, 2, 0)], (
        "expected the two copies' downstream mates to merge into a single "
        f"F1R2=2 closure; got {downstream_closures}"
    )


@pytest.mark.parametrize("boundary", [(C + B) // 2, B])
def test_chunk_boundary_between_rugged_and_majority_position_is_fixed(
    tmp_path, boundary
):
    """A chunk boundary strictly between C (rugged) and B (majority) used
    to drop copy2's downstream mate entirely (see this test's previous
    version / the module docstring for the bug). With
    bamIterateMultipleRegionWithOverflow: the earlier worker (holding U/C)
    notices copy1's mate (at B) crosses its own region_end (computed from
    the MC tag, no fixed buffer needed), does a second "overflow" fetch
    reaching B, and completes the SAME redirect+drain merge it would have
    done with no split at all. The later worker's own normal scan skips
    the read at B outright (its mate at U precedes the later worker's
    region start), so the family is counted exactly once, by the earlier
    worker, matching the no-split baseline.
    """
    bam_path = tmp_path / "family.bam"
    _build_bam(bam_path)

    assert C < boundary <= B

    worker1 = _run_family_batching(str(bam_path), [(CHROM, 0, boundary)])
    worker2 = _run_family_batching(str(bam_path), [(CHROM, boundary, CONTIG_LEN)])
    all_closures = worker1 + worker2
    downstream_closures = [c for c in all_closures if c[0] == DOWNSTREAM_LABEL]

    assert downstream_closures == [(DOWNSTREAM_LABEL, 2, 0)], (
        "expected the earlier worker's overflow pass to fully reconstruct "
        f"the family (F1R2=2), matching the no-split baseline; got {downstream_closures}"
    )
    # And it must be worker1 (the one with the overflow pass) that owns
    # this closure, not worker2 -- otherwise this "fix" would just be
    # trading the drop for a double-count spread across both workers.
    worker1_downstream = [c for c in worker1 if c[0] == DOWNSTREAM_LABEL]
    worker2_downstream = [c for c in worker2 if c[0] == DOWNSTREAM_LABEL]
    assert worker1_downstream == [(DOWNSTREAM_LABEL, 2, 0)]
    assert worker2_downstream == []


def test_unrelated_read_in_overflow_zone_is_not_dropped_or_double_counted(tmp_path):
    """A second, entirely unrelated family sitting inside the earlier
    worker's overflow window (but not itself crossing the boundary) must
    be (a) excluded from the earlier worker's overflow pass -- its mate
    doesn't precede the boundary, so bamIterateMultipleRegionWithOverflow
    skips it there -- and (b) picked up normally by the later worker's own
    scan instead, exactly once, at F1R2=1. This is the complementary
    check to the rugged-mate fix: the overflow pass must be narrowly
    scoped to genuine boundary-crossing families, not just "everything in
    range"."""
    bam_path = tmp_path / "family_with_unrelated.bam"
    boundary = (C + B) // 2
    UNRELATED_START = boundary + 15  # inside the overflow window (boundary, last_end]
    UNRELATED_MATE = 1500  # downstream, does NOT precede the boundary
    unrelated_label = "CCC+GGG+" + str((UNRELATED_MATE + READLEN) - UNRELATED_START)
    _build_bam(
        bam_path,
        extra_pairs=[
            (
                "unrelated",
                UNRELATED_START,
                UNRELATED_MATE,
                (UNRELATED_MATE + READLEN) - UNRELATED_START,
                "CCC-GGG",
            )
        ],
    )
    assert C < boundary < B < UNRELATED_START < UNRELATED_MATE

    worker1 = _run_family_batching(str(bam_path), [(CHROM, 0, boundary)])
    worker2 = _run_family_batching(str(bam_path), [(CHROM, boundary, CONTIG_LEN)])

    worker1_unrelated = [c for c in worker1 if c[0] == unrelated_label]
    worker2_unrelated = [c for c in worker2 if c[0] == unrelated_label]
    assert (
        worker1_unrelated == []
    ), f"unrelated family must not be pulled into worker1's overflow pass; got {worker1_unrelated}"
    assert worker2_unrelated == [
        (unrelated_label, 1, 0)
    ], f"unrelated family must be counted exactly once, by worker2; got {worker2_unrelated}"

    # And the rugged-mate fix itself must still hold with this extra data
    # sitting nearby.
    all_downstream = [c for c in worker1 + worker2 if c[0] == DOWNSTREAM_LABEL]
    assert all_downstream == [(DOWNSTREAM_LABEL, 2, 0)]


@pytest.mark.parametrize("boundary", [1100, (C + B) // 2, B])
def test_previous_coverage_extent_includes_mates_across_a_gap(tmp_path, boundary):
    from DupCaller_sub.funcs.misc import previous_region_coverage_start

    path = tmp_path / "extent.bam"
    _build_bam(path)
    assert previous_region_coverage_start(str(path), [(CHROM, 0, boundary)]) == (
        CHROM,
        B + READLEN,
    )


def test_prerouted_coverage_merges_without_rewriting_main_files(tmp_path, monkeypatch):
    from collections import Counter
    from Bio import bgzf
    from DupCaller_sub import Caller
    from DupCaller_sub.funcs.misc import (
        get_bed_file_for_position,
        previous_region_coverage_start,
    )

    path = tmp_path / "routing.bam"
    boundary = 1100  # no alignment overlaps the boundary; mates are farther away
    _build_bam(path, extra_pairs=[("other", 1320, 1500, 230, "CCC-GGG")])
    regions = [[(CHROM, 0, boundary)], [(CHROM, boundary, CONTIG_LEN)]]
    extent = previous_region_coverage_start(str(path), regions[0])
    assert extent == (CHROM, B + READLEN)
    expected = Counter()
    for worker in range(2):
        counts = Counter()
        for rec, _ in bamIterateMultipleRegionWithOverflow(
            str(path), regions[worker], None
        ):
            if worker and rec.next_reference_start < boundary:
                continue
            counts.update(range(rec.reference_start, rec.reference_end))
        expected.update(counts)
        prefix = str(tmp_path / f"ROUTE_{worker}_coverage")
        with bgzf.open(prefix + ".bed.gz", "wt") as main, bgzf.open(
            prefix + "_prev_region.tmp.bed.gz", "wt"
        ) as prev, bgzf.open(prefix + "_next_region.tmp.bed.gz", "wt") as next_file:
            for pos, depth in sorted(counts.items()):
                target = get_bed_file_for_position(
                    pos,
                    CHROM,
                    CHROM,
                    extent[1] if worker else 0,
                    CHROM,
                    CONTIG_LEN if worker else boundary,
                    main,
                    prev,
                    next_file,
                )
                target.write(f"{CHROM}\t{pos}\t{pos + 1}\t{depth}\n")

    def unexpected_split(*args):
        raise AssertionError("Pre-routed coverage must not rewrite main BED files")

    monkeypatch.setattr(Caller, "_split_bed_by_start", unexpected_split)
    Caller.merge_and_combine_coverage_files("ROUTE", str(tmp_path), 2)
    with bgzf.open(str(tmp_path / "ROUTE_coverage.bed.gz"), "rt") as merged:
        rows = [line.strip().split("\t") for line in merged]
    assert [int(row[1]) for row in rows] == sorted(expected)
    assert [float(row[3]) for row in rows] == [
        expected[pos] for pos in sorted(expected)
    ]
    assert (tmp_path / "ROUTE_coverage.bed.gz.tbi").exists()


def test_coverage_boundary_position_goes_to_next_file():
    from DupCaller_sub.funcs.misc import get_bed_file_for_position

    assert (
        get_bed_file_for_position(
            1100, CHROM, CHROM, 0, CHROM, 1100, "main", "prev", "next"
        )
        == "next"
    )
