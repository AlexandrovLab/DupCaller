"""Regression test for funcs/misc.py's splitBamRegions cut-site math.

Merging origin/main's fractional within-window interpolation (fixes
duplicate/zero-length cut sites when one window holds a disproportionate
share of reads) with dev's existing end-of-contig clamp (fixes cut sites
overshooting the contig when -p requests more chunks than there's read
data to fill) introduced a new bug: the interpolation indexes
total_reads_by_windows[cut_inds], but an oversaturated -p can pin
cut_inds one past the end of that array, raising IndexError before the
existing clamp ever gets a chance to run.
"""
import pysam

from DupCaller_sub.funcs.misc import splitBamRegions

CHROM1 = "chr1"
CHROM1_LEN = 500
CHROM2 = "chr2"
CHROM2_LEN = 300


def _read(header, name, chrom_tid, pos):
    rec = pysam.AlignedSegment(header)
    rec.query_name = name
    rec.query_sequence = "A" * 20
    rec.query_qualities = pysam.qualitystring_to_array("I" * 20)
    rec.reference_id = chrom_tid
    rec.reference_start = pos
    rec.mapping_quality = 60
    rec.cigar = [(0, 20)]
    return rec


def _build_bam(path):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": CHROM1, "LN": CHROM1_LEN},
            {"SN": CHROM2, "LN": CHROM2_LEN},
        ],
    }
    unsorted_path = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as bam:
        # A handful of reads only -- with a large requested thread count
        # (num) below, this makes the trailing chunk-boundary targets
        # exceed the total read count, the oversaturation case.
        for i, pos in enumerate([10, 200, 400]):
            bam.write(_read(bam.header, f"chr1_{i}", 0, pos))
        for i, pos in enumerate([50, 250]):
            bam.write(_read(bam.header, f"chr2_{i}", 1, pos))
    pysam.sort("-o", str(path), unsorted_path)
    pysam.index(str(path))


def test_oversaturated_thread_count_does_not_crash_and_stays_in_bounds(tmp_path):
    bam_path = tmp_path / "test.bam"
    _build_bam(bam_path)

    cutSite, chunkSize, contigs_sorted = splitBamRegions(
        [str(bam_path)], num=20, contigs=[CHROM1, CHROM2], step=50
    )

    contig_lens = {0: CHROM1_LEN, 1: CHROM2_LEN}
    for contig_idx, pos in cutSite:
        assert pos <= contig_lens[contig_idx]

    # Non-decreasing contig index across cut sites (regions are consumed
    # in contig order downstream).
    contig_indices = [c for c, _ in cutSite]
    assert contig_indices == sorted(contig_indices)
