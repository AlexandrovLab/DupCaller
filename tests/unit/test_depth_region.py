"""extractDepthRegion must return exactly [start, end) depth, and
prepare_reference_mats must keep a normal-depth failure masked when the
tumor maxAF depth requirement also applies."""
import numpy as np
import pysam

from DupCaller_sub.funcs import call
from DupCaller_sub.funcs.depth import extractDepthRegion

CHROM = "chr1"


def _bam(path, depth_at):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": 1000}],
    }
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        n = 0
        for pos, depth in sorted(depth_at.items()):
            for _ in range(depth):
                rec = pysam.AlignedSegment(bam.header)
                rec.query_name = f"r{n}"
                n += 1
                rec.query_sequence = "A"
                rec.query_qualities = [40]
                rec.reference_id = 0
                rec.reference_start = pos
                rec.mapping_quality = 60
                rec.cigar = [(0, 1)]
                bam.write(rec)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))


def test_extract_depth_region_is_half_open(tmp_path):
    bam = tmp_path / "t.bam"
    # 0-based: 99 just before, 100 first, 150 inside, 199 left empty so a
    # wrapped write from 99 would show up there, 200 just after.
    _bam(bam, {99: 3, 100: 1, 150: 2, 200: 5})
    depth = extractDepthRegion(
        str(bam), CHROM, 100, 200, {"mapq": 0, "maxDepth": 1000, "minBq": 30}
    )
    expected = np.zeros(100)
    expected[0] = 1
    expected[50] = 2
    np.testing.assert_array_equal(depth, expected)


def test_extract_depth_region_allele_counts_half_open(tmp_path):
    bam = tmp_path / "t.bam"
    _bam(bam, {99: 3, 100: 1})
    depth, acgt = extractDepthRegion(
        str(bam),
        CHROM,
        100,
        110,
        {"mapq": 0, "maxDepth": 1000, "minBq": 30},
        count_alleles=True,
    )
    assert depth[-1] == 0 and acgt[-1].sum() == 0
    assert depth[0] == 1


def test_maxaf_mask_does_not_overwrite_normal_mask(monkeypatch):
    # locus 0: fails normal depth only; 1: fails tumor depth only; 2: passes.
    normal = np.array([1.0, 10.0, 10.0])
    tumor = np.array([10.0, 1.0, 10.0])
    monkeypatch.setattr(
        call,
        "extractDepthRegion",
        lambda bam, *a, **k: normal if bam == "normal" else tumor,
    )
    monkeypatch.setattr(
        call, "prepareAlignMask", lambda *a, **k: np.zeros(3, dtype=bool)
    )
    params = {
        "germline_cutoff": 0.5,
        "trinuc2num_dict": {},
        "isLearn": False,
        "minNdepth": 5,
        "maxAF": 0.5,  # -> tumor min depth 2
    }
    masks = call.prepare_reference_mats(
        CHROM,
        0,
        3,
        np.zeros(3, dtype=int),
        np.zeros(3, dtype=int),
        None,
        None,
        None,
        None,
        ["normal"],
        "tumor",
        params,
    )
    n_cov_mask = masks[3]
    np.testing.assert_array_equal(n_cov_mask, [True, True, False])
