"""--skipCoveragePass round-1 deferred depth: one tumor extraction per
distinct founding-family LR at or above --minLR, each exempting exactly the
families at or above that LR (the founding families round 2 would emit at
the final threshold); indels a single extraction (cut -inf)."""

import numpy as np

from DupCaller_sub.funcs.call import _early_deferred_depth
from DupCaller_sub.funcs.misc import str_profile_is_current


def _rec(lr, tag, masked=True):
    return {
        "chrom": "chr1",
        "pos": 100,
        "ref": "C",
        "alt": "A",
        "filter": "masked" if masked else "PASS",
        "infos": {
            "LR": lr,
            "NOISEM": 1,
            "SNPM": 0,
            "TAG1": tag,
            "TAG2": tag,
            "TL": 300,
            "FS": [50],
        },
    }


def test_snv_depth_per_family_lr_level():
    recs = [_rec(8.0, "AAA"), _rec(6.0, "CCC"), _rec(2.0, "GGG")]
    seen = []

    def extract(bam, keys, params, minbq=None, call_barcodes=None):
        if call_barcodes is None:  # normal
            return {k: (0, 30, 0, 30) for k in keys}
        tags = sorted(t[0] for t in call_barcodes[keys[0]])
        seen.append(tags)
        return {k: (len(tags), 20, 0, 40) for k in keys}

    params = {"minLR": 5.0, "minBq": 18, "maxAF": 1, "minNdepth": 10, "normalVAF": 0.01}
    _early_deferred_depth(recs, "T", ["N"], params, extract, indel=False)
    # Levels: LR >= 8 (A) and LR >= 6 (A, C); the LR-2 family never.
    assert seen == [["AAA"], ["AAA", "CCC"]]
    entries = recs[0]["deferred_depth"]
    assert [cut for cut, _ in entries] == [8.0, 6.0]
    assert entries[0][1][0][0][0] == 1 and entries[1][1][0][0][0] == 2


def test_indel_single_extraction_and_normal_insertion_drop():
    recs = [_rec(8.0, "AAA"), _rec(3.0, "CCC")]

    def extract(bam, keys, params, minbq=None, call_barcodes=None):
        if call_barcodes is None:
            return {k: (0, 30, 1, 31) for k in keys}  # a normal read with the indel
        return {k: (2, 20, 0, 40) for k in keys}

    params = {"minLR": 5.0, "minBq": 18, "maxAF": 1, "minNdepth": 10, "normalVAF": 0.5}
    _early_deferred_depth(recs, "T", ["N"], params, extract, indel=True)
    assert all("deferred_depth" not in r for r in recs)  # round 2 drops ni > 0
    _early_deferred_depth(recs, "T", [], params, extract, indel=True)
    assert [cut for cut, _ in recs[0]["deferred_depth"]] == [-np.inf]


def test_str_profile_layout_check(tmp_path):
    new = tmp_path / "new.str.txt"
    new.write_text(
        "\t" + "\t".join(f"c{i}" for i in range(15)) + "\n0-1" + "\t0" * 15 + "\n"
    )
    old = tmp_path / "old.str.txt"
    old.write_text(
        "\t" + "\t".join(str(d) for d in range(-5, 6)) + "\n0-1" + "\t0" * 11 + "\n"
    )
    assert str_profile_is_current(new) and not str_profile_is_current(old)
