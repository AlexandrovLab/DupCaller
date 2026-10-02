"""SBS and DBS are mutually exclusive: the two bases of a PASS DBS leave the
SBS set as FILTER dbs_member, matched by read family, not by position."""

from DupCaller_sub.Caller import _relabel_dbs_members


def _rec(pos, tag1="AAA", filt="PASS"):
    return {
        "chrom": "chr1",
        "pos": pos,
        "filter": filt,
        "infos": {"TAG1": tag1, "TAG2": "CCC", "SP": 100, "TL": 150},
    }


def test_dbs_bases_become_dbs_member():
    sbs = [_rec(10), _rec(11), _rec(20)]
    assert _relabel_dbs_members(sbs, [_rec(10)]) == 2
    assert [m["filter"] for m in sbs] == ["dbs_member", "dbs_member", "PASS"]


def test_other_family_at_dbs_position_stays_pass():
    sbs = [_rec(10), _rec(11), _rec(10, tag1="GGG")]
    _relabel_dbs_members(sbs, [_rec(10)])
    assert sbs[2]["filter"] == "PASS"


def test_failed_dbs_leaves_its_sbs_alone():
    sbs = [_rec(10), _rec(11)]
    assert _relabel_dbs_members(sbs, [_rec(10, filt="underpowered")]) == 0
    assert [m["filter"] for m in sbs] == ["PASS", "PASS"]


def test_non_pass_sbs_keeps_its_filter():
    sbs = [_rec(10, filt="strand_independence"), _rec(11)]
    _relabel_dbs_members(sbs, [_rec(10)])
    assert sbs[0]["filter"] == "strand_independence"
