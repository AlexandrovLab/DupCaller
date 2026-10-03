"""classify_indel_channel must give the same ID83 channel as
SigProfilerMatrixGenerator (1.3.6) for left-aligned indels. Expected labels
below are SPMG's own output for these contexts (checked 2026-10-02 against
the installed package, and on 46,370 planted chr22 indels).

Each case gives the reference from the base after the anchor onward
(`after_anchor`) and the event; the anchor is chosen so the event stays
left-aligned. The `anno` column (hp.h5/str.h5 values, which the classifier
no longer reads) is kept to show that the str.h5 annotation is irrelevant.
"""

import pytest

from DupCaller_sub.Estimate import INDEL83_TO_SIGPROFILER_LABELS
from DupCaller_sub.funcs.misc import classify_indel_channel

CASES = [
    # name, after_anchor, kind, payload, anno, SPMG label
    # 1bp homopolymer
    ("del C from C1", "CTTGA", "del", 1, (1, 0, 0), "1:Del:C:0"),
    ("del A from A4", "AAAAGT", "del", 1, (4, 0, 0), "1:Del:T:3"),
    ("ins G, no G after", "TTGA", "ins", "G", (2, 0, 0), "1:Ins:C:0"),
    ("ins T into T7", "TTTTTTTGA", "ins", "T", (7, 0, 0), "1:Ins:T:5"),
    # >=2bp deletion: no MH -> Del:R:0 (fix 1), real MH -> Del:M
    ("del AGT, no MH", "AGTTTTTTT", "del", 3, (1, 0, 0), "3:Del:R:0"),
    ("del AT, MH1", "ATAGGG", "del", 2, (1, 0, 0), "2:Del:M:1"),
    ("del 12bp, MH5+", "ACGTTGCAAGCTACGTTGCGG", "del", 12, (1, 0, 0), "5:Del:M:5"),
    # STR, unit == event length
    ("del AC from (AC)5", "ACACACACACGG", "del", 2, (1, 2, 5), "2:Del:R:4"),
    ("ins AC into (AC)5", "ACACACACACGG", "ins", "AC", (1, 2, 5), "2:Ins:R:5"),
    ("ins AG before (AG)2", "AGAGCC", "ins", "AG", (1, 2, 2), "2:Ins:R:2"),
    # unannotated single copy (fix 2)
    ("ins TA before 1 copy", "TAGCC", "ins", "TA", (1, 0, 0), "2:Ins:R:1"),
    ("ins CAG before 1 copy", "CAGTTT", "ins", "CAG", (1, 0, 0), "3:Ins:R:1"),
    # event = k copies of the tract unit (fix 3)
    ("del ACAC from (AC)8", "AC" * 8 + "GG", "del", 4, (1, 2, 8), "4:Del:R:3"),
    ("ins ACAC into (AC)8", "AC" * 8 + "GG", "ins", "ACAC", (1, 2, 8), "4:Ins:R:4"),
    ("del CAGCAG from (CAG)6", "CAG" * 6 + "TT", "del", 6, (1, 3, 6), "5:Del:R:2"),
    (
        "ins CAGCAG into (CAG)6",
        "CAG" * 6 + "TT",
        "ins",
        "CAGCAG",
        (1, 3, 6),
        "5:Ins:R:3",
    ),
    # multi-base unit of one base: counted from hp_run
    ("ins TT into T5", "TTTTTGC", "ins", "TT", (5, 0, 0), "2:Ins:R:2"),
    ("del AA from A6", "AAAAAAGC", "del", 2, (6, 0, 0), "2:Del:R:2"),
    ("ins AAA into A13", "A" * 13 + "GC", "ins", "AAA", (13, 0, 0), "3:Ins:R:4"),
    # tract of the right unit length but another motif must not be used
    ("ins GG at (AG)2", "GGAGAGTC", "ins", "GG", (2, 2, 2), "2:Ins:R:1"),
    ("ins AAAT at (AAAG)2", "AAATAAAGAAAGC", "ins", "AAAT", (3, 4, 2), "4:Ins:R:1"),
    # >10bp unit, one following copy: fallback still matches SPMG
    (
        "del 12bp, 1 copy after",
        "ACGTTGCAAGCT" * 2 + "GG",
        "del",
        12,
        (1, 0, 0),
        "5:Del:R:1",
    ),
]


def _classify(after_anchor, kind, payload, anno=None):
    last = after_anchor[payload - 1] if kind == "del" else payload[-1]
    anchor = next(b for b in "GCTA" if b != last and b != after_anchor[0])
    ref_before = "TTGCA" + anchor
    if kind == "del":
        indel_seq = after_anchor[:payload]
        ref_after = after_anchor[payload:]
        indel_len = -payload
    else:
        indel_seq = payload
        ref_after = after_anchor
        indel_len = len(payload)
    label = classify_indel_channel(indel_seq, ref_before, ref_after, indel_len)
    return INDEL83_TO_SIGPROFILER_LABELS[label]


@pytest.mark.parametrize(
    "after_anchor,kind,payload,anno,expected",
    [c[1:] for c in CASES],
    ids=[c[0] for c in CASES],
)
def test_matches_sigprofiler(after_anchor, kind, payload, anno, expected):
    assert _classify(after_anchor, kind, payload, anno) == expected


def test_long_unit_multiple_copies():
    # >10bp unit (never in str.h5, PERF -M 10) with 2 following copies.
    after = "ACGTTGCAAGCT" * 3 + "GG"
    assert _classify(after, "del", 12) == "5:Del:R:2"


@pytest.mark.parametrize(
    "after_anchor,payload,expected",
    [
        # short repeats whose str.h5 tract lost to an overlapping longer one
        ("AGAGATGGGG", "AG", "2:Ins:R:2"),  # (AG)2 vs (TAG)2
        ("TTATTTATTTATTTG", "TTAT", "4:Ins:R:3"),  # (TTAT)3 vs (ATTATTT)2
        ("CACATCACATCACA", "CA", "2:Ins:R:2"),  # (CA)2 vs (CACAT)2
    ],
)
def test_overlapped_short_repeats(after_anchor, payload, expected):
    assert _classify(after_anchor, "ins", payload) == expected


def test_left_copies_counted():
    # Not left-aligned (VCFs from DupCaller always are), but SPMG also
    # counts copies ending at the anchor: ins AC after ...ACAC^ACGG.
    label = classify_indel_channel("AC", "GGACAC", "ACGG", 2)
    assert INDEL83_TO_SIGPROFILER_LABELS[label] == "2:Ins:R:3"
