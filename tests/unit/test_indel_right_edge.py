"""genotypeDSIndel near the right edge of the loaded window: context
lookups (reference_int/hp_raw/str_raw at anchor+1) must never index past
the window; a candidate without its context base is dropped, not scored."""

import numpy as np
import pysam
import pytest

from DupCaller_sub.funcs.misc import fallback_error_file, load_error_matrices
from DupCaller_sub.funcs.prob import genotypeDSIndel

HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
)
REF = "ACGTTGCAAGCTTACG"  # window 100..115
REF_INT = np.array(["ATCG".index(b) for b in REF])
N = len(REF)


@pytest.fixture(scope="module")
def params():
    pre = fallback_error_file(".amp.tn.srd.txt")[: -len(".amp.tn.srd.txt")]
    p = {
        "amperr_file": pre + ".amp.tn.srd.txt",
        "dmgerr_file": pre + ".dmg.tn.txt",
        "amperri_file": pre + ".amp.id.txt",
        "dmgerri_file": pre + ".dmg.id.txt",
        "pseudocount": 0.5,
        "minBq": 10,
    }
    load_error_matrices(p)
    return p


def _read(cigar, seq, read1):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = "r"
    rec.query_sequence = seq
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    rec.reference_id = 0
    rec.reference_start = 100
    rec.mapping_quality = 60
    rec.cigar = cigar
    rec.flag = 1 | (64 if read1 else 128)  # paired, forward
    return rec


def _run(params, cigar, seq):
    reads = [_read(cigar, seq, True) for _ in range(3)]
    reads += [_read(cigar, seq, False) for _ in range(3)]
    hp_raw = np.zeros([3, N])
    hp_raw[0] = 1
    str_raw = np.zeros([3, N])
    return genotypeDSIndel(
        reads, 100, 100 + N, REF_INT, np.ones(N, bool), hp_raw, str_raw, params
    )


def test_insertion_after_final_base_dropped(params):
    # TT doesn't left-align past the final G, so its context base is
    # beyond the window.
    out = _run(params, [(0, N), (1, 2)], REF + "TT")
    assert len(out[2]) == 0


def test_insertion_one_before_final_base_scored(params):
    out = _run(params, [(0, N - 1), (1, 2), (0, 1)], REF[:-1] + "TT" + REF[-1])
    assert out[2] == [f"{100 + N - 2}:2:TT"]
    assert out[6][0] == 3 and out[8][0] == 3  # anchored: ALT on both strands


def test_long_deletion_reaching_right_edge(params):
    # 5bp deletion ending one base before the window end: scored.
    out = _run(params, [(0, 10), (2, 5), (0, 1)], REF[:10] + REF[15])
    assert out[2] == ["109:-5"]
    # Deletion running through the last window base: in-window, no crash;
    # its only evidence is an unanchored trailing deletion, so no read
    # counts as either allele.
    out = _run(params, [(0, 11), (2, 5)], REF[:11])
    assert out[2] == ["110:-5"]
    assert out[5][0] == out[6][0] == out[7][0] == out[8][0] == 0
