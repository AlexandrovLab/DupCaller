"""call --indel-context-trim: _indel_context_hi returns the exclusive end of an indel's scored repeat context (the
anchor through the base after the HP run / STR tract), and agrees with indels.indel_context_fits_window on whether
that context ends inside the window."""
import numpy as np
import pytest

from DupCaller_sub.funcs.call import _indel_context_hi
from DupCaller_sub.funcs.indels import indel_context_fits_window

_B2N = {"A": 0, "T": 1, "C": 2, "G": 3}


def _ref_int(seq):
    return np.array([_B2N[b] for b in seq])


def _hp_raw(seq):
    ref = _ref_int(seq)
    cut = np.ones(len(ref), dtype=bool)
    cut[1:] = ref[1:] != ref[:-1]
    run_id = np.cumsum(cut) - 1
    return np.vstack((np.bincount(run_id)[run_id], cut)).astype(float)


REF = "GCAAAATG"  # A run at local 2-5, base after it at 6
NO_STR = np.zeros([3, len(REF)])


@pytest.mark.parametrize(
    "indel, anchor, expected_hi",
    [
        ("10:-1", 1, 7),  # 1 bp deletion in the A run: context = run + the base after it
        ("10:1:A", 1, 7),  # 1 bp insertion of the run's base: same context
        ("10:1:G", 1, 3),  # 1 bp insertion of another base: only the first base after the anchor
        ("10:-2", 1, 3),  # 2 bp deletion with no STR tract: only the first deleted base
        ("10:-1", len(REF) - 1, len(REF) + 1),  # context starts past the window
    ],
)
def test_indel_context_hi(indel, anchor, expected_hi):
    assert _indel_context_hi(indel, anchor, _ref_int(REF), _hp_raw(REF), NO_STR) == expected_hi


@pytest.mark.parametrize("window_end", range(3, len(REF) + 1))
@pytest.mark.parametrize("indel", ["10:-1", "10:1:A", "10:1:G", "10:-2"])
def test_agrees_with_indel_context_fits_window(indel, window_end):
    ref = REF[:window_end]
    ref_int, hp_raw, str_raw = _ref_int(ref), _hp_raw(ref), NO_STR[:, :window_end]
    hi = _indel_context_hi(indel, 1, ref_int, hp_raw, str_raw)
    assert (hi <= window_end) == indel_context_fits_window(indel, 1, ref_int, hp_raw, str_raw)
