"""Per-locus coverage carried between 1 Mb windows is re-added at its own
genomic position, even when the next window starts upstream of the key
that triggered it (a rerouted rugged mate pulls batch_min_start back)."""

import numpy as np

from DupCaller_sub.funcs.call import (
    _carry_in_coverage,
    _coverage_window_rows,
    _cut_coverage_leftover,
)


def _window(start, rows, marks):
    arr = np.zeros((rows, 4))
    for g in marks:
        arr[g - start, 0] = g
    return arr


def _positions(arr, start):
    return {start + int(i) for i in np.flatnonzero(arr[:, 0])}


def test_carry_over_keeps_genomic_position():
    old_start, new_start = 1000, 1090
    old = _window(old_start, 200, [1050, 1089, 1090, 1095, 1150])
    flush, left, gstart = _cut_coverage_leftover([old], old_start, new_start, True)
    assert _positions(old[:flush], old_start) == {1050, 1089}
    rows = _coverage_window_rows(left, gstart, new_start)
    new = np.zeros((rows, 4))
    _carry_in_coverage([new], left, gstart, new_start)
    assert _positions(new, new_start) == {1090, 1095, 1150}
    assert all(new[g - new_start, 0] == g for g in (1090, 1095, 1150))


def test_upstream_new_start_carries_whole_window():
    old_start, new_start = 1000, 990
    old = _window(old_start, 100, [1000, 1050])
    flush, left, gstart = _cut_coverage_leftover([old], old_start, new_start, True)
    assert flush == 0 and gstart == old_start
    new = np.zeros((_coverage_window_rows(left, gstart, new_start), 4))
    _carry_in_coverage([new], left, gstart, new_start)
    assert _positions(new, new_start) == {1000, 1050}


def test_chromosome_change_flushes_everything():
    old = _window(1000, 100, [1000, 1099])
    flush, left, gstart = _cut_coverage_leftover([old], 1000, 5, False)
    assert flush == 100 and left is None and gstart is None
