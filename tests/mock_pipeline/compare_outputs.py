#!/usr/bin/env python3
"""Compares a fresh mock-pipeline run against the premade expected/ outputs.

Text/VCF/tsv files are compared token-by-token; tokens that both parse as
floats are compared with a small tolerance (RTOL/ATOL below), everything
else must match exactly. .npz files are compared array-by-array with the
same tolerance.
"""
import sys
import numpy as np

RTOL = 1e-6
ATOL = 1e-9


def _is_float(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def compare_text(expected_path, actual_path):
    """Returns a list of human-readable diff messages (empty if equal)."""
    with open(expected_path) as f:
        expected_lines = f.readlines()
    with open(actual_path) as f:
        actual_lines = f.readlines()

    diffs = []
    if len(expected_lines) != len(actual_lines):
        diffs.append(
            f"line count differs: expected {len(expected_lines)}, "
            f"got {len(actual_lines)}"
        )
        return diffs

    for lineno, (exp_line, act_line) in enumerate(
        zip(expected_lines, actual_lines), start=1
    ):
        exp_toks = exp_line.rstrip("\n").split("\t")
        act_toks = act_line.rstrip("\n").split("\t")
        if len(exp_toks) != len(act_toks):
            diffs.append(
                f"line {lineno}: field count differs: {exp_line!r} vs {act_line!r}"
            )
            continue
        for col, (exp_tok, act_tok) in enumerate(zip(exp_toks, act_toks), start=1):
            if exp_tok == act_tok:
                continue
            if _is_float(exp_tok) and _is_float(act_tok):
                exp_val, act_val = float(exp_tok), float(act_tok)
                if np.isnan(exp_val) and np.isnan(act_val):
                    continue
                if abs(exp_val - act_val) <= ATOL + RTOL * abs(exp_val):
                    continue
            diffs.append(
                f"line {lineno} col {col}: expected {exp_tok!r}, got {act_tok!r}"
            )
    return diffs


def compare_npz(expected_path, actual_path):
    diffs = []
    expected = np.load(expected_path)
    actual = np.load(actual_path)
    if sorted(expected.files) != sorted(actual.files):
        diffs.append(
            f"array names differ: expected {expected.files}, got {actual.files}"
        )
        return diffs
    for key in expected.files:
        exp_arr, act_arr = expected[key], actual[key]
        if exp_arr.shape != act_arr.shape:
            diffs.append(
                f"array {key!r}: shape differs: {exp_arr.shape} vs {act_arr.shape}"
            )
            continue
        if np.issubdtype(exp_arr.dtype, np.floating):
            if not np.allclose(exp_arr, act_arr, rtol=RTOL, atol=ATOL, equal_nan=True):
                diffs.append(f"array {key!r}: values differ beyond tolerance")
        elif not np.array_equal(exp_arr, act_arr):
            diffs.append(f"array {key!r}: values differ")
    return diffs


def compare_file(expected_path, actual_path):
    if str(expected_path).endswith(".npz"):
        return compare_npz(expected_path, actual_path)
    return compare_text(expected_path, actual_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: compare_outputs.py EXPECTED ACTUAL")
    diffs = compare_file(sys.argv[1], sys.argv[2])
    for d in diffs:
        print(d)
    sys.exit(1 if diffs else 0)
