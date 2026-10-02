"""Reference indel opportunity never credits N positions, matching the
coverage side (noise_mask drops N) and reference_base_number."""

import numpy as np

from DupCaller_sub.funcs.misc import indel100_reference_bucket_indices


def _opportunity(ref):
    n = len(ref)
    ref_int = np.array(["ATCG".index(b) if b in "ATCG" else 4 for b in ref])
    nxt = np.empty(n, dtype=np.int16)
    nxt[:-1] = ref_int[1:]
    nxt[-1] = -1
    zeros = np.zeros(n, dtype=int)
    return indel100_reference_bucket_indices(
        np.ones(n, dtype=int), zeros, zeros, ref_int, nxt, np.ones(n), zeros
    )


def test_n_positions_add_no_opportunity():
    # ACGT's last base has no valid next base either way (contig end vs N).
    assert np.array_equal(_opportunity("ACGT"), _opportunity("NNNACGTNNNN"))


def test_flat_columns_count_non_n_positions():
    out = _opportunity("ACGTNNNNAC")
    # Flat STR deletion rep1 (col 44) and microhomology (col 92).
    assert out[44] == 6 and out[92] == 6
