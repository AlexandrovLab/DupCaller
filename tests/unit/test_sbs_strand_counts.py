"""_sbs_strand_counts (one bincount over all of a strand's reads) gives the
same SBS amp tallies as the per-read loop it replaced."""

import numpy as np
import pytest

from DupCaller_sub.funcs.learn import MAX_BQ, NUM_BQ, _sbs_strand_counts


def _per_read(seq_mat, qual_mat, antimask, trinuc_int):
    trinuc_masked = trinuc_int[antimask]
    count_mat = np.zeros([96, 4])
    bq_hist = np.zeros([96, 4, NUM_BQ])
    for mm in range(seq_mat.shape[0]):
        seq_masked = seq_mat[mm, antimask]
        qual_masked = qual_mat[mm, antimask]
        alt_1d = trinuc_masked + seq_masked * 96
        valid = (qual_masked > 0) & (seq_masked < 4)
        count_mat += (
            np.bincount(alt_1d, weights=valid.astype(float), minlength=96 * 4)[
                0 : 4 * 96
            ]
            .reshape([4, 96])
            .T
        )
        if valid.any():
            bq_valid = np.clip(qual_masked[valid].astype(int), 0, MAX_BQ)
            flat_idx = alt_1d[valid] * NUM_BQ + bq_valid
            bq_hist += (
                np.bincount(flat_idx, minlength=4 * 96 * NUM_BQ)
                .reshape([4, 96, NUM_BQ])
                .transpose(1, 0, 2)
            )
    return count_mat, bq_hist


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("n_reads", [0, 1, 7])
def test_matches_per_read_loop(seed, n_reads):
    rng = np.random.default_rng(seed)
    n = 150
    seq_mat = rng.integers(0, 5, size=(n_reads, n))  # 4 = N/deletion/off-read
    qual_mat = rng.choice([0, 0, 2, 11, 25, 37, 40], size=(n_reads, n)).astype(float)
    antimask = rng.random(n) < 0.7
    trinuc_int = rng.integers(0, 65, size=n)
    count_mat, bq_hist = _sbs_strand_counts(seq_mat, qual_mat, antimask, trinuc_int)
    ref_count, ref_hist = _per_read(seq_mat, qual_mat, antimask, trinuc_int)
    assert count_mat.dtype == ref_count.dtype and bq_hist.dtype == ref_hist.dtype
    assert np.array_equal(count_mat, ref_count)
    assert np.array_equal(bq_hist, ref_hist)
