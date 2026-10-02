import itertools

import numpy as np
from scipy.stats import binom

from DupCaller_sub.funcs.prob import strand_nonalt_pvalues


def test_equal_bq_matches_binomial():
    # 12 reads at BQ 37, 5 non-alt (a 7 alt / 5 ref strand); 1 position.
    qual = np.full((12, 1), 37.0)
    alt = np.array([[True]] * 7 + [[False]] * 5)
    p = strand_nonalt_pvalues(qual, alt, np.array([0.0]))
    assert np.isclose(p[0], binom.sf(4, 12, 10**-3.7), rtol=1e-6)


def test_mixed_bq_matches_bruteforce():
    quals = np.array([37.0, 25.0, 37.0, 30.0, 11.0])
    alt = np.array([True, False, True, False, True])
    pamp = 1e-4
    d = 10 ** (-quals / 10) + pamp
    k = 2
    brute = sum(
        np.prod([d[i] if e else 1 - d[i] for i, e in enumerate(err)])
        for err in itertools.product([0, 1], repeat=5)
        if sum(err) >= k
    )
    p = strand_nonalt_pvalues(quals[:, None], alt[:, None], np.array([pamp]))
    assert np.isclose(p[0], brute, rtol=1e-9)


def test_uncounted_reads_ignored_and_clean_is_one():
    # Column 0: clean strand (all alt) -> p == 1.
    # Column 1: one counted non-alt read plus a zero-quality (uncounted) one.
    qual = np.array([[37.0, 37.0], [37.0, 0.0], [37.0, 37.0]])
    alt = np.array([[True, True], [True, False], [True, False]])
    p = strand_nonalt_pvalues(qual, alt, np.array([0.0, 0.0]))
    assert p[0] == 1.0
    assert np.isclose(p[1], 1 - (1 - 10**-3.7) ** 2, rtol=1e-6)


def test_no_positions():
    assert (
        strand_nonalt_pvalues(
            np.zeros((3, 0)), np.zeros((3, 0), bool), np.zeros(0)
        ).size
        == 0
    )


def test_tiny_tail_rounds_to_zero_not_negative():
    # 5 of 5 reads non-alt at BQ 37: true p ~ 3e-19, which 1 - cdf rounds
    # to 0; it must stay a valid probability.
    qual = np.full((5, 1), 37.0)
    alt = np.zeros((5, 1), bool)
    p = strand_nonalt_pvalues(qual, alt, np.array([0.0]))
    assert 0.0 <= p[0] < 1e-15
