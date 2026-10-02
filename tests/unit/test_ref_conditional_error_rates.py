import numpy as np
import pytest

from DupCaller_sub.funcs.misc import (
    _hp_ref_conditional,
    _sbs_ref_cols,
    _str_ref_conditional,
    ref_conditional_error_rates,
)


def test_sbs_ref_is_one_minus_errors_and_errors_conditional():
    rates = np.full([64, 4], 1e-3)
    ref_cols = _sbs_ref_cols()
    rates[np.arange(64), ref_cols] = 0.5  # deliberately inconsistent ref
    rates[0, (ref_cols[0] + 1) % 4] = 2e-3
    out = ref_conditional_error_rates(rates, ref_cols, np.zeros([64, 4]))
    ref = out[np.arange(64), ref_cols]
    assert np.allclose(ref[1:], 1 - 3e-3)
    assert np.isclose(ref[0], 1 - 4e-3)
    alt = (ref_cols[0] + 1) % 4
    assert np.isclose(out[0, alt], 2e-3 / (2e-3 + 1 - 4e-3))
    other = (ref_cols[0] + 2) % 4
    assert np.isclose(out[0, other], 1e-3 / (1e-3 + 1 - 4e-3))


def test_negative_ref_row_takes_fallback_row():
    rates = np.array([[0.1, 0.6, 0.6, 0.0], [0.97, 0.01, 0.01, 0.01]])
    fb = np.array([[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25]])
    out = ref_conditional_error_rates(rates, [0, 0], fb)
    assert np.isclose(out[0, 0], 0.7)
    assert np.allclose(out[0, 1:], 0.1 / 0.8)
    assert np.isclose(out[1, 0], 0.97)
    assert np.allclose(out[1, 1:], 0.01 / 0.98)


def test_bad_fallback_raises():
    rates = np.array([[0.0, 0.6, 0.6]])
    with pytest.raises(ValueError):
        ref_conditional_error_rates(rates, 0, np.array([[0.0, 0.6, 0.6]]))


def test_hp_ref_is_one_minus_del_ins_per_base():
    rates = np.tile([0.01, 0.2, 0.02], (10, 4))  # ref col not 1-del-ins
    out = _hp_ref_conditional(rates, np.ones([10, 12]), 0.5)
    assert np.allclose(out[:, 1::3], 0.97)
    assert np.allclose(out[:, 0::3], 0.01 / 0.98)
    assert np.allclose(out[:, 2::3], 0.02 / 0.99)


def test_hp_negative_ref_uses_fallback_line():
    rates = np.tile([0.01, 0.98, 0.01], (10, 4))
    rates[4, 3:6] = [0.6, 0.0, 0.5]  # HP5 T group: del+ins > 1
    fb = np.tile([9.5, 979.5, 9.5], (10, 4))  # -> 0.01/0.98/0.01 smoothed
    out = _hp_ref_conditional(rates, fb, 0.5)
    assert np.isclose(out[4, 4], 0.98)
    assert np.isclose(out[4, 3], 0.01 / 0.99)
    assert np.array_equal(
        out[:, :3],
        _hp_ref_conditional(np.tile([0.01, 0.98, 0.01], (10, 4)), fb, 0.5)[:, :3],
    )


def test_str_ref_is_one_minus_row():
    rates = np.full([5, 11], 0.01)
    rates[:, 5] = 0.5
    rates[2, :] = 0.2  # sums past 1 -> fallback
    fb = np.full([5, 11], 9.5)
    fb[:, 5] = 899.5  # smoothed: 0.01 each, ref 0.9
    out = _str_ref_conditional(rates, fb, 0.5)
    idx = np.arange(11) != 5
    assert np.allclose(out[:, 5], 0.9)
    assert np.allclose(out[:, idx], 0.01 / 0.91)
