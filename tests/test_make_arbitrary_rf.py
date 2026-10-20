"""Tests for the make_arbitrary_rf module"""

import numpy as np
import pytest
from pypulseq import Opts, make_arbitrary_rf


def test_invalid_use():
    with pytest.raises(ValueError, match='Invalid use parameter'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, use='invalid')


def test_signal_ndim():
    with pytest.raises(ValueError, match='signal should have ndim=1'):
        make_arbitrary_rf(signal=np.zeros((2, 5)), flip_angle=1)

    with pytest.raises(ValueError, match='signal should have ndim=1'):
        make_arbitrary_rf(signal=np.zeros((2, 2, 5)), flip_angle=1)


def test_empty_signal():
    with pytest.raises(ValueError, match='signal must not be empty'):
        make_arbitrary_rf(signal=np.zeros(0), flip_angle=1)


def test_negative_dwell():
    with pytest.raises(ValueError, match='Dwell must be non-negative'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, dwell=-1e-6)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, delay=-1e-6)


def test_negative_time_bw_product():
    with pytest.raises(ValueError, match='time_bw_product must be non-negative'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, time_bw_product=-1)


def test_return_gz_requires_slice_thickness_and_bandwidth():
    with pytest.raises(ValueError, match='Slice thickness must be provided'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, return_gz=True, slice_thickness=0, bandwidth=4e3)

    with pytest.raises(ValueError, match='Slice thickness must be provided'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, return_gz=True, slice_thickness=-1e-3, bandwidth=4e3)

    with pytest.raises(ValueError, match='Bandwidth of pulse must be provided'):
        make_arbitrary_rf(signal=np.ones(10), flip_angle=1, return_gz=True, slice_thickness=5e-3, bandwidth=0)


def test_valid():
    system = Opts()
    rf = make_arbitrary_rf(signal=np.ones(10), flip_angle=1, system=system)

    assert rf.type == 'rf'
    assert rf.shape_dur == pytest.approx(10 * system.rf_raster_time)
    assert rf.use == 'undefined'

    # 1-D signal is also accepted as a row/column vector
    rf = make_arbitrary_rf(signal=np.ones((1, 10)), flip_angle=1)
    assert rf.shape_dur == pytest.approx(10 * system.rf_raster_time)
