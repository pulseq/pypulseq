"""Tests for the make_extended_trapezoid module"""

import numpy as np
import pytest
from pypulseq import Opts, make_extended_trapezoid

system = Opts()
grt = system.grad_raster_time


def test_channel_error():
    with pytest.raises(ValueError, match="Invalid channel\\. Must be one of 'x', 'y' or 'z'"):
        make_extended_trapezoid(channel='p', amplitudes=[0, 1], times=[0, grt])


def test_times_amplitudes_length_mismatch():
    with pytest.raises(ValueError, match='Times and amplitudes must have the same length'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 1], times=[0, grt, 2 * grt])


def test_too_few_time_points():
    with pytest.raises(ValueError, match='At least two time points are required'):
        make_extended_trapezoid(channel='x', amplitudes=[1], times=[grt])


def test_all_zero_times():
    with pytest.raises(ValueError, match='At least one of the given times must be non-zero'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 0], times=[0, 0])


def test_times_not_ascending():
    with pytest.raises(ValueError, match='Times must be in ascending order'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 1, 0], times=[0, 2 * grt, grt])


def test_negative_first_time():
    with pytest.raises(ValueError, match='The first time point must be non-negative'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 0], times=[-grt, grt])


def test_nonfinite_times():
    with pytest.raises(ValueError, match='All elements of `times` must be finite'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 0], times=[0, np.nan])


def test_nonfinite_amplitudes():
    with pytest.raises(ValueError, match='All elements of `amplitudes` must be finite'):
        make_extended_trapezoid(channel='x', amplitudes=[0, np.inf], times=[0, grt])


def test_last_time_off_raster():
    with pytest.raises(ValueError, match='The last time point must be on a gradient raster'):
        make_extended_trapezoid(channel='x', amplitudes=[0, 0], times=[0, 1.5 * grt])


def test_first_amplitude_must_connect():
    with pytest.raises(ValueError, match='must connect to previous block'):
        make_extended_trapezoid(channel='x', amplitudes=[1, 0], times=[grt, 2 * grt], skip_check=False)

    # skip_check=True allows a non-zero first amplitude with times[0] > 0
    g = make_extended_trapezoid(channel='x', amplitudes=[1, 0], times=[grt, 2 * grt], skip_check=True)
    assert g.type == 'grad'


def test_valid():
    g = make_extended_trapezoid(
        channel='x', amplitudes=[0, 1000, 1000, 0], times=[0, grt, 3 * grt, 4 * grt], system=system
    )

    assert g.type == 'grad'
    assert g.area == pytest.approx(1000 * 3 * grt)
    assert g.first == 0
    assert g.last == 0


def test_convert_to_arbitrary():
    g = make_extended_trapezoid(
        channel='x',
        amplitudes=[0, 1000, 0],
        times=[0, 0.5 * grt, grt],
        convert_to_arbitrary=True,
        system=system,
    )

    assert g.type == 'grad'
    assert g.shape_dur == pytest.approx(grt)
