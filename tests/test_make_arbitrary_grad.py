"""Tests for the make_arbitrary_grad module"""

import re

import numpy as np
import pytest
from pypulseq import Opts, make_arbitrary_grad

system = Opts(max_slew=100, slew_unit='T/m/s')


def oversampled_triangle(slew_fraction):
    """Oversampled triangle of 17 samples, starting and ending at 0 (`first=0`, `last=0`).

    The samples are `grad_raster_time / 2` apart, so every segment has a slope of
    `slew_fraction * max_slew`.
    """
    step = slew_fraction * system.max_slew * system.grad_raster_time / 2
    return step * np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1])


def test_oversampled_slew_within_limit():
    g = make_arbitrary_grad('x', oversampled_triangle(0.9), first=0, last=0, oversampling=True, system=system)

    # The slope of each segment of the event, from its own tt, waveform, first and last.
    t = np.concatenate([[0], g.tt, [g.shape_dur]])
    amp = np.concatenate([[g.first], g.waveform, [g.last]])
    assert np.allclose(np.abs(np.diff(amp) / np.diff(t)), 0.9 * system.max_slew)


@pytest.mark.parametrize('slew_fraction', [1.5, 3.5])
def test_oversampled_slew_violation(slew_fraction):
    with pytest.raises(ValueError, match=r'Slew rate violation') as err:
        make_arbitrary_grad('x', oversampled_triangle(slew_fraction), first=0, last=0, oversampling=True, system=system)

    # The error reports the slew rate in percent of max_slew.
    reported = float(re.search(r'Slew rate violation ([\d.]+)', str(err.value)).group(1))
    assert reported == pytest.approx(slew_fraction * 100)


def test_channel_error():
    with pytest.raises(ValueError, match='Invalid channel'):
        make_arbitrary_grad('p', np.zeros(10))


def test_waveform_ndim():
    with pytest.raises(ValueError, match='waveform should have ndim=1'):
        make_arbitrary_grad('x', np.zeros((2, 5)))


def test_waveform_empty():
    with pytest.raises(ValueError, match='waveform must not be empty'):
        make_arbitrary_grad('x', np.zeros(0), first=0, last=0)


def test_waveform_too_short():
    # A single-sample waveform requires explicit first and last values
    with pytest.raises(ValueError, match='at least 2 samples'):
        make_arbitrary_grad('x', np.zeros(1))

    g = make_arbitrary_grad('x', np.zeros(1), first=0, last=0)
    assert g.type == 'grad'


def test_oversampling_even_samples():
    with pytest.raises(ValueError, match='odd number of samples'):
        make_arbitrary_grad('x', np.zeros(10), oversampling=True)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_arbitrary_grad('x', np.zeros(10), delay=-1e-6)
