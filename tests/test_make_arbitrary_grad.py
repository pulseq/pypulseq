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
