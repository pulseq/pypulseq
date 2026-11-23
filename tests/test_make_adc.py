"""Tests for the make_adc module"""

import numpy as np
import pytest
from pypulseq import Opts, make_adc


def test_num_samples():
    for invalid in [0, -1, 1.5, np.nan, np.inf]:
        with pytest.raises(ValueError, match='num_samples must be a positive integer'):
            make_adc(num_samples=invalid, duration=1e-3)

    # Integral numpy ints work and are stored as plain int
    adc = make_adc(num_samples=np.int64(64), duration=1e-3)
    assert adc.num_samples == 64
    assert isinstance(adc.num_samples, int)


def test_dwell_duration_exclusivity():
    with pytest.raises(ValueError, match='Either dwell or duration must be defined'):
        make_adc(num_samples=64)

    with pytest.raises(ValueError, match='Either dwell or duration must be defined'):
        make_adc(num_samples=64, dwell=1e-5, duration=1e-3)


def test_negative_dwell():
    with pytest.raises(ValueError, match='Dwell must be non-negative'):
        make_adc(num_samples=64, dwell=-1e-6)


def test_negative_duration():
    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_adc(num_samples=64, duration=-1e-3)

    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_adc(num_samples=64, duration=np.nan)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_adc(num_samples=64, duration=1e-3, delay=-1e-6)


def test_phase_modulation_length():
    with pytest.raises(ValueError, match='same length as the number of samples'):
        make_adc(num_samples=64, duration=1e-3, phase_modulation=np.zeros(32))


def test_valid_generation():
    adc = make_adc(num_samples=64, duration=1e-3)
    assert adc.dwell == 1e-3 / 64
    assert adc.duration == 1e-3

    adc = make_adc(num_samples=64, dwell=1e-5)
    assert adc.duration == 64 * 1e-5


def test_dead_time_warning():
    system = Opts(adc_dead_time=100e-6)

    with pytest.warns(UserWarning, match='less than the dead time'):
        adc = make_adc(num_samples=64, duration=1e-3, delay=0, system=system)

    assert adc.delay == system.adc_dead_time
