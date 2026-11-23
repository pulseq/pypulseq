"""Tests for the make_digital_output_pulse module"""

import numpy as np
import pytest
from pypulseq import Opts, make_digital_output_pulse


def test_channel_error():
    with pytest.raises(ValueError, match='is invalid\\. Must be one of'):
        make_digital_output_pulse(channel='physio1')


def test_negative_duration():
    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_digital_output_pulse(channel='osc0', duration=-1e-3)

    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_digital_output_pulse(channel='osc0', duration=np.nan)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_digital_output_pulse(channel='osc0', delay=-1e-6)


def test_valid():
    system = Opts()

    trig = make_digital_output_pulse(channel='osc0', duration=4e-3, system=system)
    assert trig.type == 'output'
    assert trig.duration == 4e-3

    # Durations shorter than the gradient raster are rounded up
    trig = make_digital_output_pulse(channel='ext1', duration=0, system=system)
    assert trig.duration == system.grad_raster_time
