"""Tests for the make_trigger module"""

import numpy as np
import pytest
from pypulseq import Opts, make_trigger


def test_channel_error():
    with pytest.raises(ValueError, match='is invalid\\. Must be one of'):
        make_trigger(channel='osc0')


def test_negative_duration():
    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_trigger(channel='physio1', duration=-1e-3)

    with pytest.raises(ValueError, match='Duration must be non-negative'):
        make_trigger(channel='physio1', duration=np.nan)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_trigger(channel='physio1', delay=-1e-6)


def test_valid():
    system = Opts()

    trigger = make_trigger(channel='physio1', duration=4e-3, system=system)
    assert trigger.type == 'trigger'
    assert trigger.duration == 4e-3

    # Durations shorter than the gradient raster are rounded up
    trigger = make_trigger(channel='physio2', duration=0, system=system)
    assert trigger.duration == system.grad_raster_time
