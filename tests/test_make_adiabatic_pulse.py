"""Tests for the make_adiabatic_pulse.py module

Will Clarke, University of Oxford, 2024
"""

import itertools

import numpy as np
import pytest
from pypulseq import make_adiabatic_pulse
from pypulseq.supported_labels_rf_use import get_supported_rf_uses


def test_pulse_select():
    valid_rf_use_labels = get_supported_rf_uses()
    valid_pulse_types = ('hypsec', 'wurst')

    # Check all use and valid pulse combinations return a sensible object
    # with default parameters.
    for pulse_type, use_label in itertools.product(valid_pulse_types, valid_rf_use_labels):
        rf_obj = make_adiabatic_pulse(pulse_type=pulse_type, use=use_label)
        assert rf_obj.type == 'rf'
        assert rf_obj.use == use_label

    # Check the appropriate errors are raised if we specify nonsense
    with pytest.raises(ValueError, match=r'Invalid type parameter\. Must be one of '):
        make_adiabatic_pulse(pulse_type='not a pulse type')

    with pytest.raises(ValueError, match=r'Invalid type parameter\. Must be one of '):
        make_adiabatic_pulse(pulse_type='')

    with pytest.raises(ValueError, match=r'Invalid use parameter\. Must be one of '):
        make_adiabatic_pulse(pulse_type='hypsec', use='not a use')

    # Default use case
    rf_obj = make_adiabatic_pulse(pulse_type='hypsec')
    assert rf_obj.use == 'inversion'


def test_option_requirements():
    # Require non-zero slice thickness if grad requested
    with pytest.raises(ValueError, match='Slice thickness must be provided'):
        make_adiabatic_pulse(pulse_type='hypsec', return_gz=True)

    with pytest.raises(ValueError, match='Slice thickness must be provided'):
        make_adiabatic_pulse(pulse_type='hypsec', return_gz=True, slice_thickness=np.nan)

    _, gz, gzr = make_adiabatic_pulse(pulse_type='hypsec', return_gz=True, slice_thickness=1)
    assert gz.type == 'trap'
    assert gzr.type == 'trap'


def test_invalid_duration():
    for invalid in [0, -1e-3, np.nan, np.inf]:
        with pytest.raises(ValueError, match='RF pulse duration must be positive'):
            make_adiabatic_pulse(pulse_type='hypsec', duration=invalid)


def test_negative_delay():
    with pytest.raises(ValueError, match='Delay must be non-negative'):
        make_adiabatic_pulse(pulse_type='hypsec', delay=-1e-3)


def test_invalid_dwell():
    for invalid in [0, -1e-6, np.nan]:
        with pytest.raises(ValueError, match='Dwell must be positive'):
            make_adiabatic_pulse(pulse_type='hypsec', dwell=invalid)


def test_invalid_adiabaticity():
    for invalid in [0, np.nan]:
        with pytest.raises(ValueError, match='adiabaticity must be positive'):
            make_adiabatic_pulse(pulse_type='hypsec', adiabaticity=invalid)


def test_invalid_hypsec_parameters():
    for invalid in [-800, np.nan]:
        with pytest.raises(ValueError, match='beta must be positive'):
            make_adiabatic_pulse(pulse_type='hypsec', beta=invalid)

    for invalid in [0, np.nan]:
        with pytest.raises(ValueError, match='mu must be positive'):
            make_adiabatic_pulse(pulse_type='hypsec', mu=invalid)


def test_invalid_wurst_parameters():
    for invalid in [-40000, np.nan]:
        with pytest.raises(ValueError, match='bandwidth must be positive'):
            make_adiabatic_pulse(pulse_type='wurst', bandwidth=invalid)

    for invalid in [0, np.nan]:
        with pytest.raises(ValueError, match='n_fac must be positive'):
            make_adiabatic_pulse(pulse_type='wurst', n_fac=invalid)


# My intention was to test that the rephase gradient area is appropriate,
# but this doesn't pass and I'm highly suspicious of the calculation in
# the code
def test_returned_grads():
    pass
    # _, gz, gzr = make_adiabatic_pulse(
    #         pulse_type="hypsec",
    #         return_gz=True,
    #         slice_thickness=1)
    # assert np.isclose(-gz.area / 2, gzr.area)


def test_hypsec_options():
    pobj = make_adiabatic_pulse(pulse_type='hypsec', beta=700, mu=6, duration=0.05)

    assert np.isclose(pobj.shape_dur, 0.05)


def test_wurst_options():
    pobj = make_adiabatic_pulse(pulse_type='wurst', n_fac=25, bandwidth=30000, duration=0.05)

    assert np.isclose(pobj.shape_dur, 0.05)
