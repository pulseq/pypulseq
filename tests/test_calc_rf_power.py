from types import SimpleNamespace

import numpy as np
import pypulseq as pp
import pytest


@pytest.fixture
def dummy_rf():
    """
    Fixture factory for creating SimpleNamespace RF pulses.
    Returns a function you can call with parameters.
    """

    def _make_rf(
        amplitude=1.0,
        duration=1.0,
        t_override=None,
    ):
        system = pp.Opts.default

        if t_override is not None:
            t = np.asarray(t_override)
            duration = t[-1] - t[0]
            signal = amplitude * np.ones_like(t)
        else:
            n = int(duration / system.rf_raster_time) + 1
            t = np.linspace(0, duration, n, endpoint=True)
            signal = amplitude * np.ones_like(t)

        rf = SimpleNamespace()
        rf.type = 'rf'
        rf.signal = signal
        rf.t = t
        rf.shape_dur = duration
        rf.center = 0.5 * duration
        rf.freq_offset = 0.0
        rf.phase_offset = 0.0
        rf.freq_ppm = 0.0
        rf.phase_ppm = 0.0
        rf.dead_time = system.rf_dead_time
        rf.ringdown_time = system.rf_ringdown_time
        rf.delay = 0.0
        rf.use = 'excitation'

        return rf

    return _make_rf


@pytest.fixture
def dummy_sequence():
    """
    Fixture for creating a minimal pypulseq.Sequence with a helper
    to add SimpleNamespace RF pulses as blocks.
    """
    seq = pp.Sequence()

    return seq


def test_constant_rf(dummy_rf):
    """Test constant amplitude pulse against analytic solution."""
    A = 2.0
    T = 4.0e-3
    rf = dummy_rf(amplitude=A, duration=T)

    E, P, rms = pp.calc_rf_power(rf)

    assert np.isclose(E, A**2 * T)
    assert np.isclose(P, A**2)
    assert np.isclose(rms, A)


def test_nonuniform_raster(dummy_rf):
    """Test nonuniformly sampled RF returns correct energy."""
    A = 2.0
    t = [0.0e-3, 2.0e-3]  # nonuniform time samples
    rf = dummy_rf(amplitude=A, t_override=t)

    E, P, rms = pp.calc_rf_power(rf)

    assert np.isclose(E, A**2 * (t[-1] - t[0]))
    assert np.isclose(P, A**2)
    assert np.isclose(rms, A)


def test_resampled_rf(dummy_rf):
    """Test resampling path with dt parameter."""
    A = 2.0
    t = [0.0e-3, 2.0e-3]  # nonuniform time samples
    rf = dummy_rf(amplitude=A, t_override=t)

    E, P, rms = pp.calc_rf_power(rf, dt=pp.Opts.default.rf_raster_time)

    assert np.isclose(E, A**2 * (t[-1] - t[0]), rtol=1e-3)
    assert np.isclose(P, A**2, rtol=1e-3)
    assert np.isclose(rms, A, rtol=1e-3)


def test_sequence_single_block(dummy_rf, dummy_sequence):
    """Test sequence with a single RF block."""
    seq = dummy_sequence
    A = 2.0
    T = 5.0
    rf = dummy_rf(amplitude=A, duration=T)
    seq.add_block(rf)

    mean_pwr, peak, rms, energy = seq.calc_rf_power()

    assert np.isclose(energy, A**2 * T)
    assert np.isclose(mean_pwr, A**2)
    assert np.isclose(rms, A)
    assert np.isclose(peak, A**2)


def test_sequence_partial_overlap(dummy_rf, dummy_sequence):
    """Test sequence block partially overlapping time_range."""
    seq = dummy_sequence
    A = 2.0
    T = 10.0e-3
    rf = dummy_rf(amplitude=A, duration=T)
    seq.add_block(rf)

    mean_pwr, _, rms, energy = seq.calc_rf_power(time_range=(2.0e-3, 6.0e-3))

    # Only 4s of 10s contribute
    assert np.isclose(energy, A**2 * 4.0e-3)
    assert np.isclose(mean_pwr, A**2)
    assert np.isclose(rms, A)


def test_sequence_sliding_window(dummy_rf, dummy_sequence):
    """
    Test sliding window RF power calculation.

    NOTE:
    The sliding window implementation intentionally follows Pulseq's MATLAB behavior, where the accumulated energy is
    compared against the maximum *before* trimming the window back to the nominal window_duration.

    This leads to a conservative (slightly overestimated) mean power, as explicitly documented in the original Pulseq
    code and consistent with SAR-style worst-case windowing.
    """
    seq = dummy_sequence

    # Two consecutive RF blocks:
    # Block 1: A=1, duration=5s → energy = 1^2 * 5 = 5
    # Block 2: A=3, duration=5s → energy = 3^2 * 5 = 45
    rf1 = dummy_rf(amplitude=1.0, duration=5.0)
    rf2 = dummy_rf(amplitude=3.0, duration=5.0)

    seq.add_block(rf1)
    seq.add_block(rf2)

    mean_pwr, peak, rms, _ = seq.calc_rf_power(window_duration=5.0)

    # During accumulation, both blocks temporarily fall inside the window:
    # total_energy = 5 + 45 = 50
    #
    # mean_pwr is computed using the nominal window duration (5s),
    # resulting in a conservative overestimate:
    # mean_pwr = 50 / 5 = 10
    assert np.isclose(mean_pwr, 10.0)

    # RMS amplitude follows the same conservative windowing logic
    assert np.isclose(rms, np.sqrt(10.0))

    # Peak power is taken from the strongest RF block
    assert np.isclose(peak, 9.0)


def test_native_integration_dt_none(dummy_rf):
    """Test the dt=None path, which integrates the native rf.t samples."""
    A = 2.0
    T = 4.0e-3
    rf = dummy_rf(amplitude=A, duration=T)

    E, P, rms = pp.calc_rf_power(rf, dt=None)

    assert np.isclose(E, A**2 * T)
    assert np.isclose(P, A**2)
    assert np.isclose(rms, A)


def test_block_pulse_energy():
    """Block pulse energy against analytic solution (cf. MATLAB testCalcRfPower)."""
    duration = 1e-3
    rf = pp.make_block_pulse(flip_angle=np.pi / 2, duration=duration)
    amp = np.max(np.abs(rf.signal))

    total_energy, peak_pwr, rf_rms = pp.calc_rf_power(rf)

    assert np.isclose(total_energy, amp**2 * duration, rtol=0.02)
    assert np.isclose(peak_pwr, amp**2, rtol=0.02)
    assert np.isclose(rf_rms, amp, rtol=0.02)


def test_rms_relation_sinc():
    """rf_rms must equal sqrt(total_energy / shape_dur) (cf. MATLAB testCalcRfPower)."""
    rf = pp.make_sinc_pulse(flip_angle=np.pi / 2, duration=4e-3, time_bw_product=4)

    total_energy, _, rf_rms = pp.calc_rf_power(rf)

    assert np.isclose(rf_rms, np.sqrt(total_energy / rf.shape_dur), rtol=1e-6)


def test_dt_consistency():
    """Different dt values give consistent energy (cf. MATLAB testCalcRfPower)."""
    rf = pp.make_block_pulse(flip_angle=np.pi / 2, duration=1e-3)

    e1, _, _ = pp.calc_rf_power(rf, dt=1e-6)
    e2, _, _ = pp.calc_rf_power(rf, dt=0.5e-6)
    e3, _, _ = pp.calc_rf_power(rf, dt=None)

    assert np.isclose(e1, e2, rtol=0.02)
    assert np.isclose(e1, e3, rtol=0.02)


def test_sequence_window_with_delay_blocks(dummy_rf, dummy_sequence):
    """Sliding window must slide correctly past blocks without RF events."""
    seq = dummy_sequence
    A = 2.0
    T = 1.0e-3
    rf = dummy_rf(amplitude=A, duration=T)

    seq.add_block(rf)
    seq.add_block(pp.make_delay(10e-3))

    mean_pwr, peak, rms, energy = seq.calc_rf_power(window_duration=5e-3)

    # The largest window contains only the RF block:
    # energy = A**2 * T, averaged over the nominal window duration
    assert np.isclose(energy, A**2 * T)
    assert np.isclose(mean_pwr, A**2 * T / 5e-3)
    assert np.isclose(rms, np.sqrt(A**2 * T / 5e-3))
    assert np.isclose(peak, A**2)


def test_sequence_window_multiple_rf_and_delays(dummy_rf, dummy_sequence):
    """Sliding window over RF blocks separated by a delay block."""
    seq = dummy_sequence
    rf1 = dummy_rf(amplitude=2.0, duration=1e-3)  # energy 4e-3
    rf2 = dummy_rf(amplitude=4.0, duration=1e-3)  # energy 16e-3

    seq.add_block(rf1)
    seq.add_block(pp.make_delay(2e-3))
    seq.add_block(rf2)

    # window too small to contain both RF blocks: only the stronger counts
    mean_pwr, _, rms, energy = seq.calc_rf_power(window_duration=2.5e-3)
    assert np.isclose(energy, 16e-3)
    assert np.isclose(mean_pwr, 16e-3 / 2.5e-3)
    assert np.isclose(rms, np.sqrt(16e-3 / 2.5e-3))

    # window large enough for both RF blocks and the delay between them
    mean_pwr, _, rms, energy = seq.calc_rf_power(window_duration=4e-3)
    assert np.isclose(energy, 20e-3)
    assert np.isclose(mean_pwr, 20e-3 / 4e-3)
    assert np.isclose(rms, np.sqrt(20e-3 / 4e-3))


def test_sequence_block_range(dummy_rf, dummy_sequence):
    """block_range restricts the calculation to 1-based block indices, like MATLAB blockRange."""
    seq = dummy_sequence
    rf1 = dummy_rf(amplitude=1.0, duration=1e-3)  # energy 1e-3
    rf2 = dummy_rf(amplitude=2.0, duration=1e-3)  # energy 4e-3
    rf3 = dummy_rf(amplitude=4.0, duration=1e-3)  # energy 16e-3

    seq.add_block(rf1)
    seq.add_block(rf2)
    seq.add_block(rf3)

    mean_pwr, peak, rms, energy = seq.calc_rf_power(block_range=(2, 3))
    assert np.isclose(energy, 20e-3)
    assert np.isclose(mean_pwr, 20e-3 / 2e-3)
    assert np.isclose(rms, np.sqrt(20e-3 / 2e-3))
    assert np.isclose(peak, 16.0)

    # single block
    mean_pwr, peak, _, energy = seq.calc_rf_power(block_range=(2, 2))
    assert np.isclose(energy, 4e-3)
    assert np.isclose(mean_pwr, 4.0)
    assert np.isclose(peak, 4.0)

    # inf end selects all remaining blocks
    _, _, _, energy = seq.calc_rf_power(block_range=(2, np.inf))
    assert np.isclose(energy, 20e-3)

    # default covers all blocks
    _, _, _, energy = seq.calc_rf_power()
    assert np.isclose(energy, 21e-3)


def test_sequence_range_mutually_exclusive(dummy_rf, dummy_sequence):
    seq = dummy_sequence
    seq.add_block(dummy_rf(amplitude=1.0, duration=1e-3))

    with pytest.raises(ValueError, match='either time_range or block_range'):
        seq.calc_rf_power(time_range=(0, 1e-3), block_range=(1, 1))


def test_invalid_time_range(dummy_rf, dummy_sequence):
    seq = dummy_sequence
    seq.add_block(dummy_rf(amplitude=1.0, duration=1e-3))

    for bad_range in [5.0, (0.0,), (0.0, 1e-3, 2e-3), (2e-3, 1e-3), 'bad']:
        with pytest.raises(ValueError, match='time_range'):
            seq.calc_rf_power(time_range=bad_range)


def test_invalid_block_range(dummy_rf, dummy_sequence):
    seq = dummy_sequence
    seq.add_block(dummy_rf(amplitude=1.0, duration=1e-3))
    seq.add_block(dummy_rf(amplitude=1.0, duration=1e-3))

    for bad_range in [(0, 2), (2, 1), (1, 3), 'bad']:
        with pytest.raises(ValueError, match='block_range'):
            seq.calc_rf_power(block_range=bad_range)


def test_invalid_window_duration(dummy_rf, dummy_sequence):
    seq = dummy_sequence
    seq.add_block(dummy_rf(amplitude=1.0, duration=1e-3))

    for bad_window in [0.0, -1e-3]:
        with pytest.raises(ValueError, match='window_duration'):
            seq.calc_rf_power(window_duration=bad_window)


def test_sequence_no_rf(dummy_sequence):
    """A sequence without RF events has zero RF power."""
    seq = dummy_sequence
    seq.add_block(pp.make_delay(1e-3))

    mean_pwr, peak, rms, energy = seq.calc_rf_power()

    assert energy == 0
    assert mean_pwr == 0
    assert rms == 0
    assert peak == 0


def test_sequence_empty(dummy_sequence):
    """An empty sequence returns NaN for mean power and RMS, like MATLAB."""
    seq = dummy_sequence

    mean_pwr, _, rms, energy = seq.calc_rf_power()

    assert np.isnan(mean_pwr)
    assert np.isnan(rms)
    assert energy == 0
