import numpy as np
import pypulseq as pp
import pytest
from pypulseq.Sequence.sequence import Sequence


def _eig_min(a: np.ndarray) -> float:
    # a is 3x3 symmetric
    w = np.linalg.eigvalsh(0.5 * (a + a.T))
    return float(w.min())


def test_calc_moments_btensor_basic_properties():
    sys = pp.Opts(
        max_grad=np.inf,
        grad_unit='mT/m',
        max_slew=np.inf,
        slew_unit='T/m/s',
        rf_dead_time=0.0,
        rf_ringdown_time=0.0,
        adc_dead_time=0.0,
    )

    # Initialize sequence
    seq = Sequence(sys)

    # ---- Define RF pulses ----
    rf90 = pp.make_block_pulse(
        flip_angle=np.pi / 2,
        duration=1e-3,
        system=sys,
        use='excitation',
    )
    rf180 = pp.make_block_pulse(
        flip_angle=np.pi,
        duration=1e-3,
        system=sys,
        use='refocusing',
    )

    # ---- Target b-value ----
    b_target = 1000e6  # s/m^2 (1000 s/mm^2)
    te = 80e-3
    delta = 10e-3
    big_delta = 30e-3

    # Solve for gradient amplitude (Hz/m)
    g_amp = (1 / (2 * np.pi * delta)) * np.sqrt(b_target / (big_delta - delta / 3))

    # ---- Define diffusion gradient lobe ----
    g = pp.make_trapezoid(channel='z', system=sys, amplitude=g_amp, duration=delta)

    # Interval between end of first g and start of second g
    diff_wait = big_delta - delta

    # RF free time
    rf_free_time = diff_wait - pp.calc_duration(rf180)

    # Center rf between gradient lobes
    g_wait = 0.5 * rf_free_time

    # Compute te wait time
    texc_center, _ = pp.calc_rf_center(rf90)
    tref_center, _ = pp.calc_rf_center(rf180)
    t_ref = (pp.calc_duration(rf90) - texc_center) + pp.calc_duration(g) + g_wait + tref_center
    te_wait = 0.5 * te - t_ref

    # ---- Timing layout ----
    seq.add_block(rf90)
    seq.add_block(pp.make_delay(te_wait))
    seq.add_block(g)
    seq.add_block(pp.make_delay(g_wait))
    seq.add_block(rf180)
    seq.add_block(pp.make_delay(g_wait))
    seq.add_block(g)

    # ---- Run moment calculator -----
    b, m1, m2, m3 = seq.calc_moments_btensor()  # defaults: calc_b=True only

    # Shapes
    assert b.ndim == 3 and b.shape[1:] == (3, 3)
    assert m1.shape == (b.shape[0], 3)
    assert m2.shape == (b.shape[0], 3)
    assert m3.shape == (b.shape[0], 3)

    # With one excitation this will be one repetition
    assert b.shape[0] == 1

    # Symmetry
    assert np.allclose(b[0], b[0].T, atol=1e-10, rtol=1e-10)

    # PSD (numerical tolerance)
    assert _eig_min(b[0]) > -1e-8

    # Dominant zz for z-only gradients
    assert b[0, 2, 2] >= b[0, 0, 0]
    assert b[0, 2, 2] >= b[0, 1, 1]

    # Off-diagonals should be ~0 (no x/y)
    assert abs(b[0, 0, 1]) < 1e-10
    assert abs(b[0, 0, 2]) < 1e-10
    assert abs(b[0, 1, 2]) < 1e-10

    # Convert to physical b-value (radians)
    b_measured = b[0, 2, 2]

    # Check quantitative correctness
    assert np.isclose(b_measured, b_target, rtol=1e-2)


def test_calc_moments_btensor_moments_flags():
    sys = pp.Opts(max_grad=28, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s')
    seq = Sequence(sys)

    rf90 = pp.make_block_pulse(np.pi / 2, duration=1e-3, system=sys, use='excitation')
    rf180 = pp.make_block_pulse(np.pi, duration=1e-3, system=sys, use='refocusing')

    g = pp.make_trapezoid(channel='x', system=sys, amplitude=5e-3 * sys.gamma, duration=5e-3)

    seq.add_block(rf90)
    seq.add_block(pp.make_delay(2e-3))
    seq.add_block(g)
    seq.add_block(pp.make_delay(1e-3))
    seq.add_block(rf180)
    seq.add_block(pp.make_delay(1e-3))
    seq.add_block(g)
    seq.add_block(pp.make_delay(2e-3))

    b, m1, m2, m3 = seq.calc_moments_btensor(True, True, True, True, 0)

    assert np.any(np.abs(m1) > 0)
    assert np.any(np.abs(m2) > 0)
    assert np.any(np.abs(m3) > 0)
    assert np.any(np.abs(b) > 0)


def _make_dwi_se(n_rep: int = 2) -> Sequence:
    # Diffusion-weighted spin echo with realistic (slow) slew rate and oblique diffusion directions
    sys = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=50, slew_unit='T/m/s')
    seq = Sequence(sys)

    rf90 = pp.make_block_pulse(flip_angle=np.pi / 2, duration=1e-3, system=sys, use='excitation')
    rf180 = pp.make_block_pulse(flip_angle=np.pi, duration=1e-3, system=sys, use='refocusing')

    te = 90e-3
    delta = 15e-3
    big_delta = 35e-3
    directions = [(0.48, 0.6, 0.64), (0, 0.6, 0.8)]

    for r in range(n_rep):
        g = [
            pp.make_trapezoid(channel=ch, system=sys, amplitude=0.8 * sys.max_grad * c, duration=delta)
            for ch, c in zip('xyz', directions[r % len(directions)], strict=True)
            if c != 0
        ]

        g_wait = 0.5 * (big_delta - delta - pp.calc_duration(rf180))
        g_wait = np.ceil(g_wait / sys.grad_raster_time) * sys.grad_raster_time
        t_ref = (
            (pp.calc_duration(rf90) - pp.calc_rf_center(rf90)[0] - rf90.delay)
            + pp.calc_duration(*g)
            + g_wait
            + pp.calc_rf_center(rf180)[0]
            + rf180.delay
        )
        te_wait = np.round((0.5 * te - t_ref) / sys.grad_raster_time) * sys.grad_raster_time

        seq.add_block(rf90)
        seq.add_block(pp.make_delay(te_wait))
        seq.add_block(*g)
        seq.add_block(pp.make_delay(g_wait))
        seq.add_block(rf180)
        seq.add_block(pp.make_delay(g_wait))
        seq.add_block(*g)
        seq.add_block(pp.make_delay(te + 20e-3))

    return seq


def test_calc_moments_btensor_matlab_reference(tmp_path):
    # Write and read back, so that the gradients match those seen by MATLAB
    filename = tmp_path / 'dwi_se.seq'
    _make_dwi_se().write(str(filename))
    seq = Sequence()
    seq.read(filename)

    b, m1, m2, m3 = seq.calc_moments_btensor(True, True, True, True, 0)

    # Reference values from MATLAB Pulseq calcMomentsBtensor(true, true, true, true, 0) on the same .seq file
    b_ref = [
        [
            [422059133.166, 522481719.621, 555496769.027],
            [522481719.621, 646803061.202, 687675504.089],
            [555496769.027, 687675504.089, 731131326.098],
        ],
        [
            [0, 0, 0],
            [0, 646803061.202, 847769380.931],
            [0, 847769380.931, 1111201330.59],
        ],
    ]
    m1_ref = [[-4138.46605624, -5119.14483084, -5441.2202917], [0, -5119.14483084, -6700.80524078]]
    m2_ref = [[-372.461945061, -460.723034776, -489.709826253], [0, -460.723034776, -603.07247167]]
    m3_ref = [[-26.62437566, -32.9282415702, -34.9981264325], [0, -32.9282415702, -43.0905704111]]

    np.testing.assert_allclose(b, b_ref, rtol=1e-9, atol=1e-3)
    np.testing.assert_allclose(m1, m1_ref, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(m2, m2_ref, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(m3, m3_ref, rtol=1e-9, atol=1e-9)


def test_calc_moments_btensor_n_dummy():
    seq = _make_dwi_se(n_rep=3)

    ref = seq.calc_moments_btensor(True, True, True, True, 0)
    out = seq.calc_moments_btensor(True, True, True, True, 1)

    # Skipping dummy scans must return the remaining repetitions unchanged
    for o, r in zip(out, ref, strict=True):
        assert o.shape[0] == 2
        np.testing.assert_allclose(o, r[1:])

    with pytest.raises(ValueError):
        seq.calc_moments_btensor(n_dummy=3)
