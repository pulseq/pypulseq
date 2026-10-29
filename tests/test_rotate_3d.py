from copy import copy
from types import SimpleNamespace

import numpy as np
import pypulseq as pp
import pytest
from pypulseq.rotate_3d import rotate_3d
from scipy.spatial.transform import Rotation

gx_trap = pp.make_trapezoid('x', area=1000, duration=1e-3)
gy_trap = pp.make_trapezoid('y', area=-500, duration=1e-3)
gz_trap = pp.make_trapezoid('z', area=200, duration=1e-3)
gx_ext = pp.make_extended_trapezoid('x', amplitudes=[0, 1e5, 5e4, 0], times=[0, 1e-4, 4e-4, 6e-4])
gy_arb = pp.make_extended_trapezoid(
    'y', amplitudes=[0, 8e4, 8e4, 0], times=[0, 1e-4, 3e-4, 5e-4], convert_to_arbitrary=True
)
rf = pp.make_block_pulse(flip_angle=np.pi / 2, duration=1e-3)
adc = pp.make_adc(num_samples=100, duration=1e-3)
delay = pp.make_delay(2e-3)


def rot_mat(axis, angle):  # standard right-handed rotation matrix
    return Rotation.from_euler(axis, angle).as_matrix()


def by_channel(events):  # dict channel -> gradient, non-gradients ignored
    grads = {}
    for event in events:
        if getattr(event, 'type', None) not in ('grad', 'trap'):
            continue
        assert event.channel not in grads, f'duplicate channel {event.channel}'
        grads[event.channel] = event
    return grads


def assert_grads_equal(a, b, rel=1e-6, abs_=1e-6):
    assert a.type == b.type and a.channel == b.channel
    assert a.delay == pytest.approx(b.delay, abs=1e-12)
    assert a.area == pytest.approx(b.area, rel=rel, abs=abs_)
    if a.type == 'trap':
        for f in ('amplitude', 'rise_time', 'flat_time', 'fall_time', 'flat_area'):
            assert getattr(a, f) == pytest.approx(getattr(b, f), rel=rel, abs=abs_), f
    else:
        np.testing.assert_allclose(a.tt, b.tt, rtol=rel, atol=1e-12)
        np.testing.assert_allclose(a.waveform, b.waveform, rtol=rel, atol=abs_)
        for f in ('first', 'last', 'shape_dur'):
            assert getattr(a, f) == pytest.approx(getattr(b, f), rel=rel, abs=abs_), f


def test_identity():
    amplitudes = [g.amplitude for g in (gx_trap, gy_trap, gz_trap)]
    out = rotate_3d(gx_trap, gy_trap, gz_trap, rotation=np.eye(3))
    assert len(out) == 3
    assert [g.channel for g in out] == ['x', 'y', 'z']
    for rotated, original in zip(out, (gx_trap, gy_trap, gz_trap)):
        assert_grads_equal(rotated, original)
    # Inputs are not modified
    assert [g.amplitude for g in (gx_trap, gy_trap, gz_trap)] == amplitudes


def test_rotation_input_forms():
    angle = 0.3
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    reference = by_channel(rotate_3d(gx_trap, gy_trap, rotation=rot_mat('z', angle)))
    forms = [
        rot_mat('z', angle),
        [c, 0, 0, s],
        angle,
        [angle, 0.0],
        Rotation.from_euler('z', angle),
    ]
    for r in forms:
        out = by_channel(rotate_3d(gx_trap, gy_trap, rotation=r))
        assert set(out) == set(reference)
        for channel in reference:
            assert_grads_equal(out[channel], reference[channel])

    # [phi, theta] means rotation about y by theta followed by rotation about z by phi
    phi, theta = 0.4, -0.7
    a = by_channel(rotate_3d(gx_trap, gy_trap, gz_trap, rotation=[phi, theta]))
    b = by_channel(rotate_3d(gx_trap, gy_trap, gz_trap, rotation=rot_mat('z', phi) @ rot_mat('y', theta)))
    assert set(a) == set(b)
    for channel in a:
        assert_grads_equal(a[channel], b[channel])


@pytest.mark.parametrize('grads', [[g] for g in (gx_trap, gy_trap, gz_trap, gx_ext, gy_arb)] + [[gx_trap, gy_trap]])
@pytest.mark.parametrize('angle', [0.0, 0.1, np.pi / 4, np.pi / 2, np.pi, -np.pi / 3, 2 * np.pi])
@pytest.mark.parametrize('axis', ['x', 'y', 'z'])
def test_vs_rotate(axis, angle, grads):
    # pypulseq's rotate() about the y-axis currently uses the opposite handedness to a right-handed
    # rotation matrix (and to MATLAB's mr.rotate); compensate here until rotate() is fixed
    sign = -1 if axis == 'y' else 1
    a = by_channel(pp.rotate(*[copy(g) for g in grads], axis=axis, angle=angle))
    b = by_channel(rotate_3d(*grads, rotation=rot_mat(axis, sign * angle)))
    assert set(a) == set(b)
    for channel in a:
        assert_grads_equal(a[channel], b[channel], rel=1e-6, abs_=1e-3)


def test_oblique_rotation_and_inverse():
    R = Rotation.from_rotvec(0.9 * np.array([0.5, 0.3, 0.6]) / np.linalg.norm([0.5, 0.3, 0.6])).as_matrix()
    out = rotate_3d(gx_trap, rotation=R)
    assert len(out) == 3
    assert all(g.type == 'trap' for g in out)
    for j, g in enumerate(out):
        assert g.amplitude == pytest.approx(R[j, 0] * gx_trap.amplitude, rel=1e-9)
    amplitudes = np.array([g.amplitude for g in out])
    assert np.linalg.norm(amplitudes) == pytest.approx(abs(gx_trap.amplitude), rel=1e-9)

    back = rotate_3d(*out, rotation=R.T)
    assert len(back) == 1
    assert back[0].channel == 'x'
    assert_grads_equal(back[0], gx_trap, abs_=1e-6 * abs(gx_trap.amplitude))


def test_passthrough_and_id():
    gx = copy(gx_trap)
    gx.id = 7
    rf_with_id = copy(rf)
    rf_with_id.id = 3
    out = rotate_3d(rf_with_id, gx, adc, delay, rotation=rot_mat('z', np.pi / 2))
    assert out[0] is rf_with_id
    assert out[1] is adc
    assert out[2] is delay
    assert out[3].type == 'trap' and out[3].channel == 'y'
    assert not hasattr(out[3], 'id')
    assert rf_with_id.id == 3
    assert gx.id == 7
    assert gx.amplitude == pytest.approx(gx_trap.amplitude)

    # None inputs are skipped
    out = rotate_3d(None, gx_trap, None, rotation=np.eye(3))
    assert len(out) == 1
    assert_grads_equal(out[0], gx_trap)


def test_block_and_list_input():
    seq = pp.Sequence()
    seq.add_block(rf, gx_trap, gz_trap, adc)
    block = seq.get_block(1)
    R = rot_mat('x', 0.5)
    a = rotate_3d(block, rotation=R)
    b = rotate_3d(rf, gx_trap, gz_trap, adc, rotation=R)
    # The block's block_duration is preserved as one extra delay event
    delays = [e for e in a if getattr(e, 'type', None) == 'delay']
    assert len(delays) == 1
    assert delays[0].delay == pytest.approx(block.block_duration)
    a_events = [e for e in a if getattr(e, 'type', None) != 'delay']
    assert len(a_events) == len(b)
    assert sorted(getattr(e, 'type', None) for e in a_events) == sorted(getattr(e, 'type', None) for e in b)
    ca, cb = by_channel(a_events), by_channel(b)
    assert set(ca) == set(cb)
    for channel in ca:
        assert_grads_equal(ca[channel], cb[channel])

    # The rotated block can be passed back to add_block with the block duration preserved
    seq2 = pp.Sequence()
    seq2.add_block(*a)
    assert seq2.block_durations[1] == pytest.approx(seq.block_durations[1])

    # Nested call: rotate_3d(rotate_3d(...)) equals a single rotation by the combined matrix
    R1, R2 = rot_mat('z', 0.4), rot_mat('x', -0.6)
    nested = by_channel(rotate_3d(rotate_3d(gx_trap, gy_trap, rotation=R1), rotation=R2))
    combined = by_channel(rotate_3d(gx_trap, gy_trap, rotation=R2 @ R1))
    assert set(nested) == set(combined)
    for channel in nested:
        assert_grads_equal(nested[channel], combined[channel], abs_=1e-3)


def test_negligible_components_dropped():
    out = rotate_3d(gx_trap, rotation=rot_mat('z', np.pi / 2))
    assert len(out) == 1
    assert out[0].channel == 'y'
    assert out[0].amplitude == pytest.approx(gx_trap.amplitude, rel=1e-6)

    out = rotate_3d(gz_trap, rotation=rot_mat('z', 1.0))
    assert len(out) == 1
    assert_grads_equal(out[0], gz_trap)


def test_mixed_gradient_types_are_added():
    out = by_channel(rotate_3d(gx_trap, gy_arb, rotation=rot_mat('z', np.pi / 4)))
    assert set(out) == {'x', 'y'}
    assert out['x'].type == 'grad' and out['y'].type == 'grad'
    c, s = np.cos(np.pi / 4), np.sin(np.pi / 4)
    assert out['x'].area == pytest.approx(c * gx_trap.area - s * gy_arb.area, rel=1e-6)
    assert out['y'].area == pytest.approx(s * gx_trap.area + c * gy_arb.area, rel=1e-6)


def test_errors():
    with pytest.raises(ValueError, match='More than one gradient'):
        rotate_3d(gx_trap, gx_trap, rotation=np.eye(3))
    with pytest.raises(ValueError, match='Invalid gradient channel'):
        rotate_3d(SimpleNamespace(type='trap', channel='q'), rotation=np.eye(3))
    with pytest.raises(ValueError, match='must be a 3x3'):
        rotate_3d(gx_trap, rotation=[1, 2, 3])
    with pytest.raises(ValueError, match='phi'):
        rotate_3d(gx_trap, rotation=7.0)
    with pytest.raises(ValueError, match='theta'):
        rotate_3d(gx_trap, rotation=[0.0, 4.0])
    with pytest.raises(TypeError):
        rotate_3d(gx_trap)


def test_system_limits_enforced():
    sys_lim = pp.Opts(max_grad=30, grad_unit='mT/m', max_slew=100, slew_unit='T/m/s')
    g1 = pp.make_trapezoid('x', amplitude=0.95 * sys_lim.max_grad, flat_time=1e-3, system=sys_lim)
    g2 = pp.make_trapezoid('y', amplitude=0.95 * sys_lim.max_grad, flat_time=1e-3, system=sys_lim)
    with pytest.raises(ValueError):
        rotate_3d(g1, g2, rotation=rot_mat('z', np.pi / 4), system=sys_lim)

    sys_relaxed = pp.Opts(max_grad=60, grad_unit='mT/m', max_slew=200, slew_unit='T/m/s')
    out = rotate_3d(g1, g2, rotation=rot_mat('z', np.pi / 4), system=sys_relaxed)
    # Equal-amplitude x and y gradients rotated by 45 deg cancel on x and combine to a single y trap
    assert len(out) == 1
    assert out[0].type == 'trap' and out[0].channel == 'y'
    assert out[0].amplitude == pytest.approx(np.sqrt(2) * 0.95 * sys_lim.max_grad, rel=1e-6)


def test_sequence_roundtrip(tmp_path):
    seq = pp.Sequence()
    seq.add_block(*rotate_3d(rf, gx_trap, gz_trap, rotation=Rotation.from_rotvec([0.3, -0.2, 0.5]).as_matrix()))
    seq.add_block(*rotate_3d(gx_ext, gy_arb, rotation=[0.7, 0.2]))
    ok, err = seq.check_timing()
    assert ok, err

    seq.write(tmp_path / 'r.seq')
    seq2 = pp.Sequence()
    seq2.read(tmp_path / 'r.seq')

    for block_id in (1, 2):
        a = by_channel(list(seq.get_block(block_id).__dict__.values()))
        b = by_channel(list(seq2.get_block(block_id).__dict__.values()))
        assert set(a) == set(b)
        for channel in a:
            assert a[channel].area == pytest.approx(b[channel].area, rel=1e-4)
