"""
Tests for the ROTATIONS extension: make_rotation, rotation blocks in set_block/get_block,
rotation-aware gradient continuity checks, waveforms, write/read and plotting.
"""

from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pypulseq as pp
import pytest
from scipy.spatial.transform import Rotation


def rot_mat_from_quat(q):
    """Rotation matrix of a scalar-first quaternion [w, x, y, z]."""
    return Rotation.from_quat(np.asarray(q), scalar_first=True).as_matrix()


# ----------------------------------------------------------------------
# make_rotation
# ----------------------------------------------------------------------


def test_make_rotation_phi():
    phi = 1.234
    rot = pp.make_rotation(phi=phi)
    assert rot.type == 'rot3D'
    np.testing.assert_allclose(rot.rot_quaternion, [np.cos(phi / 2), 0, 0, np.sin(phi / 2)])


def test_make_rotation_phi_theta():
    phi, theta = 0.7, -0.4
    rot = pp.make_rotation(phi=phi, theta=theta)
    expected = (Rotation.from_euler('z', phi) * Rotation.from_euler('y', theta)).as_matrix()
    np.testing.assert_allclose(rot_mat_from_quat(rot.rot_quaternion), expected, atol=1e-12)


def test_make_rotation_axis_angle():
    angle = 0.9
    axis = np.array([1.0, -2.0, 2.0])
    rot = pp.make_rotation(axis=axis, angle=angle)
    axis = axis / np.linalg.norm(axis)
    np.testing.assert_allclose(rot.rot_quaternion, [np.cos(angle / 2), *(np.sin(angle / 2) * axis)])


def test_make_rotation_quaternion_normalized():
    q = np.array([2.0, 0.0, 0.0, 0.0])
    rot = pp.make_rotation(quaternion=q)
    np.testing.assert_allclose(rot.rot_quaternion, [1, 0, 0, 0])


def test_make_rotation_rot_matrix():
    rot_mat = Rotation.from_euler('yzx', [0.3, -0.5, 1.1]).as_matrix()
    rot = pp.make_rotation(rot_matrix=rot_mat)
    np.testing.assert_allclose(rot_mat_from_quat(rot.rot_quaternion), rot_mat, atol=1e-12)


@pytest.mark.parametrize(
    'kwargs_a,kwargs_b',
    [
        ({'phi': -np.pi / 2}, {'phi': 3 * np.pi / 2}),  # [-pi, pi) vs [0, 2pi) convention
        ({'phi': -0.3}, {'phi': 2 * np.pi - 0.3}),
        ({'phi': -np.pi}, {'phi': np.pi}),  # w == 0: sign decided by the first non-zero component
        ({'phi': 4.0, 'theta': 0.2}, {'phi': 4.0 - 2 * np.pi, 'theta': 0.2}),
        ({'axis': [0, 1, 0], 'angle': 0.5}, {'axis': [0, -1, 0], 'angle': -0.5}),
        ({'quaternion': [-0.5, 0.5, 0.5, 0.5]}, {'quaternion': [0.5, -0.5, -0.5, -0.5]}),
        ({'quaternion': [0, 0, 0, -1]}, {'quaternion': [0, 0, 0, 1]}),
        ({'phi': np.pi}, {'quaternion': [0, 0, 0, 1]}),
        ({'rot_matrix': Rotation.from_euler('z', 4.0).as_matrix()}, {'phi': 4.0}),
    ],
)
def test_make_rotation_canonical_sign(kwargs_a, kwargs_b):
    # Inputs describing the same rotation yield identical quaternions (w >= 0, or first non-zero component > 0)
    q_a = pp.make_rotation(**kwargs_a).rot_quaternion
    q_b = pp.make_rotation(**kwargs_b).rot_quaternion
    np.testing.assert_allclose(q_a, q_b, atol=1e-12)
    assert q_a[np.flatnonzero(np.abs(q_a) > 1e-12)[0]] > 0
    np.testing.assert_allclose(rot_mat_from_quat(q_a), rot_mat_from_quat(q_b), atol=1e-12)


def test_make_rotation_scipy_rotation():
    r = Rotation.from_euler('x', -0.6)
    rot = pp.make_rotation(rotation=r)
    np.testing.assert_allclose(rot_mat_from_quat(rot.rot_quaternion), r.as_matrix(), atol=1e-12)


def test_make_rotation_no_negative_zero():
    # Canonicalization must not introduce -0.0 components (they would show up as '-0' when formatted)
    q = pp.make_rotation(phi=3 * np.pi / 2).rot_quaternion
    assert not any(np.signbit(v) for v in q if v == 0)


def test_make_rotation_input_validation():
    with pytest.raises(ValueError):  # no input
        pp.make_rotation()
    with pytest.raises(ValueError):  # two input groups
        pp.make_rotation(phi=0.5, quaternion=[1, 0, 0, 0])
    with pytest.raises(ValueError):  # theta without phi
        pp.make_rotation(theta=0.5)
    with pytest.raises(ValueError):  # axis without angle
        pp.make_rotation(axis=[0, 0, 1])
    with pytest.raises(ValueError):  # angle without axis
        pp.make_rotation(angle=0.5)
    with pytest.raises(ValueError):  # zero axis
        pp.make_rotation(axis=[0, 0, 0], angle=0.5)
    with pytest.raises(ValueError):  # non-finite axis
        pp.make_rotation(axis=[0, np.nan, 1], angle=0.5)
    with pytest.raises(ValueError):  # wrong axis size
        pp.make_rotation(axis=[0, 0, 1, 0], angle=0.5)
    with pytest.raises(ValueError):  # phi too small
        pp.make_rotation(phi=-4.0)
    with pytest.raises(ValueError):  # phi too large
        pp.make_rotation(phi=2 * np.pi + 1e-3)
    with pytest.raises(ValueError):  # theta out of range
        pp.make_rotation(phi=0.5, theta=4.0)
    with pytest.raises(ValueError):  # angle out of range
        pp.make_rotation(axis=[0, 0, 1], angle=4.0)
    with pytest.raises(ValueError):  # non-finite angle
        pp.make_rotation(axis=[0, 0, 1], angle=np.nan)
    with pytest.raises(ValueError):  # zero quaternion
        pp.make_rotation(quaternion=[0, 0, 0, 0])
    with pytest.raises(ValueError):  # wrong quaternion size
        pp.make_rotation(quaternion=[1, 0, 0])
    with pytest.raises(ValueError):  # non-finite quaternion
        pp.make_rotation(quaternion=[1, 0, 0, np.nan])
    with pytest.raises(ValueError):  # improper rotation (det = -1)
        pp.make_rotation(rot_matrix=np.diag([1.0, 1.0, -1.0]))
    with pytest.raises(ValueError):  # non-orthogonal matrix
        pp.make_rotation(rot_matrix=2 * np.eye(3))
    with pytest.raises(ValueError):  # wrong matrix shape
        pp.make_rotation(rot_matrix=np.eye(4))
    with pytest.raises(ValueError):  # non-finite matrix
        pp.make_rotation(rot_matrix=np.full((3, 3), np.nan))
    with pytest.raises(ValueError):  # not a Rotation object
        pp.make_rotation(rotation='x')
    with pytest.raises(ValueError):  # rotation stack
        pp.make_rotation(rotation=Rotation.from_euler('z', [0.0, 1.0]))


# ----------------------------------------------------------------------
# set_block / get_block
# ----------------------------------------------------------------------


def test_rotation_block_roundtrip():
    seq = pp.Sequence()
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    gz = pp.make_trapezoid('z', area=-200, duration=1e-3)
    adc = pp.make_adc(num_samples=100, duration=1e-3)
    rot = pp.make_rotation(phi=0.5)

    seq.add_block(gx, gz, adc, rot)

    block = seq.get_block(1)
    assert block.rotation.type == 'rot3D'
    np.testing.assert_allclose(block.rotation.rot_quaternion, rot.rot_quaternion)


def test_two_rotations_per_block_raise():
    seq = pp.Sequence()
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    with pytest.raises(ValueError, match='Only one rotation extension'):
        seq.add_block(gx, pp.make_rotation(phi=0.5), pp.make_rotation(phi=0.7))


def test_rotation_deduplication():
    seq = pp.Sequence()
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    gy = pp.make_trapezoid('y', area=500, duration=1e-3)
    seq.add_block(gx, pp.make_rotation(phi=0.5))
    seq.add_block(gy, pp.make_rotation(phi=0.5))
    assert len(seq.rotation_library.data) == 1


def test_block_readd_keeps_rotation():
    seq = pp.Sequence()
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    rot = pp.make_rotation(phi=0.5)
    seq.add_block(gx, rot)
    seq.add_block(seq.get_block(1))

    assert len(seq.block_events) == 2
    block = seq.get_block(2)
    assert block.rotation.type == 'rot3D'
    np.testing.assert_allclose(block.rotation.rot_quaternion, rot.rot_quaternion)


# ----------------------------------------------------------------------
# waveforms
# ----------------------------------------------------------------------


def test_waveforms_match_rotated_events():
    phi = 0.8
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    gy = pp.make_trapezoid('y', area=-400, duration=1e-3)

    seq_rotated = pp.Sequence()
    seq_rotated.add_block(*pp.rotate(gx, gy, angle=phi, axis='z'))

    seq_extension = pp.Sequence()
    seq_extension.add_block(gx, gy, pp.make_rotation(phi=phi))

    wave_rotated = seq_rotated.waveforms()
    wave_extension = seq_extension.waveforms()
    for w_rot, w_ext in zip(wave_rotated, wave_extension, strict=True):
        np.testing.assert_allclose(w_rot, w_ext, atol=1e-6)


def test_waveforms_match_rotate_3d():
    axis = np.array([1.0, 1.0, 1.0])
    angle = 0.7
    rot = Rotation.from_rotvec(angle * axis / np.linalg.norm(axis)).as_matrix()
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    gy = pp.make_trapezoid('y', area=-400, duration=1e-3)
    gz = pp.make_trapezoid('z', area=200, duration=1e-3)

    seq_rotated = pp.Sequence()
    seq_rotated.add_block(*pp.rotate_3d(gx, gy, gz, rotation=rot))

    seq_extension = pp.Sequence()
    seq_extension.add_block(gx, gy, gz, pp.make_rotation(axis=axis, angle=angle))

    wave_rotated = seq_rotated.waveforms()
    wave_extension = seq_extension.waveforms()
    for w_rot, w_ext in zip(wave_rotated, wave_extension, strict=True):
        np.testing.assert_allclose(w_rot, w_ext, atol=1e-6)

    # get_block must return the unrotated (logical frame) gradient
    block = seq_extension.get_block(1)
    assert block.gx.type == 'trap'
    np.testing.assert_allclose(
        [block.gx.amplitude, block.gx.rise_time, block.gx.flat_time, block.gx.fall_time],
        [gx.amplitude, gx.rise_time, gx.flat_time, gx.fall_time],
    )


# ----------------------------------------------------------------------
# write / read
# ----------------------------------------------------------------------


def _make_rot_seq():
    seq = pp.Sequence()
    rf = pp.make_block_pulse(flip_angle=np.pi / 2, duration=1e-3, use='excitation')
    gx = pp.make_trapezoid('x', area=1000, duration=1e-3)
    adc = pp.make_adc(num_samples=100, duration=1e-3)
    seq.add_block(rf, gx, adc, pp.make_rotation(phi=0.5))
    seq.add_block(pp.make_delay(1e-3))
    seq.add_block(gx, pp.make_rotation(axis=[0, 0, 1], angle=-0.3))
    return seq


def test_write_read_rotation(tmp_path, compare_seq_file):
    seq = _make_rot_seq()
    file1 = tmp_path / 'rot1.seq'
    seq.write(file1)

    # RequiredExtensions definition must contain the ROTATIONS token
    assert 'ROTATIONS' in seq.get_definition('RequiredExtensions').split()

    # The written file must contain the extension section
    content = Path(file1).read_text()
    assert 'extension ROTATIONS' in content
    assert 'RequiredExtensions ROTATIONS' in content

    # Read back and compare
    seq2 = pp.Sequence()
    seq2.read(file1)
    assert len(seq2.rotation_library.data) == len(seq.rotation_library.data)
    for k in seq2.rotation_library.data:
        np.testing.assert_allclose(seq2.rotation_library.data[k], seq.rotation_library.data[k], atol=1e-6)

    block = seq2.get_block(1)
    assert block.rotation is not None
    np.testing.assert_allclose(block.rotation.rot_quaternion, seq.get_block(1).rotation.rot_quaternion, atol=1e-6)

    wave1, wave2 = seq.waveforms(), seq2.waveforms()
    for w1, w2 in zip(wave1, wave2, strict=True):
        np.testing.assert_allclose(w1, w2, rtol=1e-5)

    # Rewriting the read sequence produces an identical file
    file2 = tmp_path / 'rot2.seq'
    seq2.write(file2)
    compare_seq_file(file1, file2)


def test_read_normalizes_quaternion(tmp_path):
    seq = _make_rot_seq()
    file1 = tmp_path / 'rot.seq'
    seq.write(file1)

    # Scale the quaternion values on the ROTATIONS extension lines by 2
    lines = Path(file1).read_text().splitlines()
    idx = next(i for i, line in enumerate(lines) if line.startswith('extension ROTATIONS'))
    for i in range(idx + 1, len(lines)):
        if not lines[i] or lines[i].startswith(('[', '#', 'extension')):
            break
        tok = lines[i].split()
        lines[i] = ' '.join([tok[0], *[f'{2 * float(v):g}' for v in tok[1:]]])
    Path(file1).write_text('\n'.join(lines) + '\n')

    seq2 = pp.Sequence()
    seq2.read(file1)
    for k in seq2.rotation_library.data:
        q = np.asarray(seq2.rotation_library.data[k])
        assert np.linalg.norm(q) == pytest.approx(1.0)
        np.testing.assert_allclose(q, np.asarray(seq.rotation_library.data[k]), atol=1e-6)


def test_write_v141_raises(tmp_path):
    seq = _make_rot_seq()
    with pytest.raises(RuntimeError, match='rotation extension'):
        seq.write(tmp_path / 'rot.seq', v141_compat=True)


# ----------------------------------------------------------------------
# Gradient continuity checks
# ----------------------------------------------------------------------


def _arb_grad(channel, first, last, n=10):
    """Arbitrary gradient ramping from `first` to `last` (on the grad raster)."""
    return pp.make_arbitrary_grad(channel, np.linspace(first, last, n), first=first, last=last)


def test_grad_check_rotated_boundary():
    amp = 1e5

    # Block ending at non-zero amplitude on x, followed by a block rotated by
    # -90 deg about z whose (logical) y gradient starts with the same amplitude:
    # in the physical frame the boundary values match -> must pass.
    seq = pp.Sequence()
    seq.add_block(_arb_grad('x', 0.0, amp))
    seq.add_block(_arb_grad('y', amp, 0.0), pp.make_rotation(phi=-np.pi / 2))
    seq.add_block(_arb_grad('x', 0.0, 0.0))

    # Same without the rotation -> must fail
    seq = pp.Sequence()
    seq.add_block(_arb_grad('x', 0.0, amp))
    with pytest.raises(RuntimeError, match='same amplitude at the connection point'):
        seq.add_block(_arb_grad('y', amp, 0.0))

    # Same with the wrong rotation axis -> must fail
    seq = pp.Sequence()
    seq.add_block(_arb_grad('x', 0.0, amp))
    with pytest.raises(RuntimeError, match='same amplitude at the connection point'):
        seq.add_block(_arb_grad('y', amp, 0.0), pp.make_rotation(axis=[1, 0, 0], angle=np.pi / 2))


def test_grad_check_overwrite_block_with_rotated_neighbors():
    amp = 1e5
    rot_z90 = pp.make_rotation(phi=np.pi / 2)
    rot_zm90 = pp.make_rotation(phi=-np.pi / 2)

    seq = pp.Sequence()
    # Physical end of block 1: Rz(+90) @ (0, amp, 0) = (-amp, 0, 0)
    seq.add_block(_arb_grad('y', 0.0, amp), rot_z90)
    # Block 2 (unrotated) must start at -amp on x; ends at +amp on x
    seq.add_block(_arb_grad('x', -amp, amp))
    # Block 3 (rotated -90): logical start (0, amp, 0) matches R.T @ (amp, 0, 0) = (0, amp, 0)
    seq.add_block(_arb_grad('y', amp, 0.0), rot_zm90)

    # Overwriting block 2 with identical content must pass (boundary values of the
    # rotated neighbors are looked up through _boundary_grad_vals)
    seq.set_block(2, _arb_grad('x', -amp, amp))

    # Overwriting block 2 with a mismatching start amplitude must fail
    with pytest.raises(RuntimeError, match='same amplitude at the connection point'):
        seq.set_block(2, _arb_grad('x', amp, amp))


# ----------------------------------------------------------------------
# plot / timing / report
# ----------------------------------------------------------------------


def test_plot_and_report():
    seq = _make_rot_seq()
    with patch('matplotlib.pyplot.show'):
        seq.plot()
    plt.close('all')

    ok, error_report = seq.check_timing()
    assert ok, error_report

    assert isinstance(seq.test_report(), str)
