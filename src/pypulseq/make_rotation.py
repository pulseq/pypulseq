from types import SimpleNamespace
from typing import Union

import numpy as np
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation

from pypulseq.utils.tracing import trace, trace_enabled


def _quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton product q1 * q2 of quaternions in scalar-first format [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ]
    )


def _scalar_first_quat(rotation: Rotation) -> np.ndarray:
    """Unit quaternion [w, x, y, z] of a scipy rotation in canonical form (w >= 0), see `Rotation.as_quat`."""
    q_xyzw = rotation.as_quat()  # canonical=True requires scipy >= 1.11
    if q_xyzw[3] < 0:
        q_xyzw = -q_xyzw
    return q_xyzw[[3, 0, 1, 2]]


def make_rotation(
    *,
    phi: Union[float, None] = None,
    theta: Union[float, None] = None,
    axis: Union[ArrayLike, None] = None,
    angle: Union[float, None] = None,
    quaternion: Union[ArrayLike, None] = None,
    rot_matrix: Union[ArrayLike, None] = None,
    rotation: Union[Rotation, None] = None,
) -> SimpleNamespace:
    """
    Create a rotation extension event (`type='rot3D'`).

    The returned object stores a unit quaternion in scalar-first format [w, x, y, z] in `rot_quaternion`. It can be
    attached to a sequence block; the gradients of that block are then interpreted in the logical (rotated) frame,
    while e.g. `Sequence.waveforms()` and `Sequence.plot()` report them in the physical scanner frame. Exactly one
    of the input forms below must be provided. Only one rotation extension is allowed per block.

    Parameters
    ----------
    phi : float, optional
        Rotation angle about the z axis in radians, in [-pi, 2*pi).
    theta : float, optional
        Rotation angle about the y axis in radians, in [-pi, pi]. Only valid together with `phi`; the combined
        rotation is Rz(phi) * Ry(theta) (rotation about y by `theta` followed by rotation about z by `phi`).
    axis : ArrayLike, optional
        Rotation axis as a 3-element vector (need not be normalized). Must be given together with `angle`.
    angle : float, optional
        Rotation angle about `axis` in radians, in [-pi, pi]. Must be given together with `axis`.
    quaternion : ArrayLike, optional
        Rotation as a 4-element quaternion [w, x, y, z]. Normalized before being stored.
    rot_matrix : ArrayLike, optional
        Rotation as a 3x3 matrix. Must be a proper rotation (orthogonal with determinant +1).
    rotation : scipy.spatial.transform.Rotation, optional
        Rotation as a single scipy `Rotation` object.

    Returns
    -------
    rot : SimpleNamespace
        Rotation extension event with `type == 'rot3D'` and `rot_quaternion` a unit quaternion [w, x, y, z].

    Raises
    ------
    ValueError
        If not exactly one rotation input form is provided, if `theta` is given without `phi`, if `axis`/`angle`
        are not given together, or if the provided rotation parameters are out of range or invalid.

    Examples
    --------
    >>> rot = pp.make_rotation(phi=np.pi / 3)
    >>> seq.add_block(gx, adc, rot)

    See Also
    --------
    pypulseq.rotate_3d.rotate_3d : Rotate gradient events explicitly.
    pypulseq.Sequence.sequence.Sequence.add_block : Add blocks to a sequence.
    """
    groups = [
        phi is not None or theta is not None,
        axis is not None or angle is not None,
        quaternion is not None,
        rot_matrix is not None,
        rotation is not None,
    ]

    if sum(groups) != 1:
        raise ValueError(
            'make_rotation requires exactly one rotation input (phi/theta, axis/angle, quaternion, '
            'rot_matrix or rotation).'
        )

    rot = SimpleNamespace(type='rot3D')

    if phi is not None or theta is not None:
        if phi is None:
            raise ValueError('Parameter theta is only valid together with phi.')
        phi = float(np.asarray(phi, dtype=float))
        theta = 0.0 if theta is None else float(np.asarray(theta, dtype=float))
        if not -np.pi <= phi < 2 * np.pi:
            raise ValueError(f'Rotation angle phi ({phi:.2f}) is invalid. Must be within [-pi, 2*pi) radians.')
        if not -np.pi <= theta <= np.pi:
            raise ValueError(f'Rotation angle theta ({theta:.2f}) is invalid. Must be within [-pi, pi] radians.')

        qz = np.array([np.cos(phi / 2), 0.0, 0.0, np.sin(phi / 2)])
        qy = np.array([np.cos(theta / 2), 0.0, np.sin(theta / 2), 0.0])
        rot.rot_quaternion = _quat_multiply(qz, qy)

    elif axis is not None or angle is not None:
        if axis is None or angle is None:
            raise ValueError('Parameters axis and angle must be provided together.')
        axis = np.asarray(axis, dtype=float).ravel()
        if axis.size != 3:
            raise ValueError(f'Rotation axis must be a 3-element vector, got {axis.size} elements.')
        norm = np.linalg.norm(axis)
        if not np.isfinite(norm) or norm == 0:
            raise ValueError('Rotation axis must have a finite, non-zero norm.')
        angle = float(np.asarray(angle, dtype=float))
        if not np.isfinite(angle) or not -np.pi <= angle <= np.pi:
            raise ValueError(f'Rotation angle ({angle:.2f}) is invalid. Must be within [-pi, pi] radians.')

        rot.rot_quaternion = np.array([np.cos(angle / 2), *(np.sin(angle / 2) * axis / norm)])

    elif quaternion is not None:
        q = np.asarray(quaternion, dtype=float).ravel()
        if q.size != 4:
            raise ValueError(f'Rotation quaternion must have 4 elements, got {q.size}.')
        if not np.all(np.isfinite(q)):
            raise ValueError('Rotation quaternion must only contain finite values.')
        norm = np.linalg.norm(q)
        if norm == 0:
            raise ValueError('Rotation quaternion norm must be non-zero.')
        rot.rot_quaternion = q / norm

    elif rot_matrix is not None:
        R = np.asarray(rot_matrix, dtype=float)
        if R.shape != (3, 3):
            raise ValueError(f'Rotation matrix must have shape (3, 3), got {R.shape}.')
        if not np.all(np.isfinite(R)):
            raise ValueError('Rotation matrix must only contain finite values.')
        if not np.allclose(R @ R.T, np.eye(3), atol=1e-6):
            raise ValueError('Rotation matrix is not orthogonal.')
        if not np.isclose(np.linalg.det(R), 1.0, atol=1e-6):
            raise ValueError('Rotation matrix determinant is not +1 (improper rotation).')
        rot.rot_quaternion = _scalar_first_quat(Rotation.from_matrix(R))

    else:
        if not isinstance(rotation, Rotation) or not rotation.single:
            raise ValueError('Parameter rotation must be a single scipy.spatial.transform.Rotation object.')
        rot.rot_quaternion = _scalar_first_quat(rotation)

    if trace_enabled():
        rot.trace = trace()

    return rot
