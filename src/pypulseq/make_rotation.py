from types import SimpleNamespace

import numpy as np
from numpy.typing import ArrayLike
from scipy.spatial.transform import Rotation

from pypulseq.utils.tracing import trace, trace_enabled


def _canonical_quat(q: np.ndarray) -> np.ndarray:
    """
    Return the quaternion with canonical sign (first non-zero component > 0).

    q and -q describe the same rotation. Choosing a fixed sign ensures that equivalent inputs (e.g. phi and
    phi - 2*pi, or a quaternion and its negative) yield identical rotation events, library entries and .seq files.
    Near-zero components are snapped to 0 first (e.g. cos(pi/2) ~ 6e-17), so that e.g. phi = pi and phi = -pi
    yield the same sign.

    Parameters
    ----------
    q : np.ndarray
        Unit quaternion [w, x, y, z].

    Returns
    -------
    np.ndarray
        Copy of `q` with canonical sign and without -0.0 components.
    """
    q = q.copy()
    q[np.abs(q) < 1e-12] = 0.0
    if q[np.flatnonzero(q)[0]] < 0:
        q = -q
        q += 0.0  # normalize -0.0 introduced by the negation to 0.0
    return q


def make_rotation(
    *,
    phi: float | None = None,
    theta: float | None = None,
    axis: ArrayLike | None = None,
    angle: float | None = None,
    quaternion: ArrayLike | None = None,
    rot_matrix: ArrayLike | None = None,
    rotation: Rotation | None = None,
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
        Rotation extension event with `type == 'rot3D'` and `rot_quaternion` a unit quaternion [w, x, y, z] in
        canonical form (w >= 0), i.e. inputs describing the same rotation yield the same quaternion.

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

    if phi is not None or theta is not None:
        if phi is None:
            raise ValueError('Parameter theta is only valid together with phi.')
        phi = float(np.asarray(phi, dtype=float))
        theta = 0.0 if theta is None else float(np.asarray(theta, dtype=float))
        if not -np.pi <= phi < 2 * np.pi:
            raise ValueError(f'Rotation angle phi ({phi:.2f}) is invalid. Must be within [-pi, 2*pi) radians.')
        if not -np.pi <= theta <= np.pi:
            raise ValueError(f'Rotation angle theta ({theta:.2f}) is invalid. Must be within [-pi, pi] radians.')

        rotation = Rotation.from_euler('z', phi) * Rotation.from_euler('y', theta)

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

        rotation = Rotation.from_rotvec(angle * axis / norm)

    elif quaternion is not None:
        q = np.asarray(quaternion, dtype=float).ravel()
        if q.size != 4:
            raise ValueError(f'Rotation quaternion must have 4 elements, got {q.size}.')
        if not np.all(np.isfinite(q)):
            raise ValueError('Rotation quaternion must only contain finite values.')
        norm = np.linalg.norm(q)
        if norm == 0:
            raise ValueError('Rotation quaternion norm must be non-zero.')
        rotation = Rotation.from_quat(q, scalar_first=True)

    elif rot_matrix is not None:
        rot_mat = np.asarray(rot_matrix, dtype=float)
        if rot_mat.shape != (3, 3):
            raise ValueError(f'Rotation matrix must have shape (3, 3), got {rot_mat.shape}.')
        if not np.all(np.isfinite(rot_mat)):
            raise ValueError('Rotation matrix must only contain finite values.')
        if not np.allclose(rot_mat @ rot_mat.T, np.eye(3), atol=1e-6):
            raise ValueError('Rotation matrix is not orthogonal.')
        if not np.isclose(np.linalg.det(rot_mat), 1.0, atol=1e-6):
            raise ValueError('Rotation matrix determinant is not +1 (improper rotation).')
        rotation = Rotation.from_matrix(rot_mat)

    else:
        if not isinstance(rotation, Rotation) or not rotation.single:
            raise ValueError('Parameter rotation must be a single scipy.spatial.transform.Rotation object.')

    rot = SimpleNamespace(type='rot3D')
    rot.rot_quaternion = _canonical_quat(rotation.as_quat(scalar_first=True))

    if trace_enabled():
        rot.trace = trace()

    return rot
