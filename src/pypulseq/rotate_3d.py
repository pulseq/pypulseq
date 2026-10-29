import numbers
from types import SimpleNamespace
from typing import List, Union

import numpy as np
from scipy.spatial.transform import Rotation

from pypulseq.add_gradients import add_gradients
from pypulseq.block_to_events import block_to_events
from pypulseq.make_delay import make_delay
from pypulseq.opts import Opts
from pypulseq.scale_grad import scale_grad
from pypulseq.utils.tracing import trace, trace_enabled


def _get_grad_abs_mag(grad: SimpleNamespace) -> float:
    if grad.type == 'trap':
        return abs(grad.amplitude)
    return np.max(np.abs(grad.waveform))


def _rotation_to_matrix(rotation: Union[np.ndarray, list, tuple, float, Rotation]) -> np.ndarray:
    """Convert the supported rotation representations to a 3x3 rotation matrix."""
    if isinstance(rotation, Rotation):
        return rotation.as_matrix()

    rotation = np.asarray(rotation, dtype=float)
    if rotation.shape == (3, 3):
        return rotation

    values = rotation.ravel()
    if values.size == 4:
        # Unit quaternion with the scalar component first [w, x, y, z]; scipy expects [x, y, z, w]
        return Rotation.from_quat(values[[1, 2, 3, 0]]).as_matrix()

    if values.size in (1, 2):
        phi = values[0]
        if not -np.pi <= phi < 2 * np.pi:
            raise ValueError(f'Rotation angle phi ({phi:.2f}) is invalid. Must be within [-pi, 2*pi) radians.')
        rot = Rotation.from_euler('z', phi)
        if values.size == 2:
            theta = values[1]
            if not -np.pi <= theta <= np.pi:
                raise ValueError(f'Rotation angle theta ({theta:.2f}) is invalid. Must be within [-pi, pi] radians.')
            rot = rot * Rotation.from_euler('y', theta)
        return rot.as_matrix()

    raise ValueError(
        'The parameter `rotation` must be a 3x3 rotation matrix, a unit quaternion [w, x, y, z], '
        'a pair of angles [phi, theta], a single angle phi or a scipy Rotation object.'
    )


def rotate_3d(
    *args: SimpleNamespace,
    rotation: Union[np.ndarray, list, tuple, float, Rotation],
    system: Union[Opts, None] = None,
) -> List[SimpleNamespace]:
    """
    Rotates the corresponding gradient(s) by the provided 3D rotation. Non-gradient events are not affected.

    The rotation can be given as a 3x3 rotation matrix, a unit quaternion with the scalar component first
    [w, x, y, z], a pair of polar angles [phi, theta] in radians (rotation about y by theta followed by a rotation
    about z by phi), a single angle phi in radians (rotation about z) or a `scipy.spatial.transform.Rotation` object.
    For a single rotation about the z-axis `pypulseq.rotate()` can be faster.

    Gradients that end up on the same axis are added together with `pypulseq.add_gradients()`. Rotated gradient
    components with negligible amplitude (relative to the largest input gradient) are dropped. Events that are not
    modified are passed through unchanged, including their `id` attribute if present; the `id` attribute is removed
    from the rotated gradients.

    See also `pypulseq.rotate.rotate()` and `pypulseq.Sequence.sequence.add_block()`.

    Parameters
    ----------
    args : SimpleNamespace
        Events to be rotated. Can also be a single block (as returned by `Sequence.get_block()`), in which case
        the block duration is preserved as a delay event, or a single list of events, e.g. the output of a
        previous `rotate_3d()` call.
    rotation : numpy.ndarray, list, tuple, float or scipy.spatial.transform.Rotation
        Rotation as a 3x3 matrix, a unit quaternion [w, x, y, z], polar angles [phi, theta] in radians, a single
        angle phi in radians or a scipy `Rotation` object.
    system : Opts, default=Opts()
        System limits, used when adding gradient components on the same axis.

    Returns
    -------
    rotated_grads : list[SimpleNamespace]
        Non-gradient events followed by the rotated gradient(s). Can be passed directly to `Sequence.add_block()`.

    Raises
    ------
    ValueError
        If `rotation` has an unsupported shape or the angles are out of range.
        If more than one gradient is provided for the same axis.
        If a gradient has an invalid channel.
    """
    if system is None:
        system = Opts.default

    rot_mat = _rotation_to_matrix(rotation)

    events = block_to_events(*args)
    if len(events) == 1 and isinstance(events[0], (list, tuple)):
        events = events[0]  # Allow nested calls, e.g. rotate_3d(rotate_3d(...), rotation=...)

    axes = ['x', 'y', 'z']
    grads_in = {}
    bypass = []
    for event in events:
        if event is None:
            continue
        if isinstance(event, numbers.Real):
            bypass.append(make_delay(event))  # block_duration of a block passed as input
            continue
        if getattr(event, 'type', None) not in ('grad', 'trap'):
            bypass.append(event)
        elif event.channel not in axes:
            raise ValueError(f'Invalid gradient channel. Expected one of {axes}, got: {event.channel}')
        elif event.channel in grads_in:
            raise ValueError(f'More than one gradient on the same axis {event.channel} provided.')
        else:
            grads_in[event.channel] = event

    # Measure of the relevant amplitude, used to drop negligible rotated components
    fthresh = 1e-6
    thresh = fthresh * max((_get_grad_abs_mag(g) for g in grads_in.values()), default=0)

    rotated_grads = []
    for j, out_axis in enumerate(axes):
        components = []
        for i, in_axis in enumerate(axes):
            if in_axis not in grads_in or abs(rot_mat[j, i]) < fthresh:
                continue
            g = scale_grad(grad=grads_in[in_axis], scale=rot_mat[j, i])
            g.channel = out_axis
            components.append(g)

        if len(components) == 0:
            continue
        g = components[0] if len(components) == 1 else add_gradients(grads=components, system=system)
        if _get_grad_abs_mag(g) >= thresh:
            rotated_grads.append(g)

    if trace_enabled():
        for grad in rotated_grads:
            grad.trace = trace()

    return [*bypass, *rotated_grads]
