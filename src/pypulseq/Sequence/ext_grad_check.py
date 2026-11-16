from __future__ import annotations

import typing

import numpy as np
from scipy.spatial.transform import Rotation

from pypulseq import eps

if typing.TYPE_CHECKING:
    from pypulseq.Sequence.sequence import Sequence


def ext_grad_check(self: Sequence, block_index: int, check_g: dict, duration: float, rot_quaternion=None) -> None:
    """
    Check the gradient continuity of the block that is being set at `block_index`.

    Gradients are checked individually (no delay allowed if they start with non-zero amplitude, must be aligned to
    the block boundary if they end with non-zero amplitude) and at the connection points to the previous and, when
    an existing block is overwritten, to the next block. The boundary values of the neighboring blocks are cached in
    `self.grad_check_data_prev` / `self.grad_check_data_next` to avoid repeated lookups when blocks are added
    sequentially. The cached values are stored in the physical (rotated) frame; for blocks with a rotation
    extension (`rot_quaternion`), they are transformed into the logical frame of the current block before the
    comparison and the boundary values of the current block are rotated to the physical frame before caching.

    Parameters
    ----------
    block_index : int
        Index of the block being set.
    check_g : dict
        Maps the gradient axis index (0: x, 1: y, 2: z) to a SimpleNamespace with `start` and `stop` tuples
        (time, amplitude) of the gradient on that axis; (0, 0) if there is no gradient on that axis.
    duration : float
        Duration of the block being set.
    rot_quaternion : array-like, optional
        Unit quaternion [w, x, y, z] of the block's rotation extension, or None if the block is not rotated.

    Raises
    ------
    RuntimeError
        If a gradient starting with non-zero amplitude is in the first block or has a delay, if a gradient ending
        with non-zero amplitude is not aligned to the block boundary, or if the amplitude at the connection point to
        the previous/next block differs by more than the maximum slew rate allows within one gradient raster.
    """
    # Maximum amplitude change that is achievable within one gradient raster; used as tolerance for all checks
    max_amp_step = self.system.max_slew * self.system.grad_raster_time
    is_first_block = block_index == 1 or self.next_free_block_ID == 1

    # Rotation matrix of this block; transforms gradients from the logical frame to the physical (scanner) frame
    if rot_quaternion is None:
        rot_mat = np.eye(3)
    else:
        q = np.asarray(rot_quaternion, dtype=float)
        rot_mat = Rotation.from_quat(q, scalar_first=True).as_matrix()

    # Checks on the gradients of this block alone
    for ax in range(3):
        (start_time, start_amp), (stop_time, stop_amp) = check_g[ax].start, check_g[ax].stop
        if abs(start_amp) > max_amp_step:
            if is_first_block:
                raise RuntimeError('First gradient in the first block has to start at 0.')
            if start_time > eps:
                raise RuntimeError('No delay allowed for gradients which start with a non-zero amplitude')
        if abs(stop_amp) > max_amp_step and abs(stop_time - duration) > 1e-7:
            raise RuntimeError('A gradient that does not end at zero needs to be aligned to the block boundary.')

    # Checks at the connection points to the neighboring blocks
    has_prev_block = self.next_free_block_ID > 1
    has_next_block = _update_neighbor_grad_cache(self, block_index)

    # The cached neighbor values are stored in the physical frame; rotate them into the logical frame of this block
    prev_logical = rot_mat.T @ np.asarray(self.grad_check_data_prev.last_grad_vals, dtype=float)
    next_logical = rot_mat.T @ np.asarray(self.grad_check_data_next.first_grad_vals, dtype=float)

    stop_amps = np.array([check_g[ax].stop[1] for ax in range(3)])
    start_amps = np.array([check_g[ax].start[1] for ax in range(3)])

    for ax in range(3):
        if has_prev_block and abs(prev_logical[ax] - start_amps[ax]) > max_amp_step:
            raise RuntimeError('Two consecutive gradients need to have the same amplitude at the connection point')
        if has_next_block and abs(next_logical[ax] - stop_amps[ax]) > max_amp_step:
            raise RuntimeError('Two consecutive gradients need to have the same amplitude at the connection point')

    # Cache the boundary values of this block in the physical frame: they are the "previous" values for
    # block_index + 1 and the "next" values for block_index - 1
    self.grad_check_data_prev.last_grad_vals = rot_mat @ stop_amps
    self.grad_check_data_next.first_grad_vals = rot_mat @ start_amps

    self.grad_check_data_prev.valid_for_block_num = block_index
    self.grad_check_data_next.valid_for_block_num = block_index - 1


def _update_neighbor_grad_cache(self: Sequence, block_index: int) -> bool:
    """
    Make sure `self.grad_check_data_prev` holds the last gradient values of the block preceding `block_index` and
    `self.grad_check_data_next` the first gradient values of the block following it. The gradient library is only
    queried if the cached values do not belong to these blocks.

    Returns
    -------
    has_next_block : bool
        True if a block follows `block_index` (only possible when an existing block is overwritten), i.e. if the
        connection to the next block has to be checked.
    """
    block_ids = list(self.block_events)
    is_existing_block = block_index in self.block_events

    # Previous block: only relevant if the sequence is not empty
    if self.next_free_block_ID > 1 and self.grad_check_data_prev.valid_for_block_num != block_index - 1:
        if is_existing_block:
            idx = block_ids.index(block_index)
            prev_block_index = block_ids[idx - 1] if idx > 0 else None
        else:
            prev_block_index = block_ids[-1]  # New block is appended (also with non-contiguous numbering)
        self.grad_check_data_prev.last_grad_vals = _boundary_grad_vals(self, prev_block_index, last=True)

    # Next block: only exists if an existing block is overwritten
    has_next_block = False
    if block_index < self.next_free_block_ID and self.grad_check_data_next.valid_for_block_num != block_index:
        next_block_index = None
        if is_existing_block:
            idx = block_ids.index(block_index)
            next_block_index = block_ids[idx + 1] if idx < len(block_ids) - 1 else None
        has_next_block = next_block_index is not None
        self.grad_check_data_next.first_grad_vals = _boundary_grad_vals(self, next_block_index, last=False)
    elif is_existing_block and self.grad_check_data_next.valid_for_block_num == block_index:
        has_next_block = True  # Cached values of the next block are still valid

    return has_next_block


def _boundary_grad_vals(self: Sequence, block_id: int | None, last: bool) -> np.ndarray:
    """
    First (`last=False`) or last (`last=True`) gradient amplitude per axis of block `block_id`, in the physical
    (rotated) frame. Zero if the block does not exist or has no arbitrary gradient on that axis (trapezoids always
    start and end at zero).
    """
    vals = np.zeros(3)
    if block_id is not None:
        for ax in range(3):
            grad_id = self.block_events[block_id][ax + 2]
            if grad_id != 0 and self.grad_library.type[grad_id] == 'g':
                vals[ax] = self.grad_library.data[grad_id][2 if last else 1]  # data = (amplitude, first, last, ...)

        # Apply the rotation of that block (if any) to obtain the physical-frame boundary values
        rot_quaternion = _block_rotation_quaternion(self, block_id)
        if rot_quaternion is not None:
            vals = Rotation.from_quat(rot_quaternion, scalar_first=True).as_matrix() @ vals
    return vals


def _block_rotation_quaternion(self: Sequence, block_id: int) -> np.ndarray | None:
    """
    Rotation quaternion [w, x, y, z] of the ROTATIONS extension attached to block `block_id`, or None if the block
    has no rotation extension.
    """
    if 'ROTATIONS' not in self.extension_string_idx:
        return None

    rot_ext_type = self.get_extension_type_ID('ROTATIONS', update=False)
    next_ext_id = self.block_events[block_id][6]
    while next_ext_id != 0:
        ext_data = self.extensions_library.data[next_ext_id]
        if ext_data[0] == rot_ext_type:
            return np.asarray(self.rotation_library.data[ext_data[1]], dtype=float)
        next_ext_id = ext_data[2]

    return None
