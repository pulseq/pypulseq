from __future__ import annotations

import typing

from pypulseq import eps

if typing.TYPE_CHECKING:
    from pypulseq.Sequence.sequence import Sequence


def ext_grad_check(self: Sequence, block_index: int, check_g: dict, duration: float) -> None:
    """
    Check the gradient continuity of the block that is being set at `block_index`.

    Gradients are checked individually (no delay allowed if they start with non-zero amplitude, must be aligned to
    the block boundary if they end with non-zero amplitude) and at the connection points to the previous and, when
    an existing block is overwritten, to the next block. The boundary values of the neighboring blocks are cached in
    `self.grad_check_data_prev` / `self.grad_check_data_next` to avoid repeated lookups when blocks are added
    sequentially.

    Parameters
    ----------
    block_index : int
        Index of the block being set.
    check_g : dict
        Maps the gradient axis index (0: x, 1: y, 2: z) to a SimpleNamespace with `start` and `stop` tuples
        (time, amplitude) of the gradient on that axis; (0, 0) if there is no gradient on that axis.
    duration : float
        Duration of the block being set.

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

    # TODO (v1.5.1 rotation extension): obtain the rotation of this block to transform between logical/physical axes

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

    # TODO (v1.5.1 rotation extension): rotate the cached neighbor values into the logical frame of this block
    for ax in range(3):
        start_amp, stop_amp = check_g[ax].start[1], check_g[ax].stop[1]
        if has_prev_block and abs(self.grad_check_data_prev.last_grad_vals[ax] - start_amp) > max_amp_step:
            raise RuntimeError('Two consecutive gradients need to have the same amplitude at the connection point')
        if has_next_block and abs(self.grad_check_data_next.first_grad_vals[ax] - stop_amp) > max_amp_step:
            raise RuntimeError('Two consecutive gradients need to have the same amplitude at the connection point')

        # TODO (v1.5.1 rotation extension): rotate back to physical axes before caching
        # Cache the boundary values of this block: they are the "previous" values for block_index + 1 and the
        # "next" values for block_index - 1
        self.grad_check_data_prev.last_grad_vals[ax] = stop_amp
        self.grad_check_data_next.first_grad_vals[ax] = start_amp

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


def _boundary_grad_vals(self: Sequence, block_id: int | None, last: bool) -> list:
    """
    First (`last=False`) or last (`last=True`) gradient amplitude per axis of block `block_id`. Zero if the block
    does not exist or has no arbitrary gradient on that axis (trapezoids always start and end at zero).
    """
    vals = [0.0, 0.0, 0.0]
    if block_id is not None:
        for ax in range(3):
            grad_id = self.block_events[block_id][ax + 2]
            if grad_id != 0 and self.grad_library.type[grad_id] == 'g':
                vals[ax] = self.grad_library.data[grad_id][2 if last else 1]  # data = (amplitude, first, last, ...)
    return vals
