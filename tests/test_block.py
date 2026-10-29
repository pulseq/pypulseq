import numpy as np
import pypulseq as pp
import pytest

# Gradient definitions used in tests
gx_trap = pp.make_trapezoid('x', area=1000, duration=1e-3)
gx_extended = pp.make_extended_trapezoid('x', amplitudes=[0, 100000, 0], times=[0, 1e-4, 2e-4])
gx_extended_delay = pp.make_extended_trapezoid('x', amplitudes=[0, 100000, 0], times=[1e-4, 2e-4, 3e-4])
gx_endshigh = pp.make_extended_trapezoid('x', amplitudes=[0, 100000, 100000], times=[0, 1e-4, 2e-4])
gx_startshigh = pp.make_extended_trapezoid('x', amplitudes=[100000, 100000, 0], times=[0, 1e-4, 2e-4])
gx_startshigh2 = pp.make_extended_trapezoid('x', amplitudes=[200000, 100000, 0], times=[0, 1e-4, 2e-4])
gx_allhigh = pp.make_extended_trapezoid('x', amplitudes=[100000, 100000, 100000], times=[0, 1e-4, 2e-4])
gx_oversampled = pp.make_arbitrary_grad(
    'x', 10000 * np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1]), first=0, last=0, oversampling=True
)
delay = pp.make_delay(1e-3)


def test_add_block_none_only_noop():
    seq = pp.Sequence()
    seq.add_block(None)
    assert len(seq.block_events) == 0


def test_add_block_ignores_none():
    seq1 = pp.Sequence()
    seq1.add_block(pp.make_delay(1))

    seq2 = pp.Sequence()
    seq2.add_block(None, pp.make_delay(1), None)

    assert list(seq1.block_durations.values()) == list(seq2.block_durations.values())


def test_gradient_continuity1():
    # Trap followed by extended gradient: No error
    seq = pp.Sequence()
    seq.add_block(gx_trap)
    seq.add_block(gx_extended)
    seq.add_block(gx_trap)


def test_gradient_continuity2():
    # Trap followed by non-zero gradient
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_trap)
        seq.add_block(gx_startshigh)  # raises


def test_gradient_continuity3():
    # Gradient starts at non-zero in first block
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_startshigh)  # raises


def test_gradient_continuity4():
    # Gradient starts and ends at non-zero
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(delay)
        seq.add_block(gx_allhigh)


def test_gradient_continuity5():
    # Gradient starts at zero and has a delay: No error
    seq = pp.Sequence()
    seq.add_block(gx_extended_delay)


def test_gradient_continuity6():
    # Gradient starts at non-zero in other blocks
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(delay)
        seq.add_block(gx_startshigh)  # raises


def test_gradient_continuity7():
    # Gradient ends high and is followed by empty block
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(delay)  # raises


def test_gradient_continuity8():
    # Gradient ends high and is followed by trapezoid
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(gx_trap)  # raises


def test_gradient_continuity9():
    # Gradient ends high and is followed by connecting gradient: No error
    seq = pp.Sequence()
    seq.add_block(gx_endshigh)
    seq.add_block(gx_startshigh)


def test_gradient_continuity10():
    # Gradient in last block ends high: No error, this is caught by seq.write()
    seq = pp.Sequence()
    seq.add_block(gx_endshigh)


def test_gradient_continuity11():
    # Non-zero, but non-connecting gradients
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(gx_startshigh2)


## Test gradient continuity checks in set_block


def test_gradient_continuity_setblock1():
    # Use set_block to insert gradient
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(delay)
        seq.add_block(delay)
        seq.add_block(delay)

        seq.set_block(1, gx_startshigh)


def test_gradient_continuity_setblock2():
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(delay)
        seq.add_block(delay)
        seq.add_block(delay)

        seq.set_block(2, gx_startshigh)


def test_gradient_continuity_setblock3():
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(delay)
        seq.add_block(delay)
        seq.add_block(delay)

        seq.set_block(3, gx_startshigh)


def test_gradient_continuity_setblock4():
    # Overwrite valid gradient with empty block
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(gx_allhigh)
        seq.add_block(gx_startshigh)

        seq.set_block(2, delay)


def test_gradient_continuity_setblock5():
    # Overwrite valid gradient with gradient that is valid on one side
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(gx_allhigh)
        seq.add_block(gx_startshigh)

        seq.set_block(2, gx_startshigh)


def test_gradient_continuity_setblock6():
    # Add new gradient with non-contiguous block numbering
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(gx_allhigh)
        seq.add_block(gx_startshigh)

        seq.set_block(6, gx_startshigh)


def test_gradient_continuity_setblock7():
    # Valid sequence with non-contiguous block numbering
    seq = pp.Sequence()
    seq.set_block(10, gx_endshigh)
    seq.set_block(5, gx_allhigh)
    seq.set_block(7, gx_startshigh)

    assert list(seq.block_events.keys()) == [10, 5, 7]


def test_gradient_continuity_setblock8():
    # Overwrite the last block with a gradient that ends high: No error expected (same as add_block, see continuity10)
    seq = pp.Sequence()
    seq.add_block(gx_trap)
    seq.add_block(gx_trap)
    seq.set_block(2, gx_endshigh)


def test_gradient_continuity_setblock9():
    # Overwrite a middle block with a gradient that connects on both sides: No error expected
    seq = pp.Sequence()
    seq.add_block(gx_endshigh)
    seq.add_block(gx_allhigh)
    seq.add_block(gx_startshigh)
    seq.set_block(2, gx_allhigh)

    # Continue appending afterwards: block 3 ends at zero, so a trapezoid is fine ...
    seq.add_block(gx_trap)

    # Overwriting block 3 with a gradient that ends high is not allowed, because block 4 starts at zero
    with pytest.raises(RuntimeError):
        seq.set_block(3, gx_allhigh)


def test_gradient_continuity_zero_duration_block():
    # Gradient ends high and is followed by a block of zero duration (label only): must raise like continuity7
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.add_block(gx_endshigh)
        seq.add_block(pp.make_label(label='LIN', type='SET', value=1))  # raises


def test_gradient_continuity_first_block_noncontiguous():
    # First block of an empty sequence starts high, even with non-contiguous numbering: must raise
    with pytest.raises(RuntimeError):
        seq = pp.Sequence()
        seq.set_block(10, gx_startshigh)


## Test gradient continuity checks with the rotation extension

# Rotation extension events used in the continuity tests
rot_none = None  # no rotation extension
rot_eye = pp.make_rotation(phi=0.0)  # identity rotation
rot_z90 = pp.make_rotation(phi=np.pi / 2)  # 90 deg about z: maps a logical x gradient onto physical y

gy_startshigh = pp.make_extended_trapezoid('y', amplitudes=[100000, 100000, 0], times=[0, 1e-4, 2e-4])


def _events(*args):
    """Filter out None events (set_block does not ignore them, unlike add_block)."""
    return tuple(a for a in args if a is not None)


def _add_blocks(seq, blocks, rot):
    """Add each block (tuple of events) to seq with the rotation extension `rot` attached."""
    for events in blocks:
        seq.add_block(*events, rot)


def _set_block(seq, block_index, events, rot):
    """Replace block with events and the rotation extension `rot` attached."""
    seq.set_block(block_index, *_events(*events, rot))


# A uniform rotation (same rotation on every block) must never change the outcome of the
# gradient continuity checks
@pytest.mark.parametrize('rot', [rot_none, rot_eye, rot_z90], ids=['none', 'eye', 'z90'])
def test_gradient_continuity_uniform_rotation_ok(rot):
    # Trap followed by extended gradient
    seq = pp.Sequence()
    _add_blocks(seq, [(gx_trap,), (gx_extended,), (gx_trap,)], rot)

    # Gradient starts at zero and has a delay
    seq = pp.Sequence()
    seq.add_block(gx_extended_delay, rot)

    # Gradient ends high and is followed by connecting gradient
    seq = pp.Sequence()
    _add_blocks(seq, [(gx_endshigh,), (gx_startshigh,)], rot)

    # Gradient in last block ends high: No error, this is caught by seq.write()
    seq = pp.Sequence()
    seq.add_block(gx_endshigh, rot)

    # Valid sequence with non-contiguous block numbering
    seq = pp.Sequence()
    _set_block(seq, 10, (gx_endshigh,), rot)
    _set_block(seq, 5, (gx_allhigh,), rot)
    _set_block(seq, 7, (gx_startshigh,), rot)

    assert list(seq.block_events.keys()) == [10, 5, 7]


# Lists of blocks (tuples of events) where the last add_block must raise RuntimeError
uniform_rotation_raises = [
    [(gx_trap,), (gx_startshigh,)],  # trap followed by non-zero gradient
    [(gx_startshigh,)],  # gradient starts at non-zero in first block
    [(delay,), (gx_allhigh,)],  # gradient starts and ends at non-zero
    [(delay,), (gx_startshigh,)],  # gradient starts at non-zero in other blocks
    [(gx_endshigh,), (delay,)],  # gradient ends high followed by empty block
    [(gx_endshigh,), (gx_trap,)],  # gradient ends high followed by trapezoid
    [(gx_endshigh,), (gx_startshigh2,)],  # non-zero, but non-connecting gradients
]


@pytest.mark.parametrize('blocks', uniform_rotation_raises)
@pytest.mark.parametrize('rot', [rot_none, rot_eye, rot_z90], ids=['none', 'eye', 'z90'])
def test_gradient_continuity_uniform_rotation_raises(rot, blocks):
    seq = pp.Sequence()
    _add_blocks(seq, blocks[:-1], rot)
    with pytest.raises(RuntimeError):
        seq.add_block(*_events(*blocks[-1], rot))


@pytest.mark.parametrize('block_index', [1, 2, 3])
@pytest.mark.parametrize('rot', [rot_none, rot_eye, rot_z90], ids=['none', 'eye', 'z90'])
def test_gradient_continuity_setblock_uniform_rotation_raises(rot, block_index):
    # Use set_block to insert a gradient that starts at non-zero between delay blocks
    seq = pp.Sequence()
    _add_blocks(seq, [(delay,), (delay,), (delay,)], rot)

    with pytest.raises(RuntimeError):
        _set_block(seq, block_index, (gx_startshigh,), rot)


@pytest.mark.parametrize(
    'block_index,events',
    [(2, (delay,)), (2, (gx_startshigh,)), (6, (gx_startshigh,))],
    ids=['overwrite-mid-with-delay', 'overwrite-mid-startshigh', 'noncontiguous-startshigh'],
)
@pytest.mark.parametrize('rot', [rot_none, rot_eye, rot_z90], ids=['none', 'eye', 'z90'])
def test_gradient_continuity_setblock_uniform_rotation_raises2(rot, block_index, events):
    # Overwrite a valid gradient with a block that breaks the continuity
    seq = pp.Sequence()
    _add_blocks(seq, [(gx_endshigh,), (gx_allhigh,), (gx_startshigh,)], rot)

    with pytest.raises(RuntimeError):
        _set_block(seq, block_index, events, rot)


# A rotation that differs between connected blocks breaks the connection; the same rotation
# keeps it
def test_gradient_continuity_mixed_rotations():
    # Same non-identity rotation on both blocks keeps the connection
    seq = pp.Sequence()
    _add_blocks(seq, [(gx_endshigh,), (gx_startshigh,)], rot_z90)

    # Rotation only on the second block breaks the connection
    seq = pp.Sequence()
    seq.add_block(gx_endshigh, rot_eye)
    with pytest.raises(RuntimeError):
        seq.add_block(gx_startshigh, rot_z90)

    # Rotation only on the first block breaks the connection
    seq = pp.Sequence()
    seq.add_block(gx_endshigh, rot_z90)
    with pytest.raises(RuntimeError):
        seq.add_block(gx_startshigh, rot_eye)


def test_gradient_continuity_rotation_physical_frame():
    # Rz(90) maps a logical x gradient onto physical y: connecting on (logical) y in the next
    # block passes because the boundary values match in the physical frame
    seq = pp.Sequence()
    seq.add_block(gx_endshigh, rot_z90)
    seq.add_block(gy_startshigh)

    # Rz(90) maps a logical y gradient onto physical -x, so it does not connect to an
    # unrotated x gradient ending high
    seq = pp.Sequence()
    seq.add_block(gx_endshigh)
    with pytest.raises(RuntimeError):
        seq.add_block(gy_startshigh, rot_z90)


@pytest.mark.parametrize('block_index', [1, 2, 3])
@pytest.mark.parametrize('base_rot', [rot_none, rot_eye], ids=['none', 'eye'])
def test_gradient_continuity_setblock_mixed_rotation_raises(base_rot, block_index):
    # Overwriting a block with a differently rotated gradient breaks the connection on both sides
    seq = pp.Sequence()
    _add_blocks(seq, [(gx_endshigh,), (gx_allhigh,), (gx_startshigh,)], base_rot)

    grad = [gx_endshigh, gx_allhigh, gx_startshigh][block_index - 1]
    with pytest.raises(RuntimeError):
        seq.set_block(block_index, grad, rot_z90)


def test_gradient_continuity_setblock10():
    # Add new gradient after valid override
    seq = pp.Sequence()
    seq.add_block(gx_endshigh)
    seq.add_block(gx_allhigh)
    seq.add_block(gx_startshigh)
    seq.set_block(2, gx_allhigh)
    seq.add_block(gx_endshigh)
    seq.add_block(gx_allhigh)
    seq.add_block(gx_startshigh)


def test_get_block_oversampled_grad():
    # get_block returns the oversampled gradient that was added, with the same timing
    seq = pp.Sequence()
    seq.add_block(gx_oversampled)
    gx = seq.get_block(1).gx

    assert gx.shape_dur == pytest.approx(gx_oversampled.shape_dur)
    assert np.allclose(gx.tt, gx_oversampled.tt)
    assert gx.area == pytest.approx(gx_oversampled.area)

    ok, errors = seq.check_timing()
    assert ok, errors


# TODO: Add other block functionality tests
