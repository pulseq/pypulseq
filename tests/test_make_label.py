"""Tests for the make_label module"""

import numpy as np
import pytest
from pypulseq import make_label


def test_invalid_label():
    with pytest.raises(ValueError, match='Invalid label'):
        make_label(label='FOO', type='SET', value=1)


def test_invalid_type():
    with pytest.raises(ValueError, match="Invalid type\\. Must be one of 'SET' or 'INC'"):
        make_label(label='SLC', type='INVALID', value=1)


def test_inc_not_compatible_with_flags():
    with pytest.raises(ValueError, match='labelinc is not compatible with flags'):
        make_label(label='NAV', type='INC', value=1)


def test_invalid_value():
    for invalid in ['one', '2', None, np.array(2), 2 + 0j]:
        with pytest.raises(ValueError, match='Must supply a valid numerical or logical value'):
            make_label(label='SLC', type='SET', value=invalid)

    with pytest.raises(ValueError, match='Label value must be finite'):
        make_label(label='SLC', type='SET', value=np.nan)

    with pytest.raises(ValueError, match='Label value must be finite'):
        make_label(label='SLC', type='SET', value=np.inf)


def test_valid():
    label = make_label(label='SLC', type='SET', value=1)
    assert label.type == 'labelset'
    assert label.value == 1

    label = make_label(label='ECO', type='INC', value=1)
    assert label.type == 'labelinc'

    # Flags are compatible with SET
    label = make_label(label='NAV', type='SET', value=True)
    assert label.value == 1

    # Values are converted to int
    label = make_label(label='SLC', type='SET', value=2.7)
    assert label.value == 2
    assert isinstance(label.value, int)


@pytest.mark.parametrize(
    'value',
    [np.int64(2), np.int32(2), np.uint8(2), np.float32(2.0), np.float64(2.0), np.bool_(True), True, 2, 2.0],
)
def test_numpy_scalar_and_builtin_values(value):
    # Numpy scalars (and Python bool/int/float) are accepted and converted to int
    label = make_label(label='PAR', type='SET', value=value)
    assert label.value == int(value)
    assert isinstance(label.value, int)
