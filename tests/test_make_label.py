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
    with pytest.raises(ValueError, match='Must supply a valid numerical or logical value'):
        make_label(label='SLC', type='SET', value='one')

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
