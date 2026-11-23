import pypulseq as pp
import pytest
from pypulseq.supported_labels_rf_use import _SUPPORTED_LABELS, add_supported_label, get_supported_labels


@pytest.fixture(autouse=True)
def _restore_supported_labels():
    """Custom labels are registered in module state; restore it after each test."""
    original = list(_SUPPORTED_LABELS)
    yield
    _SUPPORTED_LABELS[:] = original


def test_add_valid_label():
    labels = add_supported_label('MYCUSTOM')
    assert 'MYCUSTOM' in labels
    assert 'MYCUSTOM' in get_supported_labels()


def test_add_duplicate_label_warns():
    add_supported_label('TESTDUP')
    with pytest.warns(UserWarning, match='already known'):
        add_supported_label('TESTDUP')
    # the label must not be appended twice (label IDs are positional)
    assert get_supported_labels().count('TESTDUP') == 1


@pytest.mark.parametrize('bad_label', [123, None, '', 'HAS SPACE', 'HAS\tTAB'])
def test_add_invalid_label(bad_label):
    with pytest.raises(ValueError):
        add_supported_label(bad_label)


def test_add_multiple_labels():
    add_supported_label('CUSTOM1')
    add_supported_label('CUSTOM2')
    labels = get_supported_labels()
    assert 'CUSTOM1' in labels
    assert 'CUSTOM2' in labels


def test_custom_label_with_make_label():
    add_supported_label('MYTEST')

    label = pp.make_label('MYTEST', 'SET', 42)
    assert label.type == 'labelset'
    assert label.label == 'MYTEST'
    assert label.value == 42

    # custom labels are counters, so 'INC' is allowed
    label = pp.make_label('MYTEST', 'INC', -3)
    assert label.type == 'labelinc'
    assert label.value == -3


def test_flag_still_rejects_inc():
    """Flags must keep rejecting 'INC' after custom labels were appended."""
    add_supported_label('MYCOUNTER')
    with pytest.raises(ValueError, match='not compatible with flags'):
        pp.make_label('NAV', 'INC', 1)
    # 'TRID' is a counter (last built-in label): INC must stay allowed
    assert pp.make_label('TRID', 'INC', 1).type == 'labelinc'


def test_unregistered_label_rejected():
    with pytest.raises(ValueError, match='Invalid label'):
        pp.make_label('NOTREGISTERED', 'SET', 1)


def test_custom_label_roundtrip(tmp_path):
    """A sequence with custom labels must survive write/read, with the label names in the .seq file."""
    add_supported_label('MYSET')
    add_supported_label('MYINC')

    seq = pp.Sequence()
    seq.add_block(
        pp.make_label('MYSET', 'SET', 7),
        pp.make_label('MYINC', 'INC', 1),
        pp.make_delay(1e-3),
    )

    path = tmp_path / 'custom_label.seq'
    seq.write(path)

    content = path.read_text()
    assert 'MYSET' in content
    assert 'MYINC' in content

    seq2 = pp.Sequence()
    seq2.read(path)
    labels = seq2.get_block(1).label
    assert labels[0].label == 'MYSET'
    assert labels[0].value == 7
    assert labels[1].label == 'MYINC'
    assert labels[1].value == 1


def test_read_unregistered_label_error(tmp_path):
    """Reading a file with an unregistered label must fail with a helpful message."""
    add_supported_label('GHOST')

    seq = pp.Sequence()
    seq.add_block(pp.make_label('GHOST', 'SET', 1), pp.make_delay(1e-3))
    path = tmp_path / 'ghost.seq'
    seq.write(path)

    # simulate a fresh session where the label was never registered
    _SUPPORTED_LABELS.remove('GHOST')

    seq2 = pp.Sequence()
    with pytest.raises(ValueError, match='add_supported_label'):
        seq2.read(path)
