from typing import List, Tuple
from warnings import warn

# List of supported labels. Order is significant: label events are stored with a 1-based numeric ID equal to the
# index in this list, and .seq files map label IDs to label names by position. Custom labels can be registered with
# add_supported_label(); they are always appended at the end to keep the IDs of existing labels stable.
_SUPPORTED_LABELS: List[str] = [
    'SLC',
    'SEG',
    'REP',
    'AVG',
    'SET',
    'ECO',
    'PHS',
    'LIN',
    'PAR',
    'ACQ',
    'NAV',
    'REV',
    'SMS',
    'REF',
    'IMA',  # For parallel imaging
    'NOISE',  # noise adjust scan, for iPAT acceleration
    'PMC',  # for MoCo/PMC Pulseq version to recognize blocks that can be prospectively corrected for motion
    'NOROT',
    'NOPOS',
    'NOSCL',  # instruct the interpreter to ignore the position, rotation or scaling of the FOV specified on the UI
    'ONCE',  # a 3-state flag that instructs the interpreter to alter the sequence when executing multiple repeats as follows: blocks with ONCE==0 are executed on every repetition; ONCE==1: only the first repetition; ONCE==2: only the last repetition
    'TRID',  # an integer ID of the TR (sequence segment) used by the GE interpreter to optimize the execution on the scanner
]

# Labels that represent flags. Per the Pulseq specification, 'INC' label events are not applicable to them.
LABEL_FLAGS = frozenset(('NAV', 'REV', 'SMS', 'REF', 'IMA', 'NOISE', 'PMC', 'NOROT', 'NOPOS', 'NOSCL', 'ONCE'))


def get_supported_labels() -> Tuple[str, ...]:
    """
    Returns
    -------
    tuple
        Supported labels.
    """
    return tuple(_SUPPORTED_LABELS)


def add_supported_label(new_label: str) -> Tuple[str, ...]:
    """
    Register a custom label, equivalent to `mr.addCustomLabel` in MATLAB Pulseq.

    The label is appended to the list returned by `get_supported_labels()`, which makes it usable in
    `pypulseq.make_label` and when reading .seq files that contain it. Custom labels are not flag labels,
    so both 'SET' and 'INC' label events are allowed for them.

    Parameters
    ----------
    new_label : str
        Name of the label to register. Must be a non-empty string without whitespace.

    Returns
    -------
    tuple
        The updated list of supported labels.
    """
    if not isinstance(new_label, str) or not new_label or any(c.isspace() for c in new_label):
        raise ValueError('New label must be a non-empty string without whitespace.')

    if new_label in _SUPPORTED_LABELS:
        warn(f'Label {new_label} is already known.', stacklevel=2)
        return tuple(_SUPPORTED_LABELS)

    _SUPPORTED_LABELS.append(new_label)
    return tuple(_SUPPORTED_LABELS)


def get_supported_rf_uses() -> Tuple[str, str, str, str, str, str, str]:
    """
    Returns
    -------
    tuple
        Supported RF use labels.
    """
    return 'excitation', 'refocusing', 'inversion', 'saturation', 'preparation', 'other', 'undefined'
