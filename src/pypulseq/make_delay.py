from types import SimpleNamespace

import numpy as np


def make_delay(d: float) -> SimpleNamespace:
    """
    Creates a delay event.

    Parameters
    ----------
    d : float
        Delay time in seconds (s).

    Returns
    -------
    delay : SimpleNamespace
        Delay event.

    Raises
    ------
    ValueError
        If delay is invalid (not finite or < 0).
    """
    delay = SimpleNamespace()
    if not np.isfinite(d) or d < 0:
        raise ValueError(f'Delay {d * 1e3:.2f} ms is invalid')
    delay.type = 'delay'
    delay.delay = d
    return delay
