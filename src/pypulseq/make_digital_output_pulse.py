from types import SimpleNamespace
from typing import Union

import numpy as np

from pypulseq.opts import Opts


def make_digital_output_pulse(
    channel: str,
    delay: float = 0.0,
    duration: float = 4e-3,
    system: Union[Opts, None] = None,
) -> SimpleNamespace:
    """
    Create a digital output pulse event. Emits an output signal on a given channel with optional
    given delay and duration. For external trigger synchronization events, see `make_trigger`.

    Parameters
    ----------
    channel : str
        Must be one of 'osc0','osc1', or 'ext1'.
    delay : float, default=0
        Delay in seconds (s).
    duration : float, default=4e-3
        Duration of trigger event in seconds (s).
    system : Opts, default=Opts()
        System limits.

    Returns
    -------
    trig : SimpleNamespace
        Trigger event.

    Raises
    ------
    ValueError
        If `channel` is invalid. Must be one of 'osc0','osc1', or 'ext1'.
        If `duration` or `delay` is negative or not finite.
    """
    if system is None:
        system = Opts.default

    if channel not in ['osc0', 'osc1', 'ext1']:
        raise ValueError(f"Channel {channel} is invalid. Must be one of 'osc0','osc1', or 'ext1'.")

    if not np.isfinite(duration) or duration < 0:
        raise ValueError('Duration must be non-negative.')

    if not np.isfinite(delay) or delay < 0:
        raise ValueError('Delay must be non-negative.')

    trig = SimpleNamespace()
    trig.type = 'output'
    trig.channel = channel
    trig.delay = delay
    trig.duration = duration
    if trig.duration <= system.grad_raster_time:
        trig.duration = system.grad_raster_time

    return trig
