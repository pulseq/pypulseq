# inserted for trigger support by mveldmann

from types import SimpleNamespace

import numpy as np

from pypulseq.opts import Opts


def make_trigger(
    channel: str,
    delay: float = 0.0,
    duration: float = 0.0,
    system: Opts | None = None,
) -> SimpleNamespace:
    """
    Create a trigger halt event for a synchronization with an external signal from a given channel with an optional
    given delay prior to the sync and duration after the sync. Possible channel values: 'physio1','physio2'
    (Siemens specific).

    See also `pypulseq.Sequence.sequence.Sequence.add_block()`.

    Parameters
    ----------
    channel : str
        Must be one of 'physio1' or 'physio2'.
    delay : float, default=0
        Delay in seconds (s).
    duration : float, default=0
        Duration in seconds (s).
    system : Opts, default=Opts()
        System limits.

    Returns
    -------
    trigger : SimpleNamespace
        Trigger event.

    Raises
    ------
    ValueError
        If invalid `channel` is passed. Must be one of 'physio1' or 'physio2'.
        If `duration` or `delay` is negative or not finite.
    """
    if system is None:
        system = Opts.default

    if channel not in ['physio1', 'physio2']:
        raise ValueError(f"Channel {channel} is invalid. Must be one of 'physio1' or 'physio2'.")

    if not np.isfinite(duration) or duration < 0:
        raise ValueError('Duration must be non-negative.')

    if not np.isfinite(delay) or delay < 0:
        raise ValueError('Delay must be non-negative.')

    trigger = SimpleNamespace()
    trigger.type = 'trigger'
    trigger.channel = channel
    trigger.delay = delay
    trigger.duration = duration
    if trigger.duration <= system.grad_raster_time:
        trigger.duration = system.grad_raster_time

    return trigger
