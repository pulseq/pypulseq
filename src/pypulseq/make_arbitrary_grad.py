from types import SimpleNamespace

import numpy as np

from pypulseq import eps
from pypulseq.opts import Opts
from pypulseq.utils.tracing import trace, trace_enabled


def make_arbitrary_grad(
    channel: str,
    waveform: np.ndarray,
    first: float | None = None,
    last: float | None = None,
    delay: float = 0.0,
    max_grad: float | None = None,
    max_slew: float | None = None,
    system: Opts | None = None,
    oversampling: bool = False,
) -> SimpleNamespace:
    """
    Creates a gradient event from an arbitrary waveform.

    Note that the sample points are assumed to be equally spaced by `system.grad_raster_time`
    and that the given waveform values are the values in the middle of each raster interval.

    The duration of the gradient is thus given by the number of samples times `system.grad_raster_time`.

    See also `pypulseq.Sequence.sequence.Sequence.add_block()`.

    Parameters
    ----------
    channel : str
        Orientation of gradient event of arbitrary shape. Must be one of `x`, `y` or `z`.
    waveform : numpy.ndarray
        Arbitrary waveform. Must be a non-empty 1-D array.
    first : float, default=None
        Gradient value at the start of the gradient event. (t=0)
        Will default to a linear extrapolated value if not provided.
    last : float, default=None
        Gradient value at the end of the gradient event. (t=duration)
        Will default to a linear extrapolated value if not provided.
    delay : float, default=0
        Delay in seconds (s).
    max_grad : float, default=None
        Maximum gradient strength (Hz/m).
        Will default to `system.max_grad` if not provided or 0.
    max_slew : float, default=None
        Maximum slew rate (Hz/m/s).
        Will default to `system.max_slew` if not provided or 0.
    system : Opts, default=Opts()
        System limits.
        Will default to `pypulseq.opts.default` if not provided.
    oversampling : bool, default=False
        Boolean flag to indicate if gradient is oversampled by a factor of 2.

    Returns
    -------
    grad : SimpleNamespace
        Gradient event with arbitrary waveform.

    Raises
    ------
    ValueError
        If invalid `channel` is passed. Must be one of x, y or z.
        If `waveform` is not a non-empty 1-D array.
        If `waveform` has less than 2 samples and `first` or `last` are not provided.
        If `waveform` has an even number of samples and `oversampling` is enabled.
        If `delay` is negative or not finite.
        If slew rate is violated.
        If gradient amplitude is violated.
    """
    if system is None:
        system = Opts.default

    if channel not in ['x', 'y', 'z']:
        raise ValueError(f'Invalid channel. Must be one of x, y or z. Passed: {channel}')

    waveform = np.atleast_1d(np.squeeze(np.asarray(waveform)))
    if waveform.ndim != 1:
        raise ValueError(f'waveform should have ndim=1. Passed ndim={waveform.ndim}')
    if waveform.size == 0:
        raise ValueError('waveform must not be empty.')
    if len(waveform) < 2 and (first is None or last is None):
        raise ValueError('waveform must contain at least 2 samples unless `first` and `last` are provided.')

    if not np.isfinite(delay) or delay < 0:
        raise ValueError('Delay must be non-negative.')

    if max_grad is None or max_grad == 0:
        max_grad = system.max_grad

    if max_slew is None or max_slew == 0:
        max_slew = system.max_slew

    if first is None or last is None:

        def extrap(a, b):
            # extrapolate by 1 gradient raster (oversampling)
            # or by 1/2 gradient of the raster (non-oversampling)
            return 2 * a - b if oversampling else 0.5 * (3 * a - b)

        if first is None:
            first = extrap(waveform[0], waveform[1])
        if last is None:
            last = extrap(waveform[-1], waveform[-2])

    # Slew rate calculation
    if oversampling:
        edge_scale = system.grad_raster_time / 2
        pre = first - waveform[0]
        post = last - waveform[-1]
    else:
        edge_scale = system.grad_raster_time
        pre = 2 * (first - waveform[0])
        post = 2 * (waveform[-1] - last)

    slew_rate = np.concatenate([[pre], np.diff(waveform), [post]]) / edge_scale

    if max(abs(slew_rate)) > max_slew * (1 + eps):
        raise ValueError(f'Slew rate violation {max(abs(slew_rate)) / max_slew * 100}')
    if max(abs(waveform)) > max_grad + eps:
        raise ValueError(f'Gradient amplitude violation {max(abs(waveform)) / max_grad * 100}')

    grad = SimpleNamespace()
    grad.type = 'grad'
    grad.channel = channel
    grad.waveform = waveform
    grad.delay = delay
    if oversampling:
        if len(waveform) % 2 == 0:
            raise ValueError('When oversampling is active, waveform must have an odd number of samples')
        grad.area = (waveform[::2] * system.grad_raster_time).sum()
        grad.tt = np.arange(1, len(waveform) + 1) * 0.5 * system.grad_raster_time
        grad.shape_dur = (len(waveform) + 1) * 0.5 * system.grad_raster_time
    else:
        grad.area = (waveform * system.grad_raster_time).sum()
        grad.tt = (np.arange(len(waveform)) + 0.5) * system.grad_raster_time
        grad.shape_dur = len(waveform) * system.grad_raster_time
    grad.first = first
    grad.last = last

    if trace_enabled():
        grad.trace = trace()

    return grad
