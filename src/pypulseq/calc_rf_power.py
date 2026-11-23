from types import SimpleNamespace

import numpy as np


def calc_rf_power(
    rf: SimpleNamespace,
    dt: float | None = 1e-6,
) -> tuple[float, float, float]:
    """
    Calculate the relative power of the RF pulse.

    Returns the (relative) energy of the pulse expressed in the units of RF amplitude squared multiplied by time,
    e.g. in Pulseq these are Hz * Hz * s = Hz. Sounds strange, but is true.

    Parameters
    ----------
    rf : SimpleNamespace
        Pulseq RF event.
    dt : float or None, default=1e-6
        Resampling step in s. The signal is resampled onto a uniform grid of `round(rf.shape_dur / dt)` bins centered
        at `(arange(nn) + 0.5) * dt`, with zero-padding outside `rf.t` (same as `mr.calcRfPower` in MATLAB Pulseq).
        If None, the signal is integrated directly on its native `rf.t` samples using trapezoid integration instead
        of resampling.

    Returns
    -------
    total_energy : float
        Relative energy of the RF pulse in Hz (= Hz**2 * s).
    peak_pwr : float
        Peak power of the RF pulse in Hz**2.
    rf_rms : float
        RMS B1 amplitude of the RF pulse in Hz.

    Notes
    -----
    The power and RF amplitude calculated by this function are relative as they are calculated in units of Hz**2 or
    Hz. The RF amplitude can be converted to T by dividing the resulting value by gamma. Correspondingly, the power
    can be converted to mT**2 * s by dividing the given value by gamma**2. Nonetheless, the absolute SAR is related to
    the electric field, so the further scaling coefficient is both tx-coil-dependent (e.g. depends on the coil design)
    and also subject-dependent (e.g. depends on the reference voltage).
    """
    if dt is None:
        # Integrate the native samples directly (manual trapezoid rule,
        # equivalent to np.trapezoid which requires numpy >= 2.0)
        rfs_sq = np.abs(rf.signal) ** 2
        total_energy = np.sum(np.diff(rf.t) * (rfs_sq[:-1] + rfs_sq[1:]) / 2)
        peak_pwr = np.max(rfs_sq)
    else:
        # Resample the pulse onto a uniform time grid
        nn = round(rf.shape_dur / dt)
        t = (np.arange(nn) + 0.5) * dt
        rfs = np.interp(t, rf.t, rf.signal, left=0, right=0)

        rfs_sq = np.abs(rfs) ** 2
        total_energy = np.sum(rfs_sq) * dt
        peak_pwr = np.max(rfs_sq)

    rf_rms = np.sqrt(total_energy / rf.shape_dur)
    return total_energy, peak_pwr, rf_rms
