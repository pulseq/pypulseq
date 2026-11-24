from __future__ import annotations

import typing

import numpy as np
from scipy.signal import convolve, windows

try:
    import sounddevice as sd

    _HAS_SOUNDDEVICE = True
except ImportError:
    _HAS_SOUNDDEVICE = False

if typing.TYPE_CHECKING:
    from pypulseq.Sequence.sequence import Sequence


def seq_sound(
    seq: Sequence,
    time_range: tuple[float, float] = (0.0, np.inf),
    channel_weights: tuple[float, float, float] = (1.0, 1.0, 1.0),
    sample_rate: int = 44100,
    only_produce_sound_data: bool = False,
) -> np.ndarray:
    """
    Play out the sequence through the system speaker and return the sound data.

    The stereo sound data is produced from the gradient waveforms: the X gradient is mapped to the left channel,
    the Y gradient to the right channel and the Z gradient is split equally between both channels. The result
    is Gauss-filtered to reduce high-frequency ringing and normalized to a peak amplitude of 0.95.

    Parameters
    ----------
    seq : Sequence
        The Pulseq sequence to be played.
    time_range : tuple[float, float], default=(0.0, np.inf)
        Start and end time in seconds of the part of the sequence to be played.
        Default is 0 to infinity (entire sequence).
    channel_weights : tuple[float, float, float], default=(1.0, 1.0, 1.0)
        Weights applied to the X, Y and Z gradient waveforms before they are mixed into the sound data.
    sample_rate : int, default=44100
        Sampling rate of the audio waveform in Hz.
    only_produce_sound_data : bool, default=False
        If True, skip the playback and only return the sound data.

    Returns
    -------
    sound_data : np.ndarray
        Sound data of shape (2, N) with rows (left, right) sampled at `sample_rate`, normalized to [-0.95, 0.95].

    Raises
    ------
    RuntimeError
        If `only_produce_sound_data` is False and the optional `sounddevice` package is not installed.
    """
    # --- get waveform data ---
    gw_data = seq.waveforms(append_RF=False, time_range=time_range)
    total_duration = sum(seq.block_durations.values())

    t_start = max(0.0, time_range[0])
    t_end = min(total_duration, time_range[1])

    dwell_time = 1.0 / sample_rate
    sound_length = int(np.floor((t_end - t_start) / dwell_time)) + 1

    sound_data = np.zeros((2, sound_length), dtype=float)

    t_out = t_start + np.arange(sound_length) * dwell_time

    def _interp(channel: np.ndarray, weight: float) -> np.ndarray:
        if channel.shape[1] == 0:
            return np.zeros(sound_length)
        return np.interp(t_out, channel[0], channel[1] * weight, left=0.0, right=0.0)

    # --- X channel → channel 1, Y channel → channel 2, Z channel → split equally ---
    sound_data[0, :] = _interp(gw_data[0], channel_weights[0])
    sound_data[1, :] = _interp(gw_data[1], channel_weights[1])
    tmp = _interp(gw_data[2], 0.5 * channel_weights[2])
    sound_data[0, :] += tmp
    sound_data[1, :] += tmp

    # --- Gaussian smoothing (same as MATLAB gausswin + conv) ---
    gw_len = int(round(sample_rate / 6000) * 2 + 1)
    gw = windows.gaussian(gw_len, std=(gw_len - 1) / 5)
    gw /= np.sum(gw)

    sound_data[0, :] = convolve(sound_data[0, :], gw, mode='same')
    sound_data[1, :] = convolve(sound_data[1, :], gw, mode='same')

    # --- normalize ---
    max_val = np.max(np.abs(sound_data))
    if max_val > 0:
        sound_data = 0.95 * sound_data / max_val

    # --- play sound ---
    if not only_produce_sound_data:
        duration = sound_length * dwell_time
        print(f'playing out the sequence waveform, duration {duration:.1f}s')

        if not _HAS_SOUNDDEVICE:
            raise RuntimeError('sounddevice not installed')

        pad = np.zeros((2, sample_rate // 2))
        play_data = np.hstack([pad, sound_data, pad]).T  # (N, 2)
        sd.play(play_data, samplerate=sample_rate)
        sd.wait()

    return sound_data
