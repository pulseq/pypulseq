# -*- coding: utf-8 -*-
"""
TODO
"""
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements

import numpy as np
import numpy.linalg as linalg
import collections
import warnings

# define Events
RFPulse = collections.namedtuple('rf', ['type', 'signal', 't', 'freq_offset',
                                        'phase_offset', 'dead_time',
                                        'ringdown_time', 'gz'])
Gradient = collections.namedtuple('grad', ['type', 'channel', 'amplitude',
                                           'rise_time', 'flat_time',
                                           'fall_time', 'area', 'flat_area',
                                           't', 'waveform'])
ADC = collections.namedtuple('adc', ['type', 'num_samples', 'dwell',
                                     'duration', 'delay', 'freq_offset',
                                     'phase_offset', 'dead_time'])
Delay = collections.namedtuple('delay', ['type', 'delay'])


def opts(grad_unit='mT/m', slew_unit='T/m/s', max_grad=40, max_slew=170,
         rise_time=None, rf_dead_time=None, rf_ringdown_time=None,
         adc_dead_time=None, rf_raster_time=None, grad_raster_time=None):
    """
    TODO
    """

    valid_grad_units = ['Hz/m', 'mT/m', 'rad/ms/mm']
    valid_slew_units = ['Hz/m/s', 'mT/m/ms', 'rad/ms/mm/ms']

    # #ifdef EXTERNAL_GRADS
    valid_grad_units.extend(['A', ''])
    valid_slew_units.extend(['A/s', '1/s'])
    # #endif

    grad_to_unit = 'Hz/m'
    if grad_unit is None:
        warnings.warn('No grad_unit given, will assume  mT/m')
        grad_unit = 'mT/m'
    # #ifdef EXTERNAL_GRADS
    if grad_unit == 'A':
        grad_to_unit = ''
    # #endif EXTERNAL_GRADS

    max_grad = convert(max_grad, grad_unit, grad_to_unit)

    slew_to_unit = 'Hz/m/s'
    if slew_unit is None:
        warnings.warn('No slew_unit given, will assume  T/m/s')
        slew_unit = 'T/m/s'

    # #ifdef EXTERNAL_GRADS
    if slew_unit == 'A/s':
        slew_to_unit = ''
    # #endif EXTERNAL_GRADS

    max_slew = convert(max_slew, slew_unit, slew_to_unit)

    if rise_time is not None:
        max_slew = None
    if rf_dead_time is None:
        rf_dead_time = 0
    if rf_ringdown_time is None:
        rf_ringdown_time = 0
    if adc_dead_time is None:
        adc_dead_time = 0
    if rf_raster_time is None:
        rf_raster_time = 1e-6
    if grad_raster_time is None:
        grad_raster_time = 10e-6

    opt = {'grad_unit':        grad_to_unit,
           'slew_unit':        slew_to_unit,
           'max_grad':         max_grad,
           'max_slew':         max_slew,
           'rise_time':        rise_time,
           'rf_dead_time':     rf_dead_time,
           'rf_ringdown_time': rf_ringdown_time,
           'adc_dead_time':    adc_dead_time,
           'rf_raster_time':   rf_raster_time,
           'grad_raster_time': grad_raster_time}

    return opt


def convert(val, from_unit, to_unit=None):
    """
    TODO
    """

    if from_unit is None and to_unit is None:
        raise ValueError('No units given.')

    valid_grad_units = ['Hz/m', 'mT/m', 'rad/ms/mm']
    valid_slew_units = ['Hz/m/s', 'mT/m/ms', 'rad/ms/mm/ms']

    # #ifdef EXTERNAL_GRADS
    valid_grad_units.extend(['A', ''])
    valid_slew_units.extend(['A/s', '1/s'])
    # #endif

    if from_unit is None:
        if to_unit in valid_grad_units:
            warnings.warn('No from_unit given, will assume gradient value \
                          with unit Hz/m')
            from_unit = 'Hz/m'
        if to_unit in valid_slew_units:
            warnings.warn('No from_unit given, will assume slew value \
                          with unit Hz/m/s')
            from_unit = 'Hz/m/s'


    if from_unit == to_unit:
        # no conversion necessary
        return val

    gamma = 42.57747892e6  # Hz/T

    # set default output unit if not given
    if to_unit is None:
        if any(from_unit in s for s in valid_grad_units):
            to_unit = valid_grad_units[0]
        elif any(from_unit in s for s in valid_slew_units):
            to_unit = valid_slew_units[0]

    # Grad units
    if from_unit == 'Hz/m':
        standard = val
    elif from_unit == 'mT/m':
        standard = val*1e-3*gamma
    elif from_unit == 'rad/ms/mm':
        standard = val*1e6/(2*np.pi)
    # Slew units
    elif from_unit == 'Hz/m/s':
        standard = val
    elif from_unit == 'mT/m/ms':
        standard = val*gamma
    elif from_unit == 'T/m/s':
        standard = val*gamma
    elif from_unit == 'rad/ms/mm/ms':
        standard = val*1e9/(2*np.pi)
    # #ifdef EXTERNAL_GRADS
    elif from_unit == '1/s':
        standard = val
    elif from_unit == '':
        standard = val
    elif from_unit == 'A':
        standard = val/150.0  # careful
    elif from_unit == 'A/s':
        standard = val/150.0  # careful
    # #endif

    # Grad units
    if to_unit == 'Hz/m':
        out = standard
    elif to_unit == 'mT/m':
        out = 1e3*standard/gamma
    elif to_unit == 'rad/ms/mm':
        out = standard*2*np.pi*1e-6
    # Slew units
    elif to_unit == 'Hz/m/s':
        out = standard
    elif to_unit == 'mT/m/ms':
        out = standard/gamma
    elif to_unit == 'T/m/s':
        out = standard/gamma
    elif to_unit == 'rad/ms/mm/ms':
        out = standard*2*np.pi*1e-9
    # #ifdef EXTERNAL_GRADS
    elif to_unit == '':
        out = standard
    elif to_unit == '1/s':
        out = standard
    elif to_unit == 'A':
        out = standard*150
    elif to_unit == 'A/s':
        out = standard*150
    # #endif EXTERNAL_GRADS

    return out


def add_ramps(k, system=None, rf=None, max_grad=None, max_slew=None):
    """
    Add segment to kspace trajectory to ramp to and from the given trajectory

    kout = add_ramps(k) Add a segment to k so that the output travels from 0 to
    k(1) and a segment so that the output goes from k(end) back to 0 without
    violating the gradient and slew constraints.

    [kx,ky,...] = add_ramps([kx,ky,...]) Add segments of the same length for
    each trajectory in the list.

    [...,rf] = add_ramps(...,rf=x) Add a segment of zeros over the ramp times
    to an RF shape.

    See also Sequence.makeAbitraryGrad

    @author Stefan Kroboth
    """

    if system is None:
        system = opts()
    if max_grad is not None:
        system['max_grad'] = max_grad
    if max_slew is not None:
        system['max_slew'] = max_slew

    kn = np.zeros((3, k[0].size))
    if type(k) is list:
        n_channels = len(k)
        for i in range(n_channels):
            kn[i, :] = k[i]
        k = kn
    else:
        n_channels = k.shape[0]
        kn[0:n_channels, :] = k
        kn[n_channels:, :] = np.zeros((3-n_channels, k.shape[1]))
        k = kn

    out = calc_ramp(np.zeros((3, 2)), k[:, 0:2], system)
    k_up = out.kout
    # ok1 = out.success
    out = calc_ramp(k[:, -2:], np.zeros((3, 2)), system)
    k_down = out.kout
    # ok2 = out.success

    k_up = np.concatenate((np.zeros((3, 2)), k_up), axis=1)
    k_down = np.concatenate((k_down, np.zeros((3, 1))), axis=1)

    k = np.concatenate((k_up, k, k_down), axis=1)
    out = []

    for i in range(n_channels):
        out.append(k[i, :])

    if rf is not None:
        rf = rf[:, None]
        out.append(np.concatenate((np.zeros((k_up.shape[1]*10, 1)),
                                   rf, np.zeros((k_down.shape[1]*10, 1))),
                                  axis=0))

    return out


def calc_ramp(k_0, kend, system, max_points=500, max_grad=None, max_slew=None):
    """
    The aim of calc_ramp is to join the points k_0 and kend
    in three-dimensional k-space in minimal time, observing the gradient and
    slew limits, and the gradient strength g_0 before k(0,2) and g_end after
    kend(:,1)

    In the context of a fixed gradient dwell time this is a discrete problem
    with a priori unknown number of discretization steps. Therefore calc_ramp
    tries out the optimization with 0 steps, then 1 step, and so on, until
    all conditions can be fulfilled, this yielding a short connection.

    N.B. The connection found this way is not necessarily always the shortest
    (there are some counterexamples) but still quite short.
    Improvements possible.

    Usage: [kout, success] = calc_ramp(k_0, kend, system, max_points,
                                      max_grad, max_slew)
    """

    if system is None:
        system = opts()
    if max_grad is None:
        max_grad = system['max_grad']
    if max_slew is None:
        max_slew = system['max_slew']

    grad_raster = system['grad_raster_time']

    g_0 = (k_0[:, 1] - k_0[:, 0])/grad_raster
    g_end = (kend[:, 1] - kend[:, 0])/grad_raster
    k_0 = k_0[:, 1]
    kend = kend[:, 0]

    # kout = np.zeros((3, 0))  # assigned but never used?

    success = False
    use_points = 0

    while success is False and use_points <= max_points:
        if (linalg.norm(g_0) > max_grad) or (linalg.norm(g_end) > max_grad):
            break
        out = joinleft0(k_0, kend, g_0, g_end, use_points, grad_raster,
                        max_grad, max_slew)
        # kout = out.kout  # assigned but never used
        success = out.success
        use_points += 1

    return out


def joinleft0(k_0, kend, g_0, g_end, use_points, grad_raster, max_grad,
              max_slew):
    """
    Add one k-space point close to k_0. Gradient and slew limits apply in
    total vector limited mode.

    Rationale:

    0. If use_points == 0 the recursion stops. If k_0 and kend can be joined
       in one gradDwell time, return success, else return "no success".

    1. Calculate optimal k-space point kopt that would lie on a straight
       line of N=use_points evenly spaced points to kend. If this kopt can be
       reached within gradient and slew limits, kopt is the solution of this
       function call.

    2. If kopt cannot be reached, calculate the gradient limited point kgl
       closest to kopt. If this piont can be reached in one gradDwell time
       without violating the slew limit, kgl is the solution of this
       function call.

    3. If kgl is not inside the slew limit, the slew limited point closest
       to kop, ksl, is calculated. If ksl is inside the gradient limit, ksl
       is the solution.

    4. If neither kgl nor ksl are possible find the point gklsl closest to
       kopt that satisfies both limits at the same time.

    5. Call joinright0 to obtain the other points starting with a point next
       to kend.
    """

    success = False
    out = collections.namedtuple('out', ['kout', 'success'])
    if use_points == 0:
        G = np.zeros((g_0.shape[0], 3))
        G[:, 0] = g_0
        G[:, 1] = (kend-k_0)/grad_raster
        G[:, 2] = g_end
        S = (G[:, 1:]-G[:, 0:-1])/grad_raster
        koutleft = np.zeros((3, 0))
        success = inside_limits(G, S, max_grad, max_slew)
        return out(koutleft, success)

    dk = (kend-k_0)/(use_points+1)
    kopt = k_0+dk
    g_opt = (kopt-k_0)/grad_raster
    s_opt = (g_opt-g_0)/grad_raster

    ok_g_opt = np.sum(np.power(g_opt, 2)) <= max_grad**2
    ok_s_opt = np.sum(np.power(s_opt, 2)) <= max_slew**2

    if ok_g_opt and ok_s_opt:
        k_left = kopt
    else:
        a = max_grad*grad_raster
        b = max_slew*grad_raster**2

        dkprol = g_0*grad_raster
        dkconn = dk-dkprol

        ksl = k_0 + dkprol + dkconn/linalg.norm(dkconn)*b
        g_sl = (ksl-k_0)/grad_raster
        ok_g_sl = (np.sum(np.power(g_sl, 2)) <= max_grad**2)

        kgl = k_0 + dk/linalg.norm(dk)*a
        g_gl = (kgl-k_0)/grad_raster
        s_gl = (g_gl-g_0)/grad_raster
        ok_s_gl = (np.sum(np.power(s_gl, 2)) <= max_slew**2)

        if ok_g_sl:
            k_left = ksl
        elif ok_s_gl:
            k_left = kgl
        else:
            c = linalg.norm(dkprol)
            c1 = (a**2-b**2+c**2)/(2*c)
            h = np.sqrt(a**2-c1**2)
            kglsl = k_0 + c1*dkprol/linalg.norm(dkprol)
            projondkprol = (kgl*dkprol.T) * dkprol/linalg.norm(dkprol)
            hdirection = kgl - projondkprol
            kglsl = kglsl + h*hdirection/linalg.norm(hdirection)
            k_left = kglsl

    k = joinright0(k_left, kend, (k_left-k_0)/grad_raster, g_end, use_points-1,
                   grad_raster, max_grad, max_slew)

    k_left = k_left[:, None]  # get dimension info back (sort of a quick hack)
    success = k.success
    return out(np.hstack((k_left, k.kout)), success)


def joinright0(k_0, kend, g_0, g_end, use_points, grad_raster, max_grad,
               max_slew):
    """
    Add one k-space point close to kend. Gradient and slew limits apply in
    total vector limited mode. Rationale see joinleft0.
    """

    success = False
    out = collections.namedtuple('out', ['kout', 'success'])
    if use_points == 0:
        G = np.zeros((g_0.shape[0], 3))
        G[:, 0] = g_0
        G[:, 1] = (kend-k_0)/grad_raster
        G[:, 2] = g_end
        S = (G[:, 1:]-G[:, 0:-1])/grad_raster
        koutright = np.zeros((3, 0))
        success = inside_limits(G, S, max_grad, max_slew)
        return out(koutright, success)

    dk = (k_0-kend)/(use_points+1)
    kopt = kend+dk
    g_opt = (kend-kopt)/grad_raster
    s_opt = (g_end-g_opt)/grad_raster

    ok_g_opt = np.sum(np.power(g_opt, 2)) <= max_grad**2
    ok_s_opt = np.sum(np.power(s_opt, 2)) <= max_slew**2

    if ok_g_opt and ok_s_opt:
        k_right = kopt
    else:
        a = max_grad*grad_raster
        b = max_slew*grad_raster**2

        dkprol = -g_end*grad_raster
        dkconn = dk-dkprol

        ksl = kend + dkprol + dkconn/linalg.norm(dkconn)*b
        g_sl = (kend-ksl)/grad_raster
        ok_g_sl = (np.sum(np.power(g_sl, 2)) <= max_grad**2)

        kgl = k_0 + dk/linalg.norm(dk)*a
        g_gl = (kend-kgl)/grad_raster
        s_gl = (g_end-g_gl)/grad_raster
        ok_s_gl = (np.sum(np.power(s_gl, 2)) <= max_slew**2)

        if ok_g_sl:
            k_right = ksl
        elif ok_s_gl:
            k_right = kgl
        else:
            c = linalg.norm(dkprol)
            c1 = (a**2-b**2+c**2)/(2*c)
            h = np.sqrt(a**2-c1**2)
            kglsl = kend + c1*dkprol/linalg.norm(dkprol)
            projondkprol = (kgl*dkprol.T) * dkprol/linalg.norm(dkprol)
            hdirection = kgl - projondkprol
            kglsl = kglsl + h*hdirection/linalg.norm(hdirection)
            k_right = kglsl

    k = joinleft0(k_0, k_right, g_0, (kend-k_right)/grad_raster, use_points-1,
                  grad_raster, max_grad, max_slew)
    k_right = k_right[:, None]  # get dimension info back (quick hack)
    success = k.success
    return out(np.hstack((k.kout, k_right)), success)


def inside_limits(grad, slew, max_grad, max_slew):
    """
    Check if both gradient and slew rates are inside the respective limits
    """

    grad2 = np.sum(np.power(grad, 2))
    slew2 = np.sum(np.power(slew, 2))
    return (np.max(grad2) <= max_grad**2) and (np.max(slew2) <= max_slew**2)


def traj2grad(k, system=None, grad_raster_time=None):
    """
    Convert a k-space trajectory to a gradient waveform using finite
    differences. The trajectory is in units of 1/m.
    The size of k is [n_channels, nTime].

    g = traj2grad(k, rasterTime=T) : Calculate gradient waveforms
    assuming the given raster time.

    See also Sequence.make_arbitrary_grad
    """

    if system is None:
        system = opts()
    if grad_raster_time is None:
        grad_raster_time = system['grad_raster_time']

    # Special case when k is a vector
    is_vec = k.ndim == 1
    if is_vec:
        k = k[:, None].T

    # compute finite difference for gradients in Hz/m
    out = np.concatenate((k[:, 1:]-k[:, 0:-1], np.zeros((k.shape[0], 1))),
                         axis=1)/grad_raster_time
    if is_vec:
        return out.flatten()
    else:
        return out


def make_arbitrary_rf(signal, flip_angle, system=None, freq_offset=0,
                      phase_offset=0, time_bw_product=None, bandwidth=None,
                      max_grad=None, max_slew=None, slice_thickness=None):
    """
    Create an RF pulse with the given pulse shape.

    If freq_offset and phase_offset are given, a block pulse with frequency
    offset and phase offset is created.

    If bandwith and slice_thickness are provided, an RF pulse and the
    corresponding slice select gradient is retuurned. The bandwith of the pulse
    must be given for the specified shape.

    See also Sequence.make_sinc_pulse, Sequence.addBlock
    """

    if system is None:
        system = opts()

    signal = signal/np.sum(signal*system['rf_raster_time'])*flip_angle/(2*np.pi)

    N = signal.size
    duration = N*system['rf_raster_time']
    t = np.arange(1, N+1)*system['rf_raster_time']

    if (slice_thickness is not None) and (bandwidth is not None):
        if max_grad is not None:
            system['max_grad'] = max_grad
        if max_slew is not None:
            system['max_slew'] = max_slew
        if time_bw_product is not None:
            bandwith = time_bw_product/duration

        amplitude = bandwith/slice_thickness
        area = amplitude*duration
        gz = make_trapezoid('z', system, flat_time=duration, flat_area=area)

        t_fill = np.arange(1, np.round(gz.rise_time/1e-6))*1e-6  # Round to mu
        t = np.concatenate((t_fill, t+t_fill[-1], t_fill+t[-1]+t_fill[-1]))
        signal = np.concatenate((np.zeros(t_fill.shape), signal,
                                 np.zeros(t_fill.shape)))
    else:
        gz = None

    if system['rf_ringdown_time'] > 0:
        # Round to mu s
        t_fill = np.arange(1, np.round(system['rf_ringdown_time']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+t_fill))
        signal = np.concatenate((signal, np.seros(t_fill.shape)))

    return RFPulse('rf', signal, t, freq_offset, phase_offset,
                   system['rf_dead_time'], system['rf_ringdown_time'], gz)


def make_sinc_pulse(flip_angle, system=None, duration=0, freq_offset=0,
                    phase_offset=0, time_bw_product=4, apodization=0,
                    max_grad=None, max_slew=None, slice_thickness=None):
    """
    Create a slice selective sinc pulse

    If duration is given: Create sinc pulse with given flip angle (rad) and
        duration(s)
    If freq_offset and phase_offset are given: Create sinc pulse with frequence
        offset (Hz) and phase offset (rad)
    If slice_thickness is given: Return the slice select gradient corresponding
        to given slice thickness (m)
    If system is given: Create slice selection gradient with the specified
        gradient limits (e.g. amplitude, slew). If not provided, default values
        will be used.

    See also Sequence.addBlock
    """

    if system is None:
        system = opts()

    BW = time_bw_product/np.float(duration)
    alpha = apodization
    N = np.round(duration/1.0e-6)
    t = (np.arange(N)+1)*system['rf_raster_time']
    tt = t - duration/2.0
    window = (1.0-alpha+alpha*np.cos(2*np.pi*tt/duration))
    signal = window * np.sinc(BW*tt)
    flip = np.sum(signal)*system['rf_raster_time']*2*np.pi
    signal = signal*flip_angle/flip

    fill_time = 0
    if slice_thickness is not None:
        if max_grad is not None:
            system['max_grad'] = max_grad
        if max_slew is not None:
            system['max_slew'] = max_slew

        amplitude = BW/np.float(slice_thickness)
        area = amplitude*duration
        gz = make_trapezoid('z', system, flat_time=duration, flat_area=area)

        # Pad RF pulse with zeros during gradient ramp up
        fill_time = gz.rise_time
        t_fill = (np.arange(np.round(fill_time/1e-6))+1)*1e-6  # Round to mu s
        t = np.concatenate((t_fill, t+t_fill[-1]))
        signal = np.concatenate((np.zeros(t_fill.shape), signal))
    else:
        gz = None

    # Add dead time to start of pulse, if required
    if fill_time < system['rf_dead_time']:
        fill_time = system['rf_dead_time'] - fill_time
        t_fill = (np.arange(np.round(fill_time/1e-6))+1)*1e-6  # round to mu s
        t = np.concatenate((t_fill, t+t_fill[-1]))
        signal = np.concatenate((np.zeros(t_fill.shape), signal))

    if system['rf_ringdown_time'] > 0:
        # Round to mu s
        t_fill = np.arange(1, np.round(system['rf_ringdown_time']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+t_fill))
        signal = np.concatenate((signal, np.seros(t_fill.shape)))

    return RFPulse('rf', signal, t, freq_offset, phase_offset,
                   system['rf_dead_time'], system['rf_ringdown_time'], gz)


def make_trapezoid(channel, system=None, duration=0, area=None, flat_time=0,
                   flat_area=None, amplitude=None, max_grad=None,
                   max_slew=None, rise_time=None):
    """
    Create trapezoid gradient with the specified gradient limits
    (e.g. amplitude, slew).

    If duration and Area are given: Create a trapezoid gradient with given
    duration (s) and a total area (1/m) including ramps.

    If flat_time and flat_area are given: Create a trapezoid gradient with
    given flat-top time and flat-top area not including ramps

    If amplitude is given: Create a trapezoid gradient with given amplitude
    (Hz/m)

    See also Sequence.addblock and mr.opts
    """

    valid_channels = ['x', 'y', 'z']
    # #ifdef EXTERNAL_GRADS
    valid_channels = valid_channels + [str(x) for x in range(1, 13)]
    # #endif

    if system is None:
        system = opts()
    if max_grad is not None:
        system['max_grad'] = max_grad
    if max_slew is not None:
        system['max_slew'] = max_slew
    if rise_time is not None:
        system['rise_time'] = rise_time
    if system['rise_time'] is not None:
        rise_time = system['rise_time']

    amplitude_none = amplitude is None

    if (area is None) and (flat_area is None) and (amplitude is None):
        raise Exception('Must supply either \'area\', \'flat_area\' or' +
                        ' \'amplitude\'')

    if flat_time > 0:
        if amplitude is None:
            amplitude = flat_area/flat_time

        if rise_time is None:
            rise_time = np.abs(amplitude)/system['max_slew']
            rise_time = float(np.ceil(rise_time/system['grad_raster_time'])) * \
                system['grad_raster_time']
            system['rise_time'] = rise_time

        fall_time = rise_time
    elif duration is not None:
        if amplitude is None:
            if rise_time is None:
                dC = 1/np.abs(2*system['max_slew']) + \
                    1/np.abs(2*system['max_slew'])
                possible = duration**2 > 4*abs(area)*dC
                amplitude = (duration -
                             np.sqrt(duration**2-4*np.abs(area)*dC))/(2*dC)
            else:
                amplitude = area/np.float64(duration-rise_time)
                possible = (duration > 2*rise_time) & \
                    (np.abs(amplitude) < system['max_grad'])

            if not possible:
                raise Exception('Requested area is too large for this ' +
                                'gradient.')
        if rise_time is None:
            rise_time = np.ceil(amplitude/system['max_slew'] /
                                system['grad_raster_time']) *\
                system['grad_raster_time']
            system['rise_time'] = rise_time

        fall_time = rise_time
        flat_time = duration - rise_time - fall_time

        if amplitude_none:
            # Adjust amplitude (after rounding) to achieve given area
            amplitude = area/(rise_time/2+fall_time/2+flat_time)
    else:
        raise Exception('Must supply a duration.')

    if np.abs(amplitude) > system['max_grad']:
        raise Exception('Amplitude violation.')

    return Gradient('trap', channel, amplitude, rise_time, flat_time,
                    fall_time, amplitude*(flat_time +
                                          system['rise_time']/2.0 +
                                          fall_time/2.0),
                    flat_area, None, None)


def make_arbitrary_grad(channel, waveform, system=None, max_grad=None,
                        max_slew=None):
    """
    Create an gradient event with arbitrary waveform satisfying gradient
    hardware constraints.

    See also Sequence.addBlock
    """

    valid_channels = ['x', 'y', 'z']
    # #ifdef EXTERNAL_GRADS
    valid_channels = valid_channels + [str(x) for x in range(1, 13)]
    # #endif

    if system is None:
        system = opts()
    if max_grad is not None:
        system['max_grad'] = max_grad
    if max_slew is not None:
        system['max_slew'] = max_slew

    slew = (waveform[1:]-waveform[0:-1])/system['grad_raster_time']

    if np.max(np.abs(slew)) > system['max_slew']:
        raise Exception('Slew rate violation (' +
                        str(np.max(np.abs(slew))/system['max_slew']*100) +
                        '%)')
    if np.max(np.abs(waveform)) > system['max_grad']:
        raise Exception('Gradient amplitude violation (' +
                        str(np.max(np.abs(waveform))/system['max_grad']*100) +
                        '%)')

    return Gradient('grad', channel, None, None, None, None, None, None,
                    np.arange(len(waveform))*system['grad_raster_time'],
                    waveform)


def make_adc(num_samples, system=None, dwell=0, duration=0, delay=0,
             freq_offset=0, phase_offset=0):
    """
    Create and ADC readout event.

    If dwell is given: Create ADC with num_samples samples with givend dwell
        time

    If duration is given: Create ADC with num_samples and specified total
        duration

    If delay is given: Create ADC with initial delay

    If system is given: Create ADC considering system properties given in
            system. For example, a dead time after sampling can be added to the
            duration

    TODO: System limits not really satisified! Apart from dead_time... (SK)

    See also Sequence.addBlock
    """

    if system is None:
        system = opts()

    if dwell < 0:
        raise Exception('dwell must be positive')

    if duration < 0:
        raise Exception('duration must be positive')

    if ((dwell == 0) and (duration == 0)) or \
       ((np.abs(dwell) > 0) and (np.abs(duration) > 0)):
        raise Exception('Either dwell or duration must be defined')

    if duration > 0:
        dwell = duration/num_samples

    if dwell > 0:
        duration = dwell*num_samples

    return ADC('adc', num_samples, dwell, duration, delay, freq_offset,
               phase_offset, system['adc_dead_time'])


def make_block_pulse(flip_angle, duration=0, system=None, freq_offset=0,
                     phase_offset=0, time_bw_product=0, bandwidth=None,
                     max_grad=None, max_slew=None, slice_thickness=None):
    """
    Create an RF pulse with the given pulse shape.

    If freq_offset and phase_offset are given, a block pulse with frequency
    offset and phase offset is created.

    If bandwith and slice_thickness are provided, an RF pulse and the
    corresponding slice select gradient is retuurned. The bandwith of the pulse
    must be given for the specified shape.

    See also Sequence.make_sinc_pulse, Sequence.addBlock
    """

    if system is None:
        system = opts()

    if duration == 0:
        if time_bw_product > 0:
            duration = time_bw_product/bandwidth
        elif bandwidth > 0:
            duration = 1/(4*bandwidth)
        else:
            raise Exception('Either bandwidth or duration must be defined.')

    BW = 1/(4*duration)
    N = np.round(duration/1e-6)
    t = (np.arange(N)+1)*system['rf_raster_time']
    signal = flip_angle/(2*np.pi)/duration*np.ones(t.shape)

    fill_time = 0
    if slice_thickness is not None:
        if max_grad > 0:
            system['max_grad'] = max_grad
        if max_slew > 0:
            system['max_slew'] = max_slew

        amplitude = BW/slice_thickness
        area = amplitude*duration
        gz = make_trapezoid('z', system, flat_time=duration, flat_area=area)

        fill_time = gz.rise_time
        t_fill = (np.arange(np.round(fill_time/1e-6))+1)*1e-6
        t = np.concatenate((t_fill, t+t_fill[-1], t_fill+t[-1]+t_fill[-1]))
        signal = np.concatenate((np.zeros(t_fill.shape), signal,
                                 np.zeros(t_fill.shape)))
    else:
        gz = None

    if fill_time < system['rf_dead_time']:
        fill_time = system['rf_dead_time'] - fill_time
        t_fill = (np.arange(fill_time/1e-6)+1)*1e-6  # round to microsecond
        t = np.concatenate((t_fill, t+t_fill[-1]))
        signal = np.concatenate((np.zeros(t_fill.shape), signal))

    if system['rf_ringdown_time'] > 0:
        # Round to mu s
        t_fill = np.arange(1, np.round(system['rf_ringdown_time']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+t_fill))
        signal = np.concatenate((signal, np.seros(t_fill.shape)))

    return RFPulse('rf', signal, t, freq_offset, phase_offset,
                   system['rf_dead_time'], system['rf_ringdown_time'], gz)


def calc_duration(blocks):
    """
    Calculate the duration of an event or block
    Determine the maximum duration of t he provided events
    """

    if blocks is not list:
        blocks = [blocks]

    duration = 0
    for block in blocks:
        if block.type == 'delay':
            duration = max(duration, block.delay)
        elif block.type == 'rf':
            duration = max(duration, block.t[-1]+block.dead_time +
                           block.ringdown_time)
        elif block.type == 'grad':
            duration = max(duration, block.t[-1])
        elif block.type == 'adc':
            duration = max(duration, block.delay +
                           block.num_samples*block.dwell + block.dead_time)
        elif block.type == 'trap':
            duration = max(duration, block.rise_time + block.flat_time +
                           block.fall_time)

    return duration


def compress_shape(w):
    """
    Compress gradient or pulse shape using a run-length compression
    scheme on the derivative. This strategy encodes constant and linear
    waveforms with very few samples. A structure is returned with the fields:
        num_samples: The number of samples in the uncompressed waveform
        data: contining the compressed waveform

    See also decompress_shape
    """
    data = np.concatenate((w[0].flatten(), np.diff(w.flatten(), axis=0)))

    # Mask is TRUE if values change
    mask_changes = np.concatenate(([True], np.abs(np.diff(data)) > 1e-8))
    vals = data[mask_changes]  # Elements without repetitions
    # Indices of changes
    k = np.where(np.concatenate((mask_changes, [True])))[0]
    n = np.diff(k.flatten())  # Number of repetitions

    # Encode in Pulseq format
    n_extra = n-2.0
    vals2 = vals.copy()
    vals2[n_extra < 0] = np.nan
    n_extra[n_extra < 0] = np.nan
    v = np.concatenate((vals[:, None], vals2[:, None], n_extra[:, None]),
                       axis=1)
    v = v[np.isfinite(v)]
    v[np.abs(v) < 1e-10] = 0

    Shape = collections.namedtuple('shape', ['num_samples', 'data'])
    return Shape(w.size, v.flatten())


def make_delay(delay):
    """
    Create a delay event with a given delay.

    See also Sequence.addBlock
    """

    if (not np.isfinite(delay)) or (delay <= 0):
        raise Exception('Delay (' + str(delay*1e3) + 'ms) is invalid.')

    return Delay('delay', delay)  # delay, delay, delaaaayyy :D
