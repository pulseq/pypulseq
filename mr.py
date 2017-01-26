# -*- coding: utf-8 -*-
"""
TODO
"""
# pylint: disable=too-many-instance-attributes
# pylint: disable=invalid-name
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments

import numpy as np
import numpy.linalg as linalg
import collections

# define Events
RFPulse = collections.namedtuple('rf', ['type', 'signal', 't', 'freqOffset',
                                        'phaseOffset', 'deadTime',
                                        'ringdownTime', 'gz'])
Gradient = collections.namedtuple('grad', ['type', 'channel', 'amplitude',
                                           'riseTime', 'flatTime', 'fallTime',
                                           'area', 'flatArea', 't',
                                           'waveform'])
ADC = collections.namedtuple('adc', ['type', 'numSamples', 'dwell', 'duration',
                                     'delay', 'freqOffset', 'phaseOffset',
                                     'deadTime'])
Delay = collections.namedtuple('delay', ['type', 'delay'])


def opts(gradUnit=None, slewUnit=None, maxGrad=None, maxSlew=None,
         riseTime=None, rfDeadTime=None, rfRingdownTime=None, adcDeadTime=None,
         rfRasterTime=None, gradRasterTime=None):
    """
    TODO
    """

    validGradUnits = ['Hz/m', 'mT/m', 'rad/ms/mm']
    validSlewUnits = ['Hz/m/s', 'mT/m/ms', 'rad/ms/mm/ms']

    # #ifdef EXTERNAL_GRADS
    validGradUnits = validGradUnits.extend(['A', ''])
    validSlewUnits = validSlewUnits.extend(['A/s', '1/s'])
    # #endif

    gradToUnit = 'Hz/m'
    if maxGrad is None:
        maxGrad = convert(40, 'mT/m', gradToUnit)  # Default: 40 mT/m
    else:
        # #ifdef EXTERNAL_GRADS
        if gradUnit == 'A/s':
            maxGrad = convert(maxGrad, gradUnit, '')
        # #endif EXTERNAL_GRADS
        maxGrad = convert(maxGrad, gradUnit, gradToUnit)
        gradUnit = gradToUnit

    slewToUnit = 'Hz/m/s'
    if maxSlew is None:
        maxSlew = convert(170, 'T/m/s', slewToUnit)
        slewUnit = slewToUnit
    else:
        # #ifdef EXTERNAL_GRADS
        if slewUnit == 'A/s':
            slewToUnit = ''
        # #endif EXTERNAL_GRADS
        maxSlew = convert(maxSlew, slewUnit, slewToUnit)
        slewUnit = slewToUnit

    if riseTime is not None:
        maxSlew = None
    if rfDeadTime is None:
        rfDeadTime = 0
    if rfRingdownTime is None:
        rfRingdownTime = 0
    if adcDeadTime is None:
        adcDeadTime = 0
    if rfRasterTime is None:
        rfRasterTime = 1e-6
    if gradRasterTime is None:
        gradRasterTime = 10e-6

    opt = {'gradUnit':       gradUnit,
           'slewUnit':       slewUnit,
           'maxGrad':        maxGrad,
           'maxSlew':        maxSlew,
           'riseTime':       riseTime,
           'rfDeadTime':     rfDeadTime,
           'rfRingdownTime': rfRingdownTime,
           'adcDeadTime':    adcDeadTime,
           'rfRasterTime':   rfRasterTime,
           'gradRasterTime': gradRasterTime}

    return opt


def convert(val, fromUnit, toUnit=None):
    """
    TODO
    """

    validGradUnits = ['Hz/m', 'mT/m', 'rad/ms/mm']
    validSlewUnits = ['Hz/m/s', 'mT/m/ms', 'rad/ms/mm/ms']

    # #ifdef EXTERNAL_GRADS
    validGradUnits = validGradUnits.extend(['A', ''])
    validSlewUnits = validSlewUnits.extend(['A/s', '1/s'])
    # #endif

    gamma = 42.576e6  # Hz/T

    # set default output unit if not given
    if toUnit is None:
        if any(fromUnit in s for s in validGradUnits):
            toUnit = validGradUnits[0]
        elif any(fromUnit in s for s in validSlewUnits):
            toUnit = validSlewUnits[0]

    # Grad units
    if fromUnit == 'Hz/m':
        standard = val
    elif fromUnit == 'mT/m':
        standard = val*1e-3*gamma
    elif fromUnit == 'rad/ms/mm':
        standard = val*1e6/(2*np.pi)
    # Slew units
    elif fromUnit == 'Hz/m/s':
        standard = val
    elif fromUnit == 'mT/m/ms':
        standard = val*gamma
    elif fromUnit == 'T/m/s':
        standard = val*gamma
    elif fromUnit == 'rad/ms/mm/ms':
        standard = val*1e9/(2*np.pi)
    # #ifdef EXTERNAL_GRADS
    elif fromUnit == '1/s':
        standard = val
    elif fromUnit == '':
        standard = val
    elif fromUnit == 'A':
        standard = val/150  # careful
    elif fromUnit == 'A/s':
        standard = val/150  # careful
    # #endif

    # Grad units
    if toUnit == 'Hz/m':
        out = standard
    elif toUnit == 'mT/m':
        out = 1e3*standard/gamma
    elif toUnit == 'rad/ms/mm':
        out = standard*2*np.pi*1e-6
    # Slew units
    elif toUnit == 'Hz/m/s':
        out = standard
    elif toUnit == 'mT/m/ms':
        out = standard/gamma
    elif toUnit == 'T/m/s':
        out = standard/gamma
    elif toUnit == 'rad/ms/mm/ms':
        out = standard*2*np.pi*1e-9
    # #ifdef EXTERNAL_GRADS
    elif toUnit == '':
        out = standard
    elif toUnit == '1/s':
        out = standard
    elif toUnit == 'A':
        out = standard*150
    elif toUnit == 'A/s':
        out = standard*150
    # #endif EXTERNAL_GRADS

    return out


def addRamps(k, system=None, rf=None, maxGrad=None, maxSlew=None):
    """
    Add segment to kspace trajectory to ramp to and from the given trajectory

    kout = addRamps(k) Add a segment to k so that the output travels from 0 to
    k(1) and a segment so that the output goes from k(end) back to 0 without
    violating the gradient and slew constraints.

    [kx,ky,...] = addRamps([kx,ky,...]) Add segments of the same length for
    each trajectory in the list.

    [...,rf] = addRamps(...,rf=x) Add a segment of zeros over the ramp times
    to an RF shape.

    See also Sequence.makeAbitraryGrad

    @author Stefan Kroboth
    """

    if system is None:
        system = opts()
    if maxGrad is not None:
        system['maxGrad'] = maxGrad
    if maxSlew is not None:
        system['maxSlew'] = maxSlew

    kn = np.zeros((3, k[0].size))
    if type(k) is list:
        nChannels = len(k)
        for i in range(nChannels):
            kn[i, :] = k[i]
        k = kn
    else:
        nChannels = k.shape[0]
        kn[0:nChannels, :] = k
        kn[nChannels:, :] = np.zeros((3-nChannels, k.shape[1]))
        k = kn

    out = calcRamp(np.zeros((3, 2)), k[:, 0:2], system)
    kUp = out.kout
    # ok1 = out.success
    out = calcRamp(k[:, -2:], np.zeros((3, 2)), system)
    kDown = out.kout
    # ok2 = out.success

    kUp = np.concatenate((np.zeros((3, 2)), kUp), axis=1)
    kDown = np.concatenate((kDown, np.zeros((3, 1))), axis=1)

    k = np.concatenate((kUp, k, kDown), axis=1)
    out = []

    for i in range(nChannels):
        out.append(k[i, :])

    if rf is not None:
        rf = rf[:, None]
        out.append(np.concatenate((np.zeros((kUp.shape[1]*10, 1)),
                                   rf, np.zeros((kDown.shape[1]*10, 1))),
                                  axis=0))

    return out


def calcRamp(k0, kend, system, maxPoints=500, maxGrad=None, maxSlew=None):
    """
    The aim of calcRamp is to join the points k0 and kend in three-dimensional
    k-space in minimal time, observing the gradient and slew limits, and the
    gradient strength G0 before k(0,2) and Gend after kend(:,1)

    In the context of a fixed gradient dwell time this is a discrete problem
    with a priori unknown number of discretization steps. Therefore calcRamp
    tries out the optimization with 0 steps, then 1 step, and so on, until
    all conditions can be fulfilled, this yielding a short connection.

    N.B. The connection found this way is not necessarily always the shortest
    (there are some counterexamples) but still quite short.
    Improvements possible.

    Usage: [kout, success] = calcRamp(k0, kend, system, maxPoints,
                                      maxGrad, maxSlew)
    """

    if system is None:
        system = opts()
    if maxGrad is None:
        maxGrad = system['maxGrad']
    if maxSlew is None:
        maxSlew = system['maxSlew']

    gradRaster = system['gradRasterTime']

    G0 = (k0[:, 1] - k0[:, 0])/gradRaster
    Gend = (kend[:, 1] - kend[:, 0])/gradRaster
    k0 = k0[:, 1]
    kend = kend[:, 0]

    # kout = np.zeros((3, 0))  # assigned but never used?

    success = False
    usePoints = 0

    while success is False and usePoints <= maxPoints:
        if (linalg.norm(G0) > maxGrad) or (linalg.norm(Gend) > maxGrad):
            break
        out = joinleft0(k0, kend, G0, Gend, usePoints, gradRaster, maxGrad,
                        maxSlew)
        # kout = out.kout  # assigned but never used
        success = out.success
        usePoints += 1

    return out


def joinleft0(k0, kend, G0, Gend, usePoints, gradRaster, maxGrad, maxSlew):
    """
    Add one k-space point close to k0. Gradient and slew limits apply in
    total vector limited mode.

    Rationale:

    0. If usePoints == 0 the recursion stops. If k0 and kend can be joined
       in one gradDwell time, return success, else return "no success".

    1. Calculate optimal k-space point kopt that would lie on a straight
       line of N=usePoints evenly spaced points to kend. If this kopt can be
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
    if usePoints == 0:
        G = np.zeros((G0.shape[0], 3))
        G[:, 0] = G0
        G[:, 1] = (kend-k0)/gradRaster
        G[:, 2] = Gend
        S = (G[:, 1:]-G[:, 0:-1])/gradRaster
        koutleft = np.zeros((3, 0))
        success = InsideLimits(G, S, maxGrad, maxSlew)
        return out(koutleft, success)

    dk = (kend-k0)/(usePoints+1)
    kopt = k0+dk
    Gopt = (kopt-k0)/gradRaster
    Sopt = (Gopt-G0)/gradRaster

    okGopt = np.sum(np.power(Gopt, 2)) <= maxGrad**2
    okSopt = np.sum(np.power(Sopt, 2)) <= maxSlew**2

    if okGopt and okSopt:
        kLeft = kopt
    else:
        a = maxGrad*gradRaster
        b = maxSlew*gradRaster**2

        dkprol = G0*gradRaster
        dkconn = dk-dkprol

        ksl = k0 + dkprol + dkconn/linalg.norm(dkconn)*b
        Gsl = (ksl-k0)/gradRaster
        okGsl = (np.sum(np.power(Gsl, 2)) <= maxGrad**2)

        kgl = k0 + dk/linalg.norm(dk)*a
        Ggl = (kgl-k0)/gradRaster
        Sgl = (Ggl-G0)/gradRaster
        okSgl = (np.sum(np.power(Sgl, 2)) <= maxSlew**2)

        if okGsl:
            kLeft = ksl
        elif okSgl:
            kLeft = kgl
        else:
            c = linalg.norm(dkprol)
            c1 = (a**2-b**2+c**2)/(2*c)
            h = np.sqrt(a**2-c1**2)
            kglsl = k0 + c1*dkprol/linalg.norm(dkprol)
            projondkprol = (kgl*dkprol.T) * dkprol/linalg.norm(dkprol)
            hdirection = kgl - projondkprol
            kglsl = kglsl + h*hdirection/linalg.norm(hdirection)
            kLeft = kglsl

    k = joinright0(kLeft, kend, (kLeft-k0)/gradRaster, Gend, usePoints-1,
                   gradRaster, maxGrad, maxSlew)

    kLeft = kLeft[:, None]  # get dimension info back (sort of a quick hack)
    success = k.success
    return out(np.hstack((kLeft, k.kout)), success)


def joinright0(k0, kend, G0, Gend, usePoints, gradRaster, maxGrad, maxSlew):
    """
    Add one k-space point close to kend. Gradient and slew limits apply in
    total vector limited mode. Rationale see joinleft0.
    """

    success = False
    out = collections.namedtuple('out', ['kout', 'success'])
    if usePoints == 0:
        G = np.zeros((G0.shape[0], 3))
        G[:, 0] = G0
        G[:, 1] = (kend-k0)/gradRaster
        G[:, 2] = Gend
        S = (G[:, 1:]-G[:, 0:-1])/gradRaster
        koutright = np.zeros((3, 0))
        success = InsideLimits(G, S, maxGrad, maxSlew)
        return out(koutright, success)

    dk = (k0-kend)/(usePoints+1)
    kopt = kend+dk
    Gopt = (kend-kopt)/gradRaster
    Sopt = (Gend-Gopt)/gradRaster

    okGopt = np.sum(np.power(Gopt, 2)) <= maxGrad**2
    okSopt = np.sum(np.power(Sopt, 2)) <= maxSlew**2

    if okGopt and okSopt:
        kRight = kopt
    else:
        a = maxGrad*gradRaster
        b = maxSlew*gradRaster**2

        dkprol = -Gend*gradRaster
        dkconn = dk-dkprol

        ksl = kend + dkprol + dkconn/linalg.norm(dkconn)*b
        Gsl = (kend-ksl)/gradRaster
        okGsl = (np.sum(np.power(Gsl, 2)) <= maxGrad**2)

        kgl = k0 + dk/linalg.norm(dk)*a
        Ggl = (kend-kgl)/gradRaster
        Sgl = (Gend-Ggl)/gradRaster
        okSgl = (np.sum(np.power(Sgl, 2)) <= maxSlew**2)

        if okGsl:
            kRight = ksl
        elif okSgl:
            kRight = kgl
        else:
            c = linalg.norm(dkprol)
            c1 = (a**2-b**2+c**2)/(2*c)
            h = np.sqrt(a**2-c1**2)
            kglsl = kend + c1*dkprol/linalg.norm(dkprol)
            projondkprol = (kgl*dkprol.T) * dkprol/linalg.norm(dkprol)
            hdirection = kgl - projondkprol
            kglsl = kglsl + h*hdirection/linalg.norm(hdirection)
            kRight = kglsl

    k = joinleft0(k0, kRight, G0, (kend-kRight)/gradRaster, usePoints-1,
                  gradRaster, maxGrad, maxSlew)
    kRight = kRight[:, None]  # get dimension info back (sort of a quick hack)
    success = k.success
    return out(np.hstack((k.kout, kRight)), success)


def InsideLimits(grad, slew, maxGrad, maxSlew):
    """
    Check if both gradient and slew rates are inside the respective limits
    """

    grad2 = np.sum(np.power(grad, 2))
    slew2 = np.sum(np.power(slew, 2))
    return (np.max(grad2) <= maxGrad**2) and (np.max(slew2) <= maxSlew**2)


def traj2grad(k, system=None, gradRasterTime=None):
    """
    Convert a k-space trajectory to a gradient waveform using finite
    differences. The trajectory is in units of 1/m.
    The size of k is [nChannels, nTime].

    g = traj2grad(k, rasterTime=T) : Calculate gradient waveforms
    assuming the given raster time.

    See also Sequence.makeArbitraryGrad
    """

    if system is None:
        system = opts()
    if gradRasterTime is None:
        gradRasterTime = system['gradRasterTime']

    # Special case when k is a vector
    isVec = k.ndim == 1
    if isVec:
        k = k[:, None].T

    # compute finite difference for gradients in Hz/m
    out = np.concatenate((k[:, 1:]-k[:, 0:-1], np.zeros((k.shape[0], 1))),
                         axis=1)/gradRasterTime
    if isVec:
        return out.flatten()
    else:
        return out


def makeArbitraryRf(signal, flipAngle, system=None, freqOffset=0,
                    phaseOffset=0, timeBwProduct=None, bandwidth=None,
                    maxGrad=None, maxSlew=None, sliceThickness=None):
    """
    Create an RF pulse with the given pulse shape.

    If freqOffset and phaseOffset are given, a block pulse with frequency
    offset and phase offset is created.

    If bandwith and sliceThickness are provided, an RF pulse and the
    corresponding slice select gradient is retuurned. The bandwith of the pulse
    must be given for the specified shape.

    See also Sequence.makeSincPulse, Sequence.addBlock
    """

    if system is None:
        system = opts()

    signal = signal/np.sum(signal*system['rfRasterTime'])*flipAngle/(2*np.pi)

    N = signal.size
    duration = N*system['rfRasterTime']
    t = np.arange(1, N+1)*system['rfRasterTime']

    if (sliceThickness is not None) and (bandwidth is not None):
        if maxGrad is not None:
            system['maxGrad'] = maxGrad
        if maxSlew is not None:
            system['maxSlew'] = maxSlew
        if timeBwProduct is not None:
            bandwith = timeBwProduct/duration

        amplitude = bandwith/sliceThickness
        area = amplitude*duration
        gz = makeTrapezoid('z', system, flatTime=duration, flatArea=area)

        tFill = np.arange(1, np.round(gz.riseTime/1e-6))*1e-6  # Round to mu s
        t = np.concatenate((tFill, t+tFill[-1], tFill+t[-1]+tFill[-1]))
        signal = np.concatenate((np.zeros(tFill.shape), signal,
                                 np.zeros(tFill.shape)))
    else:
        gz = None

    if system['rfRingdownTime'] > 0:
        # Round to mu s
        tFill = np.arange(1, np.round(system['rfRingdownTime']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+tFill))
        signal = np.concatenate((signal, np.seros(tFill.shape)))

    return RFPulse('rf', signal, t, freqOffset, phaseOffset,
                   system['rfDeadTime'], system['rfRingdownTime'], gz)


def makeSincPulse(flipAngle, system=None, duration=0, freqOffset=0,
                  phaseOffset=0, timeBwProduct=4, apodization=0, maxGrad=None,
                  maxSlew=None, sliceThickness=None):
    """
    Create a slice selective sinc pulse

    If duration is given: Create sinc pulse with given flip angle (rad) and
        duration(s)
    If freqOffset and phaseOffset are given: Create sinc pulse with frequence
        offset (Hz) and phase offset (rad)
    If sliceThickness is given: Return the slice select gradient corresponding
        to given slice thickness (m)
    If system is given: Create slice selection gradient with the specified
        gradient limits (e.g. amplitude, slew). If not provided, default values
        will be used.

    See also Sequence.addBlock
    """

    if system is None:
        system = opts()

    BW = timeBwProduct/np.float(duration)
    alpha = apodization
    N = np.round(duration/1.0e-6)
    t = (np.arange(N)+1)*system['rfRasterTime']
    tt = t - duration/2.0
    window = (1.0-alpha+alpha*np.cos(2*np.pi*tt/duration))
    signal = window * np.sinc(BW*tt)
    flip = np.sum(signal)*system['rfRasterTime']*2*np.pi
    signal = signal*flipAngle/flip

    fillTime = 0
    if sliceThickness is not None:
        if maxGrad is not None:
            system['maxGrad'] = maxGrad
        if maxSlew is not None:
            system['maxSlew'] = maxSlew

        amplitude = BW/np.float(sliceThickness)
        area = amplitude*duration
        gz = makeTrapezoid('z', system, flatTime=duration, flatArea=area)

        # Pad RF pulse with zeros during gradient ramp up
        fillTime = gz.riseTime
        tFill = (np.arange(np.round(fillTime/1e-6))+1)*1e-6  # Round to mu s
        t = np.concatenate((tFill, t+tFill[-1]))
        signal = np.concatenate((np.zeros(tFill.shape), signal))
    else:
        gz = None

    # Add dead time to start of pulse, if required
    if fillTime < system['rfDeadTime']:
        fillTime = system['rfDeadTime'] - fillTime
        tFill = (np.arange(np.round(fillTime/1e-6))+1)*1e-6  # round to mu s
        t = np.concatenate((tFill, t+tFill[-1]))
        signal = np.concatenate((np.zeros(tFill.shape), signal))

    if system['rfRingdownTime'] > 0:
        # Round to mu s
        tFill = np.arange(1, np.round(system['rfRingdownTime']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+tFill))
        signal = np.concatenate((signal, np.seros(tFill.shape)))

    return RFPulse('rf', signal, t, freqOffset, phaseOffset,
                   system['rfDeadTime'], system['rfRingdownTime'], gz)


def makeTrapezoid(channel, system=None, duration=0, area=None, flatTime=0,
                  flatArea=None, amplitude=None, maxGrad=None, maxSlew=None,
                  riseTime=None):
    """
    Create trapezoid gradient with the specified gradient limits
    (e.g. amplitude, slew).

    If duration and Area are given: Create a trapezoid gradient with given
    duration (s) and a total area (1/m) including ramps.

    If flatTime and flatArea are given: Create a trapezoid gradient with given
    flat-top time and flat-top area not including ramps

    If amplitude is given: Create a trapezoid gradient with given amplitude
    (Hz/m)

    See also Sequence.addblock and mr.opts
    """

    validChannels = ['x', 'y', 'z']
    # #ifdef EXTERNAL_GRADS
    validChannels = validChannels + [str(x) for x in range(1, 13)]
    # #endif

    if system is None:
        system = opts()
    if maxGrad is not None:
        system['maxGrad'] = maxGrad
    if maxSlew is not None:
        system['maxSlew'] = maxSlew
    if riseTime is not None:
        system['riseTime'] = riseTime
    if system['riseTime'] is not None:
        riseTime = system['riseTime']

    amplitude_none = amplitude is None

    if (area is None) and (flatArea is None) and (amplitude is None):
        raise Exception('Must supply either \'area\', \'flatArea\' or' +
                        ' \'amplitude\'')

    if flatTime > 0:
        if amplitude is None:
            amplitude = flatArea/flatTime

        if riseTime is None:
            riseTime = np.abs(amplitude)/system['maxSlew']
            riseTime = float(np.ceil(riseTime/system['gradRasterTime'])) * \
                system['gradRasterTime']
            system['riseTime'] = riseTime

        fallTime = riseTime
    elif duration is not None:
        if amplitude is None:
            if riseTime is None:
                dC = 1/np.abs(2*system['maxSlew']) + \
                    1/np.abs(2*system['maxSlew'])
                possible = duration**2 > 4*abs(area)*dC
                amplitude = (duration -
                             np.sqrt(duration**2-4*np.abs(area)*dC))/(2*dC)
            else:
                amplitude = area/np.float64(duration-riseTime)
                possible = (duration > 2*riseTime) & \
                    (np.abs(amplitude) < system['maxGrad'])

            if not possible:
                raise Exception('Requested area is too large for this ' +
                                'gradient.')
        if riseTime is None:
            riseTime = np.ceil(amplitude/system['maxSlew'] /
                               system['gradRasterTime']) *\
                system['gradRasterTime']
            system['riseTime'] = riseTime

        fallTime = riseTime
        flatTime = duration - riseTime - fallTime

        if amplitude_none:
            # Adjust amplitude (after rounding) to achieve given area
            amplitude = area/(riseTime/2+fallTime/2+flatTime)
    else:
        raise Exception('Must supply a duration.')

    if np.abs(amplitude) > system['maxGrad']:
        raise Exception('Amplitude violation.')

    return Gradient('trap', channel, amplitude, riseTime, flatTime, fallTime,
                    amplitude*(flatTime+system['riseTime']/2.0+fallTime/2.0),
                    flatArea, None, None)


def makeArbitraryGrad(channel, waveform, system=None, maxGrad=None,
                      maxSlew=None):
    """
    Create an gradient event with arbitrary waveform satisfying gradient
    hardware constraints.

    See also Sequence.addBlock
    """

    validChannels = ['x', 'y', 'z']
    # #ifdef EXTERNAL_GRADS
    validChannels = validChannels + [str(x) for x in range(1, 13)]
    # #endif

    if system is None:
        system = opts()
    if maxGrad is not None:
        system['maxGrad'] = maxGrad
    if maxSlew is not None:
        system['maxSlew'] = maxSlew

    slew = (waveform[1:]-waveform[0:-1])/system['gradRasterTime']

    if np.max(np.abs(slew)) > system['maxSlew']:
        raise Exception('Slew rate violation (' +
                        str(np.max(np.abs(slew))/system['maxSlew']*100) + '%)')
    if np.max(np.abs(waveform)) > system['maxGrad']:
        raise Exception('Gradient amplitude violation (' +
                        str(np.max(np.abs(waveform))/system['maxGrad']*100) +
                        '%)')

    return Gradient('grad', channel, None, None, None, None, None, None,
                    np.arange(len(waveform))*system['gradRasterTime'],
                    waveform)


def makeAdc(numSamples, system=None, dwell=0, duration=0, delay=0,
            freqOffset=0, phaseOffset=0):
    """
    Create and ADC readout event.

    If dwell is given: Create ADC with numSamples samples with givend dwell
        time

    If duration is given: Create ADC with numSamples and specified total
        duration

    If delay is given: Create ADC with initial delay

    If system is given: Create ADC considering system properties given in
            system. For example, a dead time after sampling can be added to the
            duration

    TODO: System limits not really satisified! Apart from deadTime... (SK)

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
        dwell = duration/numSamples

    if dwell > 0:
        duration = dwell*numSamples

    return ADC('adc', numSamples, dwell, duration, delay, freqOffset,
               phaseOffset, system['adcDeadTime'])


def makeBlockPulse(flipAngle, duration=0, system=None, freqOffset=0,
                   phaseOffset=0, timeBwProduct=0, bandwidth=None,
                   maxGrad=None, maxSlew=None, sliceThickness=None):
    """
    Create an RF pulse with the given pulse shape.

    If freqOffset and phaseOffset are given, a block pulse with frequency
    offset and phase offset is created.

    If bandwith and sliceThickness are provided, an RF pulse and the
    corresponding slice select gradient is retuurned. The bandwith of the pulse
    must be given for the specified shape.

    See also Sequence.makeSincPulse, Sequence.addBlock
    """

    if system is None:
        system = opts()

    if duration == 0:
        if timeBwProduct > 0:
            duration = timeBwProduct/bandwidth
        elif bandwidth > 0:
            duration = 1/(4*bandwidth)
        else:
            raise Exception('Either bandwidth or duration must be defined.')

    BW = 1/(4*duration)
    N = np.round(duration/1e-6)
    t = (np.arange(N)+1)*system['rfRasterTime']
    signal = flipAngle/(2*np.pi)/duration*np.ones(t.shape)

    fillTime = 0
    if sliceThickness is not None:
        if maxGrad > 0:
            system['maxGrad'] = maxGrad
        if maxSlew > 0:
            system['maxSlew'] = maxSlew

        amplitude = BW/sliceThickness
        area = amplitude*duration
        gz = makeTrapezoid('z', system, flatTime=duration, flatArea=area)

        fillTime = gz.riseTime
        tFill = (np.arange(np.round(fillTime/1e-6))+1)*1e-6
        t = np.concatenate((tFill, t+tFill[-1], tFill+t[-1]+tFill[-1]))
        signal = np.concatenate((np.zeros(tFill.shape), signal,
                                 np.zeros(tFill.shape)))
    else:
        gz = None

    if fillTime < system['rfDeadTime']:
        fillTime = system['rfDeadTime'] - fillTime
        tFill = (np.arange(fillTime/1e-6)+1)*1e-6  # round to microsecond
        t = np.concatenate((tFill, t+tFill[-1]))
        signal = np.concatenate((np.zeros(tFill.shape), signal))

    if system['rfRingdownTime'] > 0:
        # Round to mu s
        tFill = np.arange(1, np.round(system['rfRingdownTime']/1e-6))*1e-6
        t = np.concatenate((t, t[-1]+tFill))
        signal = np.concatenate((signal, np.seros(tFill.shape)))

    return RFPulse('rf', signal, t, freqOffset, phaseOffset,
                   system['rfDeadTime'], system['rfRingdownTime'], gz)


def calcDuration(blocks):
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
            duration = max(duration, block.t[-1]+block.deadTime +
                           block.ringdownTime)
        elif block.type == 'grad':
            duration = max(duration, block.t[-1])
        elif block.type == 'adc':
            duration = max(duration, block.delay + block.numSamples*block.dwell
                           + block.deadTime)
        elif block.type == 'trap':
            duration = max(duration, block.riseTime + block.flatTime +
                           block.fallTime)

    return duration


def compressShape(w):
    """
    Compress gradient or pulse shape using a run-length compression
    scheme on the derivative. This strategy encodes constant and linear
    waveforms with very few samples. A structure is returned with the fields:
        numSamples: The number of samples in the uncompressed waveform
        data: contining the compressed waveform

    See also decompressShape
    """
    data = np.concatenate((w[0].flatten(), np.diff(w.flatten(), axis=0)))

    # Mask is TRUE if values change
    maskChanges = np.concatenate(([True], np.abs(np.diff(data)) > 1e-8))
    vals = data[maskChanges]  # Elements without repetitions
    # Indices of changes
    k = np.where(np.concatenate((maskChanges, [True])))[0]
    n = np.diff(k.flatten())  # Number of repetitions

    # Encode in Pulseq format
    nExtra = n-2.0
    vals2 = vals.copy()
    vals2[nExtra < 0] = np.nan
    nExtra[nExtra < 0] = np.nan
    v = np.concatenate((vals[:, None], vals2[:, None], nExtra[:, None]),
                       axis=1)
    v = v[np.isfinite(v)]
    v[np.abs(v) < 1e-10] = 0

    Shape = collections.namedtuple('shape', ['numSamples', 'data'])
    return Shape(w.size, v.flatten())


def makeDelay(delay):
    """
    Create a delay event with a given delay.

    See also Sequence.addBlock
    """

    if (not np.isfinite(delay)) or (delay <= 0):
        raise Exception('Delay (' + str(delay*1e3) + 'ms) is invalid.')

    return Delay('delay', delay)  # delay, delay, delaaaayyy :D
