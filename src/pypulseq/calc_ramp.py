import numpy as np

from pypulseq.opts import Opts


def calc_ramp(
    k0: np.ndarray,
    k_end: np.ndarray,
    max_grad: np.ndarray | None = None,
    max_points: int = 500,
    max_slew: np.ndarray | None = None,
    system: Opts | None = None,
    oversampling: bool = False,
) -> tuple[np.ndarray, bool]:
    """
    Join the points `k0` and `k_end` in three-dimensional  k-space in minimal time, observing the gradient and slew
    limits (`max_grad` and `max_slew` respectively), and the gradient strength `g0` before `k0[:, 1]` and `Gend` after
    `k_end[:, 1]`. In the context of a fixed gradient dwell time this is a discrete problem with an a priori unknown
    number of discretization steps. Therefore this method tries out the optimization with 0 steps, then 1 step, and so
    on, until  all conditions can be fulfilled, thus yielding a short connection.

    Parameters
    ----------
    k0 : numpy.ndarray
        Two preceding points in k-space. Shape is `[3, 2]`. From these points, the starting gradient will be calculated.
    k_end : numpy.ndarray
        Two following points in k-space. Shape is `[3, 2]`. From these points, the target gradient will be calculated.
    max_grad : float or array_like, default=0
        Maximum total gradient strength (Hz/m). Either a single value or one value for each coordinate, of shape `[3, 1]`.
    max_points : int, default=500
        Maximum number of k-space points to be used in connecting `k0` and `k_end`.
    max_slew : float or array_like, default=0
        Maximum total slew rate (Hz/m/s). Either a single value or one value for each coordinate, of shape `[3, 1]`.
    system : Opts, default=Opts()
        System limits.
    oversampling : bool, default=False
        Boolean flag to indicate if gradient is oversampled by a factor of 2.

    Returns
    -------
    k_out : numpy.ndarray
        Connected k-space trajectory.
    success : bool
        Boolean flag indicating if `k0` and `k_end` were successfully joined.
    """
    if system is None:
        system = Opts.default

    if max_grad is None:
        max_grad = np.zeros(0)

    if max_slew is None:
        max_slew = np.zeros(0)

    def __inside_limits(grad, slew):
        if mode == 0:
            grad2 = np.sum(np.square(grad), axis=1)
            slew2 = np.sum(np.square(slew), axis=1)
            ok = np.all(np.max(grad2) <= np.square(max_grad)) and np.all(np.max(slew2) <= np.square(max_slew))
        else:
            ok = (np.sum(np.max(np.abs(grad), axis=1) <= max_grad) == 3) and (
                np.sum(np.max(np.abs(slew), axis=1) <= max_slew) == 3
            )

        return ok

    def __joinleft0(k0, k_end, use_points, g0, g_end):
        if use_points == 0:
            g = np.stack((g0, (k_end - k0) / grad_raster, g_end)).T
            s = (g[:, 1:] - g[:, :-1]) / grad_raster

            k_out_left = np.zeros((3, 0))
            success = __inside_limits(g, s)

            return success, k_out_left

        dk = (k_end - k0) / (use_points + 1)
        kopt = k0 + dk
        g_opt = (kopt - k0) / grad_raster
        s_opt = (g_opt - g0) / grad_raster

        ok_g_opt = np.sum(np.square(g_opt)) <= np.square(max_grad)
        ok_s_opt = np.sum(np.square(s_opt)) <= np.square(max_slew)

        if ok_g_opt and ok_s_opt:
            k_left = kopt
        else:
            a = np.multiply(max_grad, grad_raster)
            b = np.multiply(max_slew, grad_raster**2)

            dkprol = g0 * grad_raster
            dkconn = dk - dkprol

            ksl = k0 + dkprol + dkconn / np.linalg.norm(dkconn) * b
            g_sl = (ksl - k0) / grad_raster
            ok_g_sl = np.sum(np.square(g_sl)) <= np.square(max_grad)

            kgl = k0 + np.multiply(dk / np.linalg.norm(dk), a)
            g_gl = (kgl - k0) / grad_raster
            s_gl = (g_gl - g0) / grad_raster
            ok_s_gl = np.sum(np.square(s_gl)) <= np.square(max_slew)

            if ok_g_sl:
                k_left = ksl
            elif ok_s_gl:
                k_left = kgl
            else:
                c = np.linalg.norm(dkprol)
                c1 = np.divide(np.square(a) - np.square(b) + np.square(c), (2 * c))
                h = np.sqrt(np.square(a) - np.square(c1))
                kglsl = k0 + np.multiply(c1, np.divide(dkprol, np.linalg.norm(dkprol)))
                projondkprol = (kgl * dkprol.T) * (dkprol / np.linalg.norm(dkprol))
                hdirection = kgl - projondkprol
                kglsl = kglsl + h * hdirection / np.linalg.norm(hdirection)
                k_left = kglsl

        success, k = __joinright0(k_left, k_end, (k_left - k0) / grad_raster, g_end, use_points - 1)
        if len(k) != 0:
            if len(k.shape) == 1:
                k = k.reshape((len(k), 1))
            if len(k_left.shape) == 1:
                k_left = k_left.reshape((len(k_left), 1))
            k_out_left = np.hstack((k_left, k))
        else:
            k_out_left = k_left

        return success, k_out_left

    def __joinleft1(k0, k_end, use_points, g0, g_end):
        if use_points == 0:
            g = np.stack((g0, (k_end - k0) / grad_raster, g_end))
            s = (g[:, 1:] - g[:, :-1]) / grad_raster

            k_out_left = np.zeros((3, 0))
            success = __inside_limits(g, s)

            return success, k_out_left

        k_left = np.zeros(3)

        dk = (k_end - k0) / (use_points + 1)
        kopt = k0 + dk
        g_opt = (kopt - k0) / grad_raster
        s_opt = (g_opt - g0) / grad_raster

        ok_g_opt = np.abs(g_opt) <= max_grad
        ok_s_opt = np.abs(s_opt) <= max_slew

        dkprol = g0 * grad_raster
        dkconn = dk - dkprol

        ksl = k0 + dkprol + np.multiply(np.sign(dkconn), max_slew) * grad_raster**2
        g_sl = (ksl - k0) / grad_raster
        ok_g_sl = np.abs(g_sl) <= max_grad

        kgl = k0 + np.multiply(np.sign(dk), max_grad) * grad_raster**2
        g_gl = (kgl - k0) / grad_raster
        s_gl = (g_gl - g0) / grad_raster
        ok_s_gl = np.abs(s_gl) <= max_slew

        for ii in range(3):
            if ok_g_opt[ii] == 1 and ok_s_opt[ii] == 1:
                k_left[ii] = kopt[ii]
            elif ok_g_sl[ii] == 1:
                k_left[ii] = ksl[ii]
            elif ok_s_gl[ii] == 1:
                k_left[ii] = kgl[ii]
            else:
                print('Unknown error')

        success, k = __joinright1(
            k0=k_left,
            k_end=k_end,
            use_points=use_points - 1,
            g0=(k_left - k0) / grad_raster,
            g_end=g_end,
        )
        if len(k) != 0:
            if len(k.shape) == 1:
                k = k.reshape((len(k), 1))
            if len(k_left.shape) == 1:
                k_left = k_left.reshape((len(k_left), 1))
            k_out_left = np.hstack((k_left, k))
        else:
            k_out_left = k_left

        return success, k_out_left

    def __joinright0(k0, k_end, use_points, g0, g_end):
        if use_points == 0:
            g = np.stack((g0, (k_end - k0) / grad_raster, g_end)).T
            s = (g[:, 1:] - g[:, :-1]) / grad_raster

            k_out_right = np.zeros((3, 0))
            success = __inside_limits(g, s)

            return success, k_out_right

        dk = (k0 - k_end) / (use_points + 1)
        kopt = k_end + dk
        g_opt = (k_end - kopt) / grad_raster
        s_opt = (g_end - g_opt) / grad_raster

        ok_g_opt = np.sum(np.square(g_opt)) <= np.square(max_grad)
        ok_s_opt = np.sum(np.square(s_opt)) <= np.square(max_slew)

        if ok_g_opt and ok_s_opt:
            k_right = kopt
        else:
            a = np.multiply(max_grad, grad_raster)
            b = np.multiply(max_slew, grad_raster**2)

            dkprol = -g_end * grad_raster
            dkconn = dk - dkprol

            ksl = k_end + dkprol + dkconn / np.linalg.norm(dkconn) * b
            g_sl = (k_end - ksl) / grad_raster
            ok_g_sl = np.sum(np.square(g_sl)) <= np.square(max_grad)

            kgl = k_end + np.multiply(dk / np.linalg.norm(dk), a)
            g_gl = (k_end - kgl) / grad_raster
            s_gl = (g_end - g_gl) / grad_raster
            ok_s_gl = np.sum(np.square(s_gl)) <= np.square(max_slew)

            if ok_g_sl:
                k_right = ksl
            elif ok_s_gl:
                k_right = kgl
            else:
                c = np.linalg.norm(dkprol)
                c1 = np.divide(np.square(a) - np.square(b) + np.square(c), (2 * c))
                h = np.sqrt(np.square(a) - np.square(c1))
                kglsl = k_end + np.multiply(c1, np.divide(dkprol, np.linalg.norm(dkprol)))
                projondkprol = (kgl * dkprol.T) * (dkprol / np.linalg.norm(dkprol))
                hdirection = kgl - projondkprol
                kglsl = kglsl + h * hdirection / np.linalg.norm(hdirection)
                k_right = kglsl

        success, k = __joinleft0(
            k0=k0,
            k_end=k_right,
            g0=g0,
            g_end=(k_end - k_right) / grad_raster,
            use_points=use_points - 1,
        )
        if len(k) != 0:
            if len(k.shape) == 1:
                k = k.reshape((len(k), 1))
            if len(k_right.shape) == 1:
                k_right = k_right.reshape((len(k_right), 1))
            k_out_right = np.hstack((k, k_right))
        else:
            k_out_right = k_right

        return success, k_out_right

    def __joinright1(k0, k_end, use_points, g0, g_end):
        if use_points == 0:
            g = np.stack((g0, (k_end - k0) / grad_raster, g_end))
            s = (g[:, 1:] - g[:, :-1]) / grad_raster

            k_out_right = np.zeros((3, 0))
            success = __inside_limits(g, s)

            return success, k_out_right

        k_right = np.zeros(3)

        dk = (k0 - k_end) / (use_points + 1)
        kopt = k_end + dk
        g_opt = (k_end - kopt) / grad_raster
        s_opt = (g_end - g_opt) / grad_raster

        ok_g_opt = np.abs(g_opt) <= max_grad
        ok_s_opt = np.abs(s_opt) <= max_slew

        dkprol = -g_end * grad_raster
        dkconn = dk - dkprol

        ksl = k_end + dkprol + np.multiply(np.sign(dkconn), max_slew) * grad_raster**2
        g_sl = (k_end - ksl) / grad_raster
        ok_g_sl = np.abs(g_sl) <= max_grad

        kgl = k_end + np.multiply(np.sign(dk), max_grad) * grad_raster
        g_gl = (k_end - kgl) / grad_raster
        s_gl = (g_end - g_gl) / grad_raster
        ok_s_gl = np.abs(s_gl) <= max_slew

        for ii in range(3):
            if ok_g_opt[ii] == 1 and ok_s_opt[ii] == 1:
                k_right[ii] = kopt[ii]
            elif ok_g_sl[ii] == 1:
                k_right[ii] = ksl[ii]
            elif ok_s_gl[ii] == 1:
                k_right[ii] = kgl[ii]
            else:
                print('Unknown error')

        success, k = __joinleft1(
            k0=k0,
            k_end=k_right,
            use_points=use_points - 1,
            g0=g0,
            g_end=(k_end - k_right) / grad_raster,
        )
        if len(k) != 0:
            if len(k.shape) == 1:
                k = k.reshape((len(k), 1))
            if len(k_right.shape) == 1:
                k_right = k_right.reshape((len(k_right), 1))
            k_out_right = np.hstack((k, k_right))
        else:
            k_out_right = k_right

        return success, k_out_right

    # =========
    # MAIN FUNCTION
    # =========
    if np.all(np.where(max_grad <= 0)):
        max_grad = [system.max_grad]
    if np.all(np.where(max_slew <= 0)):
        max_slew = [system.max_slew]

    if oversampling:
        grad_raster = 0.5 * system.grad_raster_time
    else:
        grad_raster = system.grad_raster_time

    if len(max_grad) == 1 and len(max_slew) == 1:
        mode = 0
    elif len(max_grad) == 3 and len(max_slew) == 3:
        mode = 1
    else:
        raise ValueError('Input value max grad or max slew in invalid format.')

    g0 = (k0[:, 1] - k0[:, 0]) / grad_raster
    g_end = (k_end[:, 1] - k_end[:, 0]) / grad_raster
    k0 = k0[:, 1]
    k_end = k_end[:, 0]

    success = 0
    k_out = np.zeros((3, 0))
    use_points = 0

    while success == 0 and use_points <= max_points:
        if mode == 0:
            if np.linalg.norm(g0) > max_grad or np.linalg.norm(g_end) > max_grad:
                break
            success, k_out = __joinleft0(k0=k0, k_end=k_end, g0=g0, g_end=g_end, use_points=use_points)
        else:
            if abs(g0) > abs(max_grad) or abs(g_end) > abs(max_grad):
                break
            success, k_out = __joinleft1(k0=k0, k_end=k_end, use_points=use_points, g0=g0, g_end=g_end)
        use_points += 1

    return k_out, success
