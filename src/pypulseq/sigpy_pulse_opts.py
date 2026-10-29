class SigpyPulseOpts:
    """Pulse configuration options for `pypulseq.make_sigpy_pulse.sigpy_n_seq`.

    Parameters
    ----------
    pulse_type : str, default='slr'
        Pulse type. Must be one of 'slr' or 'sms'.
    ptype : str, default='st'
        Pulse design method. Must be one of 'st', 'ex', 'inv', 'sat', 'se', 'fi', 'fs'.
    ftype : str, default='ls'
        Filter type. Must be one of 'ls', 'pm', 'min', 'max', 'ap'.
    d1 : float, default=0.01
        Passband ripple.
    d2 : float, default=0.01
        Stopband ripple.
    cancel_alpha_phs : bool, default=False
        Cancel alpha phase.
    n_bands : int, default=3
        Number of bands. SMS only.
    band_sep : float, default=20
        Band separation. SMS only.
    phs_0_pt : str, default='None'
        Phase 0 point. SMS only.

    Raises
    ------
    ValueError
        If `pulse_type` is not 'slr' or 'sms'.
    """

    def __init__(
        self,
        pulse_type: str = 'slr',
        ptype: str = 'st',
        ftype: str = 'ls',
        d1: float = 0.01,
        d2: float = 0.01,
        cancel_alpha_phs: bool = False,
        n_bands: int = 3,
        band_sep: int = 20,
        phs_0_pt: str = 'None',
    ):
        if pulse_type not in ('slr', 'sms'):
            raise ValueError(f"Invalid pulse_type parameter. Must be one of 'slr' or 'sms'. Passed: {pulse_type}")

        self.pulse_type = pulse_type

        if pulse_type == 'slr':
            self.ptype = ptype
            self.ftype = ftype
            self.d1 = d1
            self.d2 = d2
            self.cancel_alpha_phs = cancel_alpha_phs

        if pulse_type == 'sms':
            self.ptype = ptype
            self.ftype = ftype
            self.d1 = d1
            self.d2 = d2
            self.cancel_alpha_phs = cancel_alpha_phs
            self.n_bands = n_bands
            self.band_sep = band_sep
            self.phs_0_pt = phs_0_pt

    def __str__(self) -> str:
        s = 'Pulse options:'
        s += '\nptype: ' + str(self.ptype)
        s += '\nftype: ' + str(self.ftype)
        s += '\nd1: ' + str(self.d1)
        s += '\nd2: ' + str(self.d2)
        s += '\ncancel_alpha_phs: ' + str(self.cancel_alpha_phs)

        return s
