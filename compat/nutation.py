# -*- coding: utf-8 -*-

import ctypes
from novas import novaslib
from novas import _neg_err


def iau2000a(jd_high, jd_low):
    """
    Computes the forced nutation of the non-rigid Earth based on the
    IAU 2000A nutation model.

    Parameters
    ----------
    jd_high : float
        High-order part of TT Julian date.
    jd_low : float
        Low-order part of TT Julian date.

    Returns
    -------
    (dpsi, deps) : tuple of floats
        Nutation (luni-solar + planetary) in (longitude, obliquity)
        in radians.

    Notes
    -----
    .. [N1] The IAU 2000A nutation model is MHB_2000 without the
        free core nutation and without the corrections to Lieske
        precession.

    References
    ----------
    .. [R1] IERS Conventions (2003), Chapter 5.
    .. [R2] Simon et al. (1994) Astronomy and Astrophysics 282,
        663-683, esp. Sections 3.4-3.5.

    """

    if jd_high + jd_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_high+jd_low'))

    _iau2000a = novaslib.iau2000a
    _iau2000a.argtypes = (ctypes.c_double, ctypes.c_double,
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double))
    _iau2000a.restype = None

    dpsi = ctypes.c_double()
    deps = ctypes.c_double()

    _iau2000a(jd_high, jd_low, ctypes.byref(dpsi), ctypes.byref(deps))

    return (dpsi.value, deps.value)


def iau2000b(jd_high, jd_low):
    """
    Computes the forced nutation of the non-rigid Earth based on the
    IAU 2000B precession/nutation model.

    Parameters
    ----------
    jd_high : float
        High-order part of TT Julian date.
    jd_low : float
        Low-order part of TT Julian date.

    Returns
    -------
    (dpsi, deps) : tuple of floats
        Nutation (luni-solar + planetary) in (longitude, obliquity)
        in radians.

    Notes
    -----
    .. [N1] IAU 2000B reproduces the IAU 2000A model to a precision
        of 1 milliarcsecond in the interval 1995-2020.

    References
    ----------
    .. [R1] McCarthy, D. and Luzum, B. (2003). "An Abridged Model of
        the Precession & Nutation of the Celestial Pole," Celestial
        Mechanics and Dynamical Astronomy, Volume 85, Issue 1,
        Jan. 2003, p. 37. (IAU 2000B)
    .. [R2] IERS Conventions (2003), Chapter 5.

    """

    if jd_high + jd_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_high+jd_low'))

    _iau2000b = novaslib.iau2000b
    _iau2000b.argtypes = (ctypes.c_double, ctypes.c_double,
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double))
    _iau2000b.restype = None

    dpsi = ctypes.c_double()
    deps = ctypes.c_double()

    _iau2000b(jd_high, jd_low, ctypes.byref(dpsi), ctypes.byref(deps))

    return (dpsi.value, deps.value)


def nu2000k(jd_high, jd_low):
    """
    Computes the forced nutation of the non-rigid Earth: Model NU2000K.
    This model is a modified version of IAU 2000A, which has been
    truncated for speed of execution, and uses Simon et al. (1994)
    fundamental arguments throughout.NU2000K agrees with IAU 2000A at
    the 0.1 milliarcsecond level from 1700 to 2300.

    Parameters
    ----------
    jd_high : float
        High-order part of TT Julian date.
    jd_low : float
        Low-order part of TT Julian date.

    Returns
    -------
    (dpsi, deps) : tuple of floats
        Nutation (luni-solar + planetary) in (longitude, obliquity)
        in radians.

    Notes
    -----
    .. [N1] NU2000K was compared to IAU 2000A over six centuries
        (1700-2300). The average error in dpsi is 20
        microarcseconds, with 98% of the errors < 60
        microarcseconds;the average error in deleps is 8
        microarcseconds, with 100% of the errors < 60
        microarcseconds.
    .. [N2] NU2000K was developed by G. Kaplan (USNO) in March 2004.

    References
    ----------
    .. [R1] IERS Conventions (2003), Chapter 5.
    .. [R2] Simon et al. (1994) Astronomy and Astrophysics 282,
        663-683, esp. Sections 3.4-3.5.

    """

    if jd_high + jd_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_high+jd_low'))

    _nu2000k = novaslib.nu2000k
    _nu2000k.argtypes = (ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double))
    _nu2000k.restype = None

    dpsi = ctypes.c_double()
    deps = ctypes.c_double()

    _nu2000k(jd_high, jd_low, ctypes.byref(dpsi), ctypes.byref(deps))

    return (dpsi.value, deps.value)
