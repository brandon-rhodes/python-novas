# -*- coding: utf-8 -*-

import ctypes
from novas import novaslib, _check_c_errors
from novas import (_neg_err, _hour_range_err, _elev_range_err,
                   _az180_range_err, _az360_range_err, _vector_len_err,
                   _option_err, _jd_year_err, _month_range_err,
                   _day_range_err, _days_in_month)
from novas import NonConvergentError, IndeterminateError, InitializationError


class CatEntry(ctypes.Structure):
    _fields_ = [
        ('starname', ctypes.c_char*51),
        ('catalog', ctypes.c_char*4),
        ('starnumber', ctypes.c_long),
        ('ra', ctypes.c_double),
        ('dec', ctypes.c_double),
        ('promora', ctypes.c_double),
        ('promodec', ctypes.c_double),
        ('parallax', ctypes.c_double),
        ('radialvelocity', ctypes.c_double)
    ]


class Object(ctypes.Structure):
    _fields_ = [
        ('type', ctypes.c_short),
        ('number', ctypes.c_short),
        ('name', ctypes.c_char*51),
        ('star', CatEntry)
    ]


class OnSurface(ctypes.Structure):
    _fields_ = [
        ('latitude', ctypes.c_double),
        ('longitude', ctypes.c_double),
        ('height', ctypes.c_double),
        ('temperature', ctypes.c_double),
        ('pressure', ctypes.c_double)
    ]


class InSpace(ctypes.Structure):
    _fields_ = [
        ('_sc_pos_array', ctypes.c_double*3),
        ('_sc_vel_array', ctypes.c_double*3)
    ]

    def _get_sc_pos(self):
        return tuple([i for i in self._sc_pos_array])

    def _set_sc_pos(self, sc_pos_tuple):
        self._sc_pos_array = (ctypes.c_double*3)(*sc_pos_tuple)

    def _get_sc_vel(self):
        return tuple([i for i in self._sc_vel_array])

    def _set_sc_vel(self, sc_vel_tuple):
        self._sc_vel_array = (ctypes.c_double*3)(*sc_vel_tuple)

    sc_pos = property(_get_sc_pos, _set_sc_pos)
    sc_vel = property(_get_sc_vel, _set_sc_vel)


class Observer(ctypes.Structure):
    _fields_ = [
        ('where', ctypes.c_short),
        ('on_surf', OnSurface),
        ('near_earth', InSpace)
    ]


class SkyPos(ctypes.Structure):
    _fields_ = [
        ('_r_hat_array', ctypes.c_double*3),
        ('ra', ctypes.c_double),
        ('dec', ctypes.c_double),
        ('dis', ctypes.c_double),
        ('rv', ctypes.c_double)
    ]

    def _get_r_hat(self):
        return tuple([i for i in self._r_hat_array])

    def _set_r_hat(self, r_hat_tuple):
        self._r_hat_array = (ctypes.c_double*3)(*r_hat_tuple)

    r_hat = property(_get_r_hat, _set_r_hat)


class RAofCIO(ctypes.Structure):
    _fields_ = [
        ('jd_tdb', ctypes.c_double),
        ('ra_cio', ctypes.c_double)
    ]


def app_star(jd_tt, star, accuracy=0):
    """
    Computes the apparent place of a star at 'date', given its catalog
    mean place, proper motion, parallax, and radial velocity.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for apparent place.
    star : CatEntry
        Instance of CatEntry type object containing catalog data for the
        object in the ICRS.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        Apparent (right ascension in hours, declination in degrees),
        referred to true equator and equinox of date 'jd_tt'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
        C3.1', C61.
    .. [R2] Explanatory Supplement to the Astronomical Almanac (1992),
        Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _app_star = novaslib.app_star
    _app_star.argtypes = (ctypes.c_double, ctypes.POINTER(CatEntry),
                          ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double))
    _app_star.restype = ctypes.c_short
    _app_star.errcheck = _check_c_errors
    _app_star.c_errors = {
        11: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        12: (ValueError, "from C function 'make_object': 'number' out of range"),
        13: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        14: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        15: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        21: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        22: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        23: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _app_star(jd_tt, ctypes.byref(star), accuracy,
              ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def virtual_star(jd_tt, star, accuracy=0):
    """
    Computes the virtual place of a star at 'date', given its catalog
    mean place, proper motion, parallax, and radial velocity.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for virtual place.
    star : CatEntry
        Instance of CatEntry type object containing catalog data for the
        object in the ICRS.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        Virtual (right ascension in hours, declination in degrees),
        referred to the GCRS.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
        C3.1', C64.
    .. [R2] Explanatory Supplement to the Astronomical Almanac (1992),
        Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _virtual_star = novaslib.virtual_star
    _virtual_star.argtypes = (ctypes.c_double, ctypes.POINTER(CatEntry),
                              ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double))
    _virtual_star.restype = ctypes.c_short
    _virtual_star.errcheck = _check_c_errors
    _virtual_star.c_errors = {
        11: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        12: (ValueError, "from C function 'make_object': 'number' out of range"),
        13: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        14: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        15: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        21: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        22: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        23: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _virtual_star(jd_tt, ctypes.byref(star), accuracy,
                  ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def astro_star(jd_tt, star, accuracy=0):
    """
    Computes the astrometric place of a star at 'date', given its
    catalog mean place, proper motion, parallax, and radial velocity.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for astrometric place.
    star : CatEntry
        Instance of CatEntry type object containing catalog data for the
        object in the ICRS.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        Astrometric (right ascension in hours, declination in degrees),
        referred to the ICRS without light deflection or aberration.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
        C3.1', C66.
    .. [R2] Explanatory Supplement to the Astronomical Almanac (1992),
        Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _astro_star = novaslib.astro_star
    _astro_star.argtypes = (ctypes.c_double, ctypes.POINTER(CatEntry),
                            ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double))
    _astro_star.restype = ctypes.c_short
    _astro_star.errcheck = _check_c_errors
    _astro_star.c_errors = {
        11: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        12: (ValueError, "from C function 'make_object': 'number' out of range"),
        13: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        14: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        15: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        21: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        22: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        23: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _astro_star(jd_tt, ctypes.byref(star), accuracy,
                ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def app_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the apparent place of a planet or other solar system body.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for apparent place.
    ss_body : Object
        Instance of Object type object containing the body designation
        for the solar system body.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec, dis) : tuple of floats
        Apparent (right ascension in hours, declination in degrees,
        ...), referred to true equator and equinox of date, and true
        (..., ..., distance in AU) from Earth to solar system body at
        'jd_tt'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
        C3.1', C67.
    .. [R2] Explanatory Supplement to the Astronomical Almanac (1992),
        Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _app_planet = novaslib.app_planet
    _app_planet.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                            ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double))
    _app_planet.restype = ctypes.c_short
    _app_planet.errcheck = _check_c_errors
    _app_planet.c_errors = {
        1: (ValueError, "from C function 'app_planet': Invalid value of 'type' in ctypes.Structure 'ss_body'"),
        11: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        12: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        13: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dis = ctypes.c_double()

    _app_planet(jd_tt, ctypes.byref(ss_body), accuracy,
                ctypes.byref(ra), ctypes.byref(dec), ctypes.byref(dis))

    return (ra.value, dec.value, dis.value)


def virtual_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the virtual place of a planet or other solar system body.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for virtual place.
    ss_body : Object
        Instance of Object type object containing the body
        designation for the solar system body.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec, dis) : tuple of floats
        Virtual (right ascension in hours, declination in degrees,
        ...), referred to the GCRS, and true (..., ..., distance in
        AU) from Earth to solar system body.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C70.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _virtual_planet = novaslib.virtual_planet
    _virtual_planet.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                                ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double))
    _virtual_planet.restype = ctypes.c_short
    _virtual_planet.errcheck = _check_c_errors
    _virtual_planet.c_errors = {
        1: (ValueError, "from C function 'virtual_planet': Invalid value of 'type' in ctypes.Structure 'ss_body'"),
        11: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        12: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        13: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dis = ctypes.c_double()

    _virtual_planet(jd_tt, ctypes.byref(ss_body), accuracy,
                    ctypes.byref(ra), ctypes.byref(dec), ctypes.byref(dis))

    return (ra.value, dec.value, dis.value)


def astro_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the astrometric place of a planet or other solar system
    body.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for astrometric place.
    ss_body : Object
        Instance of Object type object containing the body
        designation for the solar system body.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec, dis) : tuple of floats
        Astrometric (right ascension in hours, declination in
        degrees, ...), referred to the ICRS without light deflection
        or aberration, and true (..., ..., distance in AU) from
        Earth to solar system body.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C73.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _astro_planet = novaslib.astro_planet
    _astro_planet.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                              ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double))
    _astro_planet.restype = ctypes.c_short
    _astro_planet.errcheck = _check_c_errors
    _astro_planet.c_errors = {
        1: (ValueError, "from C function 'astro_planet': Invalid value of 'type' in ctypes.Structure 'ss_body'"),
        11: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        12: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        13: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dis = ctypes.c_double()

    _astro_planet(jd_tt, ctypes.byref(ss_body), accuracy,
                  ctypes.byref(ra), ctypes.byref(dec), ctypes.byref(dis))

    return (ra.value, dec.value, dis.value)


def topo_star(jd_tt, delta_t, star, position, accuracy=0):
    """
    Computes the topocentric place of a star at 'date', given its
    catalog mean place, proper motion, parallax, and radial velocity.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for topocentric place.
    delta_t : float
        Difference TT-UT1 at 'date', in seconds of time.
    star : CatEntry
        Instance of CatEntry type object containing catalog data for
        the object in the ICRS.
    position : OnSurface
        Instance of OnSurface type object specifying the position of
        the observer.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        Topocentric (right ascension in hours, declination in
        degrees), referred to true equator and equinox of date
        'jd_tt'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C62-C63.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _topo_star = novaslib.topo_star
    _topo_star.argtypes = (ctypes.c_double, ctypes.c_double,
                           ctypes.POINTER(CatEntry), ctypes.POINTER(OnSurface),
                           ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double))
    _topo_star.restype = ctypes.c_short
    _topo_star.errcheck = _check_c_errors
    _topo_star.c_errors = {
        1: (ValueError, "from C function 'topo_star': Invalid value of 'where' in ctypes.Structure 'location'"),
        11: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        12: (ValueError, "from C function 'make_object': 'number' out of range"),
        13: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        14: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        15: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        21: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        22: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        23: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _topo_star(jd_tt, delta_t, ctypes.byref(star), ctypes.byref(position),
               accuracy, ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def local_star(jd_tt, delta_t, star, position, accuracy=0):
    """
    Computes the local place of a star at date 'date', given its
    catalog mean place, proper motion, parallax, and radial velocity.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for local place.
    delta_t : float
        Difference TT-UT1 at 'date', in seconds of time.
    star : CatEntry
        Instance of CatEntry type object containing catalog data for
        the object in the ICRS.
    position : OnSurface
        Instance of OnSurface type object specifying the position of
        the observer.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        Local (right ascension in hours, declination in degrees),
        referred to the 'local GCRS'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C65.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _local_star = novaslib.local_star
    _local_star.argtypes = (ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(CatEntry),
                            ctypes.POINTER(OnSurface), ctypes.c_short,
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double))
    _local_star.restype = ctypes.c_short
    _local_star.errcheck = _check_c_errors
    _local_star.c_errors = {
        1: (ValueError, "from C function 'local_star': Invalid value of 'where' in ctypes.Structure 'location'"),
        11: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        12: (ValueError, "from C function 'make_object': 'number' out of range"),
        13: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        14: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        15: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        21: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        22: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        23: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _local_star(jd_tt, delta_t, ctypes.byref(star), ctypes.byref(position),
                accuracy, ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def topo_planet(jd_tt, delta_t, ss_body, position, accuracy=0):
    """
    Computes the topocentric place of a solar system body.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for topocentric place.
    delta_t : float
        Difference TT-UT1 at 'date', in seconds of time.
    ss_body : Object
        Instance of Object type object containing the body designation
        for solar system body.
    position : OnSurface
        Instance of OnSurface type object  specifying the position of
        the observer.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec, dis) : tuple of floats
        Topocentric (right ascension in hours, declination in
        degrees, ...), referred to true equator and equinox of date,
        and true (..., ..., distance in AU) from Earth to solar
        system body at 'jd_tt'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C68-C69.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """



    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _topo_planet = novaslib.topo_planet
    _topo_planet.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                             ctypes.c_double, ctypes.POINTER(OnSurface),
                             ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                             ctypes.POINTER(ctypes.c_double),
                             ctypes.POINTER(ctypes.c_double))
    _topo_planet.restype = ctypes.c_short
    _topo_planet.errcheck = _check_c_errors
    _topo_planet.c_errors = {
        1: (ValueError, "from C function 'topo_planet': Invalid value of 'where' in ctypes.Structure 'location'"),
        11: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        12: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        13: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dis = ctypes.c_double()

    _topo_planet(jd_tt, ctypes.byref(ss_body), delta_t,
                 ctypes.byref(position), accuracy, ctypes.byref(ra),
                                ctypes.byref(dec), ctypes.byref(dis))

    return (ra.value, dec.value, dis.value)


def local_planet(jd_tt, delta_t, ss_body, position, accuracy=0):
    """
    Computes the local place of a solar system body.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for local place.
    delta_t : float
        Difference TT-UT1 at 'date', in seconds of time.
    ss_body : Object
        Instance of Object type object containing the body
        designation for solar system body.
    position : OnSurface
        Instance of OnSurface type object specifying the position of
        the observer.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec, dis) : tuple of floats
        Local (right ascension in hours, declination in degrees,
        ...), referred to the 'local GCRS', and true (..., ...,
        distance in AU) from Earth to solar system body at 'jd_tt'.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C72-C72.
    .. [R2] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _local_planet = novaslib.local_planet
    _local_planet.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                              ctypes.c_double, ctypes.POINTER(OnSurface),
                              ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double))
    _local_planet.restype = ctypes.c_short
    _local_planet.errcheck = _check_c_errors
    _local_planet.c_errors = {
        1: (ValueError, "from C function 'local_planet': Invalid value of 'where' in ctypes.Structure 'location'"),
        11: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        12: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        13: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dis = ctypes.c_double()

    _local_planet(jd_tt, ctypes.byref(ss_body), delta_t,
                  ctypes.byref(position), accuracy, ctypes.byref(ra),
                  ctypes.byref(dec), ctypes.byref(dis))

    return (ra.value, dec.value, dis.value)


def mean_star(jd_tt, ra, dec, accuracy=0):
    """
    Computes the ICRS position of a star, given its apparent place
    at 'date'. Proper motion, parallax, and radial velocity are assumed
    to be zero.

    Parameters
    ----------
    jd_tt : float
        TT Julian date of apparent place.
    ra : float
        Apparent right ascension in hours, referred to true equator
        and equinox of date.
    dec : float
        Apparent declination in degrees, referred to true equator
        and equinox of date.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ira, idec) : tuple of floats
        ICRS (right ascension in hours, declination in degrees).

    References
    ----------
    .. [R1] Explanatory Supplement to the Astronomical Almanac
        (1992), Chapter 3.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err.format(name='ra'))
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err.format(name='dec'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _mean_star = novaslib.mean_star
    _mean_star.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_short,
                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
    _mean_star.restype = ctypes.c_short
    _mean_star.errcheck = _check_c_errors
    _mean_star.c_errors = {
        1: (NonConvergentError, "from C function 'mean_star': Iterative process did not converge after 30 iterations."),
        2: (ValueError, "from C function 'mean_star': Length of dummy star name out of bounds."),
        3: (ValueError, "from C function 'mean_star': Length of dummy catalog name out of bounds."),
        11: (IndeterminateError, "from C function 'vector2radec': All vector components are zero; 'ra' and 'dec' are indeterminate."),
        12: (IndeterminateError, "from C function 'vector2radec': Both pos[0] and pos[1] are zero, but pos[2] is nonzero; 'ra' is indeterminate."),
        31: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        32: (ValueError, "from C function 'make_object': 'number' out of range"),
        33: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        34: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        35: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        41: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        42: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        43: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)")
        }

    ira = ctypes.c_double()
    idec = ctypes.c_double()

    _mean_star(jd_tt, ra, dec, accuracy,
               ctypes.byref(ira), ctypes.byref(idec))

    return (ira.value, idec.value)


def place(jd_tt, delta_t, cel_object, location, coord_sys, accuracy=0):
    """
    This function computes the apparent direction of a star or solar
    system body at a specified time and in a specified coordinate
    system.

    Parameters
    ----------
    jd_tt : float
        TT Julian date for place.
    delta_t : float
        Difference TT-UT1 at 'jd_tt', in seconds of time.
    cel_object : Object
        Instance of Object type object specifying the celestial
        object of interest.
    location : Observer
        Instance of Observer type object specifying the location of
        the observer.
    coord_sys : {0, 1, 2, 3}, optional
        Code specifying coordinate system of the output position.
            = 0 ... GCRS or "local GCRS"
            = 1 ... true equator and equinox of date
            = 2 ... true equator and CIO of date
            = 3 ... astrometric coordinates, i.e., without light
                deflection of aberration
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    output : SkyPos
        Instance of SkyPos type object specifying object's place on
        the sky at time 'jd_tt', with respect to the specified
        output coordinate system.

    Notes
    -----
    .. [N1] Values of 'location.where' and 'coord_sys' dictate the
        various standard kinds of place:
            location.where = 0 and system = 1: apparent place
            location.where = 1 and system = 1: topocentric place
            location.where = 0 and system = 0: virtual place
            location.where = 1 and system = 0: local place
            location.where = 0 and system = 3: astrometric place
            location.where = 1 and system = 3: topocentric
                                               astrometric place
    .. [N2] Input value of 'delta_t' is used only when
        'location.where' equals 1 or 2 (observer is on surface of
        Earth or in a near-Earth satellite).

    References
    ----------
    .. [R1] Kaplan, G. et al. (1989), Astronomical Journal 97,
        1197-1210.
    .. [R2] Klioner, S. (2003), Astronomical Journal 125, 1580-1597.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if coord_sys not in [0, 1, 2, 3]:
        raise ValueError(_option_err.format(name='coord_sys',
                                            allowed=[0, 1, 2, 3]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _place = novaslib.place
    _place.argtypes = (ctypes.c_double, ctypes.POINTER(Object), ctypes.POINTER(Observer), ctypes.c_double,
                       ctypes.c_short, ctypes.c_short, ctypes.POINTER(SkyPos))
    _place.restype = ctypes.c_short
    _place.errcheck = _check_c_errors
    _place.c_errors = {
        1: (ValueError, "from C function 'place': invalid value of 'coord_sys'"),
        2: (ValueError, "from C function 'place': invalid value of 'accuracy'"),
        3: (ValueError, "from C function 'place': Earth is the observed object, and the observer is either at the geocenter or on the Earth's surface (not permitted)"),
        # 10 + error from function 'epehemris'
        11: (ValueError, "from C function 'ephemeris': Invalid value of 'origin'."),
        12: (ValueError, "from C function 'ephemeris': Invalid value of 'type' in 'cel_obj'."),
        13: (MemoryError, "from C function 'ephemeris': Unable to allocate memory."),
        21: (ValueError, "from C function 'solarsystem': Invalid value of body or origin."),
        22: (RuntimeError, "from C function 'solarsystem': Error detected by JPL software."),
        31: (MemoryError, "from C function 'readeph': memory allocation error"),
        32: (ValueError, "from C function 'readeph': mismatch between asteroid name and number"),
        33: (ValueError, "from C function 'readeph': Julian date out of bounds"),
        34: (IOError, "from C function 'readeph': can not find Chebyshev polynomial file"),
        39: (RuntimeError, "from C function 'readeph': dummy 'readeph' function called (see note 1 in readeph0.c)"),
        # 40 + error from function 'geo_posvel'
        41: (ValueError, "from C function 'geo_posvel': invalid value of 'accuracy'"),
        # 50 + error from function 'light_time'
        51: (NonConvergentError, "from C function 'light_time': algorithm failed to converge after 10 iterations."),
        61: (ValueError, "from C function 'solarsystem': Invalid value of body or origin."),
        62: (RuntimeError, "from C function 'solarsystem': Error detected by JPL software."),
        # 70 + error from function 'grav_def'
        71: (ValueError, "from C function 'ephemeris': Invalid value of 'origin'."),
        72: (ValueError, "from C function 'ephemeris': Invalid value of 'type' in 'cel_obj'."),
        73: (MemoryError, "from C function 'ephemeris': Unable to allocate memory."),
#        81: (ValueError, "from C function 'solarsystem': Invalid value of body or origin."),
        82: (RuntimeError, "from C function 'solarsystem': Error detected by JPL software."),
#        91: (MemoryError, "from C function 'readeph': memory allocation error"),
#        92: (ValueError, "from C function 'readeph': mismatch between asteroid name and number"),
#        93: (ValueError, "from C function 'readeph': Julian date out of bounds"),
#        94: (IOError, "from C function 'readeph': can not find Chebyshev polynomial file"),
        99: (RuntimeError, "from C function 'readeph': dummy 'readeph' function called (see note 1 in readeph0.c)"),
        101: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        102: (ValueError, "from C function 'make_object': 'number' out of range"),
        103: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        104: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        105: (ValueError, "from C function 'make_object': 'name' is out of string bounds."),
        # 80 + error from function 'cio_location'
        81: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
#        91: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        92: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        93: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        94: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        95: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        96: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
        # 90 + error from function 'cio_basis'
        91: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    output = SkyPos()

    _place(jd_tt, ctypes.byref(cel_object), ctypes.byref(location), delta_t,
           coord_sys, accuracy, ctypes.byref(output))

    return output


def equ2gal(rai, deci):
    """
    Convert ICRS equatorial position to galactic position

    Parameters
    ----------
    ra : float
        ICRS right ascension in hours.
    dec : float
        ICRS declination in degrees.

    Returns
    -------
    (glon, glat) : tuple of floats
        Galactic (longitude, latitude) in degrees.

    References
    ----------
    .. [R1] Hipparcos and Tycho Catalogues, Vol. 1, Section 1.5.3.

    """

    if not 0.0 <= rai < 24.0:
        raise ValueError(_hour_range_err.format(name='rai'))
    if not -90.0 <= deci <= 90.0:
        raise ValueError(_elev_range_err.format(name='deci'))

    _equ2gal = novaslib.equ2gal
    _equ2gal.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double))
    _equ2gal.restype = None

    glon = ctypes.c_double()
    glat = ctypes.c_double()

    _equ2gal(rai, deci, ctypes.byref(glon), ctypes.byref(glat))

    return (glon.value, glat.value)


def equ2ecl(jd_tt, ra, dec, coord_sys=2, accuracy=0):
    """
    Convert equatorial position to ecliptic position

    Parameters
    ----------
    jd_tt : float
        TT Julian date of equator, equinox, and ecliptic used for
        coordinates.
    ra : float
        Right ascension in hours, referred to specified equator and
        equinox of date.
    dec : float
        Declination in degrees, referred to specified equator and
        equinox of date.
    coord_sys : {0, 1, 2}, optional
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date
            = 2 ... ICRS (ecliptic is always the mean plane)
                (default)
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (elon, elat) : tuple of floats
        Ecliptic (longitude, latitude) in degrees, referred to
        specified ecliptic and equinox of date.

    Notes
    -----
    .. [N1] To convert ICRS RA and dec to ecliptic coordinates (mean
        ecliptic and equinox of J2000.0), set 'system' = 2; the
        value of 'jd_tt' can be set to anything, since J2000.0 is
        assumed. Except for the input to this case, all input
        coordinates are dynamical.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err.format(name='ra'))
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err.format(name='dec'))
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='coord_sys',
                                            allowed=[0, 1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _equ2ecl = novaslib.equ2ecl
    _equ2ecl.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short, ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
    _equ2ecl.restype = ctypes.c_short
    _equ2ecl.errcheck = _check_c_errors
    _equ2ecl.c_errors = {
        1: (ValueError, "from C function 'equ2ecl': invalid value of 'coord_sys'")
        }

    elon = ctypes.c_double()
    elat = ctypes.c_double()

    _equ2ecl(jd_tt, coord_sys, accuracy, ra, dec,
             ctypes.byref(elon), ctypes.byref(elat))

    return elon.value, elat.value


def equ2ecl_vec(jd_tt, position, coord_sys=2, accuracy=0):
    """
    Convert an equatorial position vector to an ecliptic position
    vector.

    Parameters
    ----------
    jd_tt : float
        TT Julian date of equator, equinox, and ecliptic used for
        coordinates.
    position : tuple of floats, of length 3
        Position vector, referred to specified equator and equinox
        of date.
    coord_sys : {0, 1, 2}, optional
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date
            = 2 ... ICRS (ecliptic is always the mean plane)
                (default)
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, referred to specified ecliptic and equinox
        of date.

    Notes
    -----
    .. [N1] To convert an ICRS vector to an ecliptic vector (mean
        ecliptic and equinox of J2000.0 only), set 'system' = 2; the
        value of 'jd_tt' can be set to anything, since J2000.0 is
        assumed. Except for the input to this case, all vectors are
        assumed to be with respect to a dynamical system.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='coord_sys',
                                            allowed=[0, 1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _equ2ecl_vec = novaslib.equ2ecl_vec
    _equ2ecl_vec.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short, ctypes.POINTER(ctypes.c_double*3),
                             ctypes.POINTER(ctypes.c_double*3))
    _equ2ecl_vec.restype = ctypes.c_short
    _equ2ecl_vec.errcheck = _check_c_errors
    _equ2ecl_vec.c_errors = {
        1: (ValueError, "from C function 'equ2ecl_vec': invalid value of 'coord_sys'")
        }

    pos2 = (ctypes.c_double*3)()

    _equ2ecl_vec(jd_tt, coord_sys, accuracy,
                 ctypes.byref((ctypes.c_double*3)(*position)),
                 ctypes.byref(pos2))

    return tuple([i for i in pos2])


def ecl2equ_vec(jd_tt, position, coord_sys=2, accuracy=0):
    """
    Converts an ecliptic position vector to an equatorial position
    vector.

    Parameters
    ----------
    jd_tt : float
        TT Julian date of equator, equinox, and ecliptic used for
        coordinates.
    position : tuple of floats, of length 3
        Position vector, referred to specified ecliptic and equinox
        of date. If 'system' = 2, 'position' must be on mean
        ecliptic and equinox of J2000.0; see Note [N1]_ below.
    coord_sys : {0, 1, 2}, optional
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date
            = 2 ... ICRS (ecliptic is always the mean plane)
                (default)
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    position2 : tuple of floats, of length 3
        Position vector, referred to specified ecliptic and equinox
        of date.

    Notes
    -----
    .. [N1] To convert an ecliptic vector (mean ecliptic and equinox
        of J2000.0 only) to an ICRS vector, set 'system' = 2; the
        value of 'jd_tt' can be set to anything, since J2000.0 is
        assumed. Except for the output from this case, all vectors
        are assumed to be with respect to a dynamical system.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='coord_sys',
                                            allowed=[0, 1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _ecl2equ_vec = novaslib.ecl2equ_vec
    _ecl2equ_vec.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short,
                             ctypes.POINTER(ctypes.c_double*3),
                             ctypes.POINTER(ctypes.c_double*3))
    _ecl2equ_vec.restype = ctypes.c_short
    _ecl2equ_vec.errcheck = _check_c_errors
    _ecl2equ_vec.c_errors = {
        1: (ValueError, "from C function 'ecl2equ_vec': invalid value of 'coord_sys'")
        }

    pos2 = (ctypes.c_double*3)()

    _ecl2equ_vec(jd_tt, coord_sys, accuracy,
                 ctypes.byref((ctypes.c_double*3)(*position)),
                 ctypes.byref(pos2))

    return tuple([i for i in pos2])


def equ2hor(jd_ut1, delta_t, xp, yp, location, ra, dec, ref_option=0,
            accuracy=0):
    """
    This function transforms topocentric right ascension and declination
    to zenith distance and azimuth. It uses a method that properly
    accounts for polar motion, which is significant at the sub-arcsecond
    level. This function can also adjust coordinates for atmospheric
    refraction.

    Parameters
    ----------
    jd_ut1 : float
        UT1 Julian date.
    delta_t : float
        Difference TT-UT1 at 'jd_ut1', in seconds.
    xp : float
        Conventionally-defined x coordinate of celestial
        intermediate pole with respect to ITRS reference pole, in
        arcseconds.
    yp : float
        Conventionally-defined y coordinate of celestial
        intermediate pole with respect to ITRS reference pole, in
        arcseconds.
    location : OnSurface
        Instance of OnSurface type object containing observer's
        location.
    ra : float
        Topocentric right ascension of object of interest, in hours,
        referred to true equator and equinox of date.
    dec : float
        Topocentric declination of object of interest, in degrees,
        referred to true equator and equinox of date.
    ref_option : {0, 1, 2}, optional
            = 0 ... no refraction (default)
            = 1 ... include refraction, using 'standard' atmospheric
                conditions.
            = 2 ... include refraction, using atmospheric parameters
                input in the 'location' object.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (zd, az) : tuple of floats
        Topocentric (zenith distance, azimuth) in degrees. Zenith
        distance is affected by refraction if 'ref_option', is
        non-zero. Azimuth is measured east from north.
    (rar, decr) : tuple of floats
        Topocentric (right ascension in hours, declination in
        degrees) of object of interest, referred to true equator and
        equinox of date, affected by refraction if 'ref_option' is
        non-zero.

    Notes
    -----
    .. [N1] 'xp' and 'yp' can be set to zero if sub-arcsecond
        accuracy is not needed. 'ra' and 'dec' can be obtained from
        functions 'topo_star' or 'topo_planet'.
    .. [N2] The directions 'zd'= 0 (zenith) and 'az'= 0 (North) are
        here considered fixed in the terrestrial system.
        Specifically, the zenith is along the geodetic normal, and
        North is toward the ITRS pole.
    .. [N3] If 'ref_option'= 0, then 'rar'='ra' and 'decr'='dec'.

    References
    ----------
    .. [R1] Kaplan, G. (2008). USNO/AA Technical Note of 28 Apr
        2008, "Refraction as a Vector."

    """

    if jd_ut1 < 0.0:
        raise ValueError(_neg_err.format(name='jd_ut1'))
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err.format(name='ra'))
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err.format(name='dec'))
    if ref_option not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='ref_option',
                                            allowed=[0, 1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _equ2hor = novaslib.equ2hor
    _equ2hor.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_short, ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(OnSurface), ctypes.c_double, ctypes.c_double, ctypes.c_short,
                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
    _equ2hor.restype = None

    zd = ctypes.c_double()
    az = ctypes.c_double()
    rar = ctypes.c_double()
    decr = ctypes.c_double()

    _equ2hor(jd_ut1, delta_t, accuracy, xp, yp, ctypes.byref(location), ra, dec,
             ref_option, ctypes.byref(zd), ctypes.byref(az), ctypes.byref(rar),
             ctypes.byref(decr))

    return (zd.value, az.value), (rar.value, decr.value)


def gcrs2equ(jd_tt, rag, decg, coord_sys=1, accuracy=0):
    """
    Convert GCRS right ascension and declination to coordinates with
    respect to the equator of date (mean or true). For coordinates with
    respect to the true equator of date, the origin of right ascension
    can be either the true equinox or the celestial intermediate origin
    (CIO).

    Parameters
    ----------
    jd_tt : float
        TT Julian date of equator to be used for output coordinates.
    rag : float
        GCRS right ascension in hours.
    decg : float
        GCRS declination in degrees.
    coord_sys : {0, 1, 2}, optional
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date (default)
            = 2 ... true equator and CIO of date
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (ra, dec) : tuple of floats
        (Right ascension in hours, declination in degrees), referred
        to specified equator (and right ascension origin, for 'ra')
        of date.

    Notes
    -----
    .. [N1] Set input value of 'accuracy' equal to any short int if
        'coord_sys' equals 0 or 1. It is not used in these cases.
    .. [N2] This function only supports the CIO-based method.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if not 0.0 <= rag < 24.0:
        raise ValueError(_hour_range_err.format(name='rag'))
    if not -90.0 <= decg <= 90.0:
        raise ValueError(_elev_range_err.format(name='decg'))
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='coord_sys',
                                            allowed=[0, 1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _gcrs2equ = novaslib.gcrs2equ
    _gcrs2equ.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short,
                          ctypes.c_double, ctypes.c_double,
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double))
    _gcrs2equ.restype = ctypes.c_short
    _gcrs2equ.errcheck = _check_c_errors
    _gcrs2equ.c_errors = {
        -1: (IndeterminateError, "from C function 'vector2radec': All vector components are zero; 'ra' and 'dec' are indeterminate."),
        -2: (IndeterminateError, "from C function 'vector2radec': Both pos[0] and pos[1] are zero, but pos[2] is nonzero; 'ra' is indeterminate."),
        11: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        21: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        22: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        23: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        24: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        25: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        26: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
#        21: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _gcrs2equ(jd_tt, coord_sys, accuracy, rag, decg,
              ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def sidereal_time(jd_high, jd_low, delta_t, gst_type=1, method=0, accuracy=0):
    """
    Computes the Greenwich sidereal time, either mean or apparent, at
    Julian date 'jd_high' + 'jd_low'.

    Parameters
    ----------
    jd_high : float
        High-order part of UT1 Julian date.
    jd_low : float
        Low-order part of UT1 Julian date.
    delta_t  float
        Difference TT-UT1 at 'jd_high' + 'jd_low', in seconds of
        time.
    gst_type : {0, 1}, optional
            = 0 ... compute Greenwich mean sidereal time
            = 1 ... compute Greenwich apparent sidereal time
                (default)
    method : {0, 1}, optional
        Selection for method.
            = 0 ... CIO-based method (default)
            = 1 ... equinox-based method
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output time.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    gst : float
        Greenwich (mean or apparent) sidereal time, in hours.

    Notes
    -----
    .. [N1] The Julian date may be split at any point, but for
        highest precision, set 'jd_high' to be the integral part of
        the Julian date, and set 'jd_low' to be the fractional part.

    References
    ----------
    .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if jd_high + jd_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_high+jd_low'))
    if gst_type not in [0, 1]:
        raise ValueError(_option_err.format(name='gst_type', allowed=[0, 1]))
    if method not in [0, 1]:
        raise ValueError(_option_err.format(name='method', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _sidereal_time = novaslib.sidereal_time
    _sidereal_time.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_short, ctypes.c_short,
                               ctypes.c_short, ctypes.POINTER(ctypes.c_double))
    _sidereal_time.restype = ctypes.c_short
    _sidereal_time.errcheck = _check_c_errors
    _sidereal_time.c_errors = {
        1: (ValueError, "from C function 'sidereal_time': invalid value of 'accuracy'"),
        2: (ValueError, "from C function 'sidereal_time': invalid value of 'method'"),
        11: (ValueError, "from C function 'cio_ra': invalid value of 'accuracy'."),
        21: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        31: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        32: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        33: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        34: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        35: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        36: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
#        31: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    gst = ctypes.c_double()

    _sidereal_time(jd_high, jd_low, delta_t, gst_type, method, accuracy,
                   ctypes.byref(gst))

    return gst.value


def era(jd_ut1_high, jd_ut1_low=0.0):
    """
    Compute the Earth Rotation Angle (theta) for a given UT1 Julian
    date. The expression used is taken from the note to IAU Resolution
    B1.8 of 2000[R1]_.

    Parameters
    ----------
    jd_ut1_high : float
        High-order part of UT1 Julian date.
    jd_ut1_low : float (optional)
        Low-order part of UT1 Julian date.

    Returns
    -------
    theta : float
        The Earth Rotation Angle in degrees.

    Notes
    -----
    .. [N1] The algorithm used here is equivalent to the canonical
        theta = 0.7790572732640 + 1.00273781191135448 * t,
        where t is the time in days from J2000 (t = jd_ut1_high +
        jd_ut1_low - T0), but it avoids many two-PI 'wraps' that
        decrease precision (adopted from SOFA Fortran routine
        iau_era00; see also expression at top of page 35 of IERS
        Conventions (1996)).

    References
    ----------
    .. [R1] IAU Resolution B1.8, adopted at the 2000 IAU General
        Assembly, Manchester, UK.
    .. [R2] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_ut1_high+jd_ut1_low'))

    _era = novaslib.era
    _era.argtypes = (ctypes.c_double, ctypes.c_double)
    _era.restype = ctypes.c_double

    theta = _era(jd_ut1_high, jd_ut1_low)

    return theta


def ter2cel(jd_ut1_high, jd_ut1_low, delta_t, xp, yp, vec1, method=0, option=0,
            accuracy=0):
    """
    This function rotates a vector from the terrestrial to the
    celestial system. Specifically, it transforms a vector in the
    ITRS (rotating earth-fixed system) to the GCRS (a local space-
    fixed system) by applying rotations for polar motion, Earth
    rotation, nutation, precession, and the dynamical-to-GCRS
    frame tie.

    Parameters
    ----------
    jd_ut1_high : float
        High-order part of UT1 Julian date.
    jd_ut1_low : float
        Low-order part of UT1 Julian date.
    delta_t : float
        Value of Delta T (=TT-UT1) at the input UT1 Julian date.
    xp : float
        Conventionally-defined X coordinate of celestial
        intermediate pole with respect to ITRS pole, in arcseconds.
    yp : float
        Conventionally-defined Y coordinate of celestial
        intermediate pole with respect to ITRS pole, in arcseconds.
    vec1 : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to ITRS axes (terrestrial system) in
        the normal case where 'option' = 0.
    method : {0, 1}, optional
        Selection for method.
            = 0 ... CIO-based method (default)
            = 1 ... equinox-based method
    option : {0, 1}, optional
        Selection for reference axes
            = 0 ... The output vector is referred to GCRS axes.
                (default)
            = 1 ... The output vector is produced with respect to
                the equator and equinox of date. See Note [N2]_
                below.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    vec2 : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to GCRS axes (celestial system) or
        with respect to the equator and equinox of date, depending
        on 'option'.

    Notes
    -----
    .. [N1] 'xp' = 'yp' = 0 means no polar motion transformation.
    .. [N2] The 'option' flag only works for the equinox-based
        method.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C53-C54.
    .. [R2] Kaplan, G. H. (2003), 'Another Look at Non-Rotating
        Origins', Proceedings of IAU XXV Joint Discussion 16.

    """

    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_ut1_high+jd_ut1_low'))
    if len(vec1) != 3:
        raise ValueError(_vector_len_err.format(name='vec1'))
    if method not in [0, 1]:
        raise ValueError(_option_err.format(name='method', allowed=[0, 1]))
    if option not in [0, 1]:
        raise ValueError(_option_err.format(name='option', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _ter2cel = novaslib.ter2cel
    _ter2cel.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double,
                         ctypes.c_short, ctypes.c_short, ctypes.c_short,
                         ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(ctypes.c_double*3),
                         ctypes.POINTER(ctypes.c_double*3))
    _ter2cel.restype = ctypes.c_short
    _ter2cel.errcheck = _check_c_errors
    _ter2cel.c_errors = {
        1: (ValueError, "from C function 'ter2cel': invalid value of 'accuracy'"),
        2: (ValueError, "from C function 'ter2cel': invalid value of 'method'"),
        11: (ValueError, "from C function 'cio_ra': invalid value of 'accuracy'."),
        21: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        31: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        32: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        33: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        34: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        35: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        36: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
#        31: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    vec2 = (ctypes.c_double*3)()

    _ter2cel(jd_ut1_high, jd_ut1_low, delta_t, method, accuracy, option,
             xp, yp, ctypes.byref((ctypes.c_double*3)(*vec1)),
             ctypes.byref(vec2))

    return tuple([i for i in vec2])


def cel2ter(jd_ut1_high, jd_ut1_low, delta_t, xp, yp, vec1, method=0, option=0,
            accuracy=0):
    """
    Rotates a vector from the celestial to the terrestrial system.
    Specifically, transforms a vector in the GCRS (a local space-fixed
    system) to the ITRS (a rotating earth-fixed system) by applying
    rotations for the GCRS-to-dynamical frame tie, precession, nutation,
    Earth rotation, and polar motion.

    Parameters
    ----------
    jd_ut1_high : float
        High-order part of UT1 Julian date.
    jd_ut1_low : float
        Low-order part of UT1 Julian date.
    delta_t : float
        Value of Delta T (= TT - UT1) at the input UT1 Julian date.
    xp : float
        Conventionally-defined X coordinate of celestial
        intermediate pole with respect to ITRS pole, in arcseconds.
    yp : float
        Conventionally-defined Y coordinate of celestial
        intermediate pole with respect to ITRS pole, in arcseconds.
    vec1 : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to GCRS axes (celestial system) or
        with respect to the equator and equinox of date, depending
        on 'option'.
    method : {0, 1}, optional
        Selection for method
            = 0 ... CIO-based method (default)
            = 1 ... equinox-based method
    option : {0, 1}, optional
        Selection for reference axes
            = 0 ... The output vector is referred to GCRS axes.
                (default)
            = 1 ... The output vector is produced with respect to
                the equator and equinox of date. See Note 2 below.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    vec2 : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular coordinates,
        referred to ITRS axes (terrestrial system) in the normal case
        where 'option' = 0.

    Notes
    -----
    .. [N1] 'xp' = 'yp' = 0 means no polar motion transformation.
    .. [N2] The 'option' flag only works for the equinox-based method.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C54.
    .. [R2] Kaplan, G. H. (2003), 'Another Look at Non-Rotating
        Origins', Proceedings of IAU XXV Joint Discussion 16.

    """



    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_ut1_high+jd_ut1_low'))
    if len(vec1) != 3:
        raise ValueError(_vector_len_err.format(name='vec1'))
    if method not in [0, 1]:
        raise ValueError(_option_err.format(name='method', allowed=[0, 1]))
    if option not in [0, 1]:
        raise ValueError(_option_err.format(name='option', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _cel2ter = novaslib.cel2ter
    _cel2ter.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double,
                         ctypes.c_short, ctypes.c_short, ctypes.c_short,
                         ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(ctypes.c_double*3),
                         ctypes.POINTER(ctypes.c_double*3))
    _cel2ter.restype = ctypes.c_short
    _cel2ter.errcheck = _check_c_errors
    _cel2ter.c_errors = {
        1: (ValueError, "from C function 'cel2ter': invalid value of 'accuracy'"),
        2: (ValueError, "from C function 'cel2ter': invalid value of 'method'"),
        11: (ValueError, "from C function 'cio_ra': invalid value of 'accuracy'."),
        21: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        31: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        32: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        33: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        34: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        35: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        36: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
#        31: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    vec2 = (ctypes.c_double*3)()

    _cel2ter(jd_ut1_high, jd_ut1_low, delta_t, method, accuracy,
             option, xp, yp, ctypes.byref((ctypes.c_double*3)(*vec1)),
             ctypes.byref(vec2))

    return tuple([i for i in vec2])


def spin(angle, pos1):
    """
    This function transforms a vector from one coordinate system to
    another with same origin and axes rotated about the z-axis.

    Parameters
    ----------
    angle : float
        Angle of coordinate system rotation, positive
        counterclockwise when viewed from +z, in degrees.
    position : tuple of floats, of length 3
        Position vector.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector expressed in new coordinate system rotated
        about z by 'angle'.

    """

    if not 0.0 <= angle < 360.0:
        raise ValueError(_az360_range_err.format(name='angle'))
    if len(pos1) != 3:
        raise ValueError(_vector_len_err.format(name='pos1'))

    _spin = novaslib.spin
    _spin.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double*3),
                      ctypes.POINTER(ctypes.c_double*3))
    _spin.restype = None

    pos2 = (ctypes.c_double*3)()

    _spin(angle, ctypes.byref((ctypes.c_double*3)(*pos1)), ctypes.byref(pos2))

    return tuple([i for i in pos2])


def wobble(tjd, x, y, position, direction=0):
    """
    Corrects a vector in the ICRS (rotating Earth-fixed system) for
    polar motion, and also corrects the longitude origin (by a tiny
    amount) to the Terrestrial Intermediate Origin (TIO). The ICRS
    vector is thereby transformed to the system based on the true
    (rotational) equator and TIO. Since the true equator is the plane
    orthogonal to the direction of the Celestial Intermediate Pole
    (CIP), the components of the output vector are referred to Z and X
    axes toward the CIP and TIO, respectively.

    Parameters
    ----------
    tjd : float
        TT or UT1 Julian date.
    xp : float
        Conventionally-defined X coordinate of Celestial
        Intermediate Pole with respect to ICRS pole, in arcseconds.
    yp : float
        Conventionally defined Y coordinate of Celestial
        Intermediate Pole with respect to ICRS pole, in arcseconds.
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to ITRF axes.
    direction: short int
        Flag determining 'direction' of transformation;
            = 0 ... transformation applied, ICRS to terrestrial
                intermediate system (default)
            != 0 ... inverse transformation applied, terrestrial
                intermediate system to ICRS

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to true equator and TIO.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C23, C110, C111, C113, C114.
    .. [R2] Lambert & Bizouard (2002), Astronomy and Astrophysics
        394, 317-321.

    """

    if tjd < 0.0:
        raise ValueError(_neg_err.format(name='tjd'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if direction not in [0, 1]:
        raise ValueError(_option_err.format(name='direction', allowed=[0, 1]))

    _wobble = novaslib.wobble
    _wobble.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_double,
                        ctypes.c_double, ctypes.POINTER(ctypes.c_double*3),
                        ctypes.POINTER(ctypes.c_double*3))
    _wobble.restype = None

    pos2 = (ctypes.c_double*3)()

    _wobble(tjd, direction, x, y, ctypes.byref((ctypes.c_double*3)(*position)),
            ctypes.byref(pos2))

    return tuple([i for i in pos2])


def terra(location, time):
    """
    Computes the position and velocity vectors of a terrestrial observer
    with respect to the center of the Earth.

    Parameters
    ----------
    location : OnSurface
        Instance of OnSurface type object containing observer's
        location.
    time : float
        Local apparent sidereal time at reference meridian in hours.

    Returns
    -------
    (position, velocity) : tuple of tuple of floats, of length 3
        Position and velocity vector of observer with respect to
        center of Earth in equatorial rectangular coordinates,
        referred to true equator and equinox of date. Components in
        AU and AU/day.

    Notes
    -----
    .. [N1] If reference meridian is Greenwich and st=0, 'pos' is
        effectively referred to equator and Greenwich.
    .. [N2] This function ignores polar motion, unless the
        observer's longitude and latitude have been corrected for
        it, and variation in the length of day (angular velocity of
        earth).
    .. [N3] The true equator and equinox of date do not form an
        inertial system. Therefore, with respect to an inertial
        system, the very small velocity component (several
        meters/day) due to the precession and nutation of the
        Earth's axis is not accounted for here.

    """

    if not 0.0 <= time < 24.0:
        raise ValueError(_hour_range_err.format(name='time'))

    _terra = novaslib.terra
    _terra.argtypes = (ctypes.POINTER(OnSurface), ctypes.c_double,
                       ctypes.POINTER(ctypes.c_double*3),
                       ctypes.POINTER(ctypes.c_double*3))
    _terra.restype = None

    pos = (ctypes.c_double*3)()
    vel = (ctypes.c_double*3)()

    _terra(ctypes.byref(location), time, ctypes.byref(pos), ctypes.byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])


def e_tilt(jd_tdb, accuracy=0):
    """
    Computes quantities related to the orientation of the Earth's
    rotation axis at the given Julian day.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        quantities.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (mobl, tobl, ee, dpsi, deps) : tuple of five floats
        mobl = mean obliquity of the ecliptic in degrees.
        tobl = true obliquity of the ecliptic in degrees.
        ee   = equation of the equinoxes in seconds of time.
        dpsi = nutation in longitude in arcseconds.
        deps = nutation in obliquity in arcseconds.

    Notes
    -----
    .. [N1] Values of the celestial pole offsets 'PSI_COR' and
        'EPS_COR' are set using function 'cel_pole', if desired. See
        the docstring for 'cel_pole' for details.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _e_tilt = novaslib.e_tilt
    _e_tilt.argtypes = (ctypes.c_double, ctypes.c_short,
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double))
    _e_tilt.restype = None

    mobl = ctypes.c_double()
    tobl = ctypes.c_double()
    ee = ctypes.c_double()
    dpsi = ctypes.c_double()
    deps = ctypes.c_double()

    _e_tilt(jd_tdb, accuracy, ctypes.byref(mobl), ctypes.byref(tobl),
            ctypes.byref(ee), ctypes.byref(dpsi), ctypes.byref(deps))

    return mobl.value, tobl.value, ee.value, dpsi.value, deps.value


def cel_pole(tjd, type, dpole1, dpole2):
    """
    This function allows for the specification of celestial pole
    offsets for high-precision applications. Each set of offsets is
    a correction to the modeled position of the pole for a specific
    date, derived from observations and published by the IERS.

    Parameters
    ----------
    tjd : float
        TDB or TT Julian day for pole offsets.
    type : {1, 2}
        Type of pole offset.
            = 1 for corrections to angular coordinates of modeled
                pole referred to mean ecliptic of date, that is,
                delta-delta-psi and delta-delta-epsilon.
            = 2 for corrections to components of modeled pole unit
                vector referred to GCRS axes, that is, dx and dy.
    dpole1 : float
        Value of celestial pole offset in first coordinate,
        (delta-delta-psi or dx) in milliarcseconds.
    dpole2 : float
        Value of celestial pole offset in second coordinate,
        (delta-delta-epsilon or dy) in milliarcseconds.

    Notes
    -----
    .. [N1] This function sets the values of global variables
        'PSI_COR' and 'EPS_COR' declared at the top of file
        'novas.c'. These global variables are used only in NOVAS
        function 'e_tilt'.
    .. [N2] This function, if used, should be called before any
        other NOVAS functions for a given date. Values of the pole
        offsets specified via a call to this function will be used
        until explicitly changed.
    .. [N3] 'tjd' is used only for 'type' = 2, to transform dx and
        dy to the equivalent delta-delta-psi and delta-delta-epsilon
        values.
    .. [N4] For 'type' = 2, dx and dy are unit vector component
        corrections, but are expressed in milliarcseconds simply by
        multiplying by 206264806, the number of milliarcseconds in
        one radian.

    References
    ----------
    .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
    .. [R2] Kaplan, G. (2003), USNO/AA Technical Note 2003-03.

    """

    if tjd < 0.0:
        raise ValueError(_neg_err.format(name='tjd'))
    if type not in [1, 2]:
        raise ValueError(_option_err.format(name='type', allowed=[1, 2]))

    _cel_pole = novaslib.cel_pole
    _cel_pole.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_double, ctypes.c_double)
    _cel_pole.restype = ctypes.c_short
    _cel_pole.errcheck = _check_c_errors
    _cel_pole.c_errors = {
        1: (ValueError, "from C function 'cel_pole': Invalid value of 'type'")
        }

    _cel_pole(tjd, type, dpole1, dpole2)

    return None


def ee_ct(jd_tt_high, jd_tt_low, accuracy=0):
    """
    To compute the \"complementary terms\" of the equation of the
    equinoxes.

    Parameters
    ----------
    jd_tt_high : float
        High-order (integer) part of TT Julian day.
    jd_tt_low : float
        Low-order (fractional) part of TT Julian day.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the terms.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    comp_terms : float
        Complementary terms, in radians.

    Notes
    -----
    .. [N1] The series used in this function was derived from the
        first reference.  This same series was also adopted for use
        in the IAU's Standards of Fundamental Astronomy (SOFA)
        software (i.e., subroutine eect00.for and function
        eect00.c).
    .. [N2] The low-accuracy series used in this function is a
        simple implementation derived from the first reference, in
        which terms smaller than 2 microarcseconds have been
        omitted.

    References
    ----------
    .. [R1] Capitaine, N., Wallace, P.T., and McCarthy, D.D. (2003).
        Astron. & Astrophys. 406, p. 1135-1149. Table 3.
    .. [R2] IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
        (Table 5.2e presented in the printed publication is a
        truncated series. The full series, which is used in NOVAS,
        is available on the IERS Conventions Center website in file
        tab5.2e.txt.) ftp://tai.bipm.org/iers/conv2010/chapter5/

    """

    if jd_tt_high + jd_tt_low < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt_high+jd_tt_low'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _ee_ct = novaslib.ee_ct
    _ee_ct.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_short)
    _ee_ct.restype = ctypes.c_double

    comp_terms = _ee_ct(jd_tt_high, jd_tt_low, accuracy)

    return comp_terms


def frame_tie(position, direction=0):
    """
    Transform a vector from the dynamical reference system to the
    International Celestial Reference System (ICRS), or vice versa. The
    dynamical reference system is based on the dynamical mean equator
    and equinox of J2000.0. The ICRS is based on the space-fixed ICRS
    axes defined by the radio catalog positions of several hundred
    extragalactic objects.

    Parameters
    ----------
    position : tuple of floats, of length 3
        Position vector in equatorial rectangular coordinates.
    direction : {0, -1}, optional
        = 0 ... ICRS to dynamical (default)
        = -1 ... dynamical to ICRS

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector in equatorial rectangular coordinates.

    Notes
    -----
    .. [N1] For geocentric coordinates, the same transformation is
        used between the dynamical reference system and the GCRS.

    References
    ----------
    .. [R1] Hilton, J. and Hohenkerk, C. (2004), Astronomy and
        Astrophysics 413, 765-770, eq. (6) and (8).
    .. [R2] IERS (2003) Conventions, Chapter 5.

    """

    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if direction not in [0, -1]:
        raise ValueError(_option_err.format(name='direction', allowed=[0, -1]))

    _frame_tie = novaslib.frame_tie
    _frame_tie.argtypes = (ctypes.POINTER(ctypes.c_double*3), ctypes.c_short,
                           ctypes.POINTER(ctypes.c_double*3))
    _frame_tie.restype = None

    pos2 = (ctypes.c_double*3)()

    _frame_tie(ctypes.byref((ctypes.c_double*3)(*position)), direction,
               ctypes.byref(pos2))

    return tuple([i for i in pos2])


def proper_motion(jd_tdb1, position, velocity, jd_tdb2):
    """
    Applies proper motion, including foreshortening effects, to a star's
    position.

    Parameters
    ----------
    jd_tdb1 : float
        TDB Julian date of first epoch.
    position : tuple of floats, of length 3
        Position vector at first epoch.
    velocity : tuple of floats, of length 3
        Velocity vector at first epoch.
    jd_tdb2 : float
        TDB Julian date of second epoch.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector at second epoch.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C16.

    """

    if jd_tdb1 < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb1'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if len(velocity) != 3:
        raise ValueError(_vector_len_err.format(name='velocity'))
    if jd_tdb2 < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb2'))

    _proper_motion = novaslib.proper_motion
    _proper_motion.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double*3),
                               ctypes.POINTER(ctypes.c_double*3), ctypes.c_double,
                               ctypes.POINTER(ctypes.c_double*3))
    _proper_motion.restype = None

    pos2 = (ctypes.c_double*3)()

    _proper_motion(jd_tdb1, ctypes.byref((ctypes.c_double*3)(*position)),
                   ctypes.byref((ctypes.c_double*3)(*velocity)), jd_tdb2,
                   ctypes.byref(pos2))

    return tuple([i for i in pos2])


def bary2obs(pos, pos_obs):
    """
    Transform the origin from the solar system barycenter to the
    observer (or the geocenter); i.e., correct for parallax
    (annual+geocentric or just annual).

    Parameters
    ----------
    pos : tuple of floats, of length 3
        Position vector, referred to origin at solar system
        barycenter, components in AU.
    pos_obs : tuple of floats, of length 3
        Position vector of observer (or the geocenter), with respect
        to origin at solar system barycenter, components in AU.

    Returns
    -------
    (position, lighttime) : tuple of tuple and float
        Position vector, referred to origin at center of mass of the
        Earth, components in AU, and lighttime of object from Earth
        in days.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C17, C109.

    """

    if len(pos) != 3:
        raise ValueError(_vector_len_err.format(name='pos'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))

    _bary2obs = novaslib.bary2obs
    _bary2obs.argtypes = (ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double))
    _bary2obs.restype = None

    pos2 = (ctypes.c_double*3)()
    lighttime = ctypes.c_double()

    _bary2obs(ctypes.byref((ctypes.c_double*3)(*pos)),
              ctypes.byref((ctypes.c_double*3)(*pos_obs)),
              ctypes.byref(pos2), ctypes.byref(lighttime))

    return tuple([i for i in pos2]), lighttime.value


def geo_posvel(jd_tt, delta_t, observer, accuracy=0):
    """
    Compute the geocentric position and velocity of an observer on the
    surface of the earth or on a near-earth spacecraft. The final
    vectors are expressed in the GCRS.

    Parameters
    ----------
    jd_tt : float
        TT Julian date.
    delta_t : float
        Value of Delta T (=TT-UT1) at 'date'.
    observer : Observer
        Instance of Observer type object specifying the location of
        the observer.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position
        and velocity.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (position, velocity) : tuple of tuple of floats, of length 3
        Position vector of observer, with respect to origin at
        geocenter, referred to GCRS axes, components in AU and
        AU/day, respectively.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _geo_posvel = novaslib.geo_posvel
    _geo_posvel.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_short,
                            ctypes.POINTER(Observer),
                            ctypes.POINTER(ctypes.c_double*3),
                            ctypes.POINTER(ctypes.c_double*3))
    _geo_posvel.restype = ctypes.c_short
    _geo_posvel.errcheck = _check_c_errors
    _geo_posvel.c_errors ={
        1: (ValueError, "from C function 'geo_posvel': invalid value of 'accuracy'")
        }

    pos = (ctypes.c_double*3)()
    vel = (ctypes.c_double*3)()

    _geo_posvel(jd_tt, delta_t, accuracy, ctypes.byref(observer),
                ctypes.byref(pos), ctypes.byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])


def light_time(jd_tdb, ss_object, pos_obs, estimate=0.0, accuracy=0):
    """
    Compute the geocentric position of a solar system body, as antedated
    for light-time.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date of observation.
    ss_object : Object
        Instance of Object type object containing the designation
        for the solar system body.
    pos_obs : tuple of floats, of length 3
        Position vector of observer (or the geocenter), with respect
        to origin at solar system barycenter, referred to ICRS axes,
        components in AU.
    estimate : float (optional)
        First approximation to light-time, in days (can be set to
        0.0 if unknown).
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position
        and time.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    tuple of tuple of floats, of length 3, and float
        Position vector of body, with respect to origin at observer
        (or the geocenter), referred to ICRS axes, components in AU
        and final light-time, in days.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))
    if estimate < 0.0:
        raise ValueError(_neg_err.format(name='estimate'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _light_time = novaslib.light_time
    _light_time.argtypes = (ctypes.c_double, ctypes.POINTER(Object),
                            ctypes.POINTER(ctypes.c_double*3),
                            ctypes.c_double, ctypes.c_short,
                            ctypes.POINTER(ctypes.c_double*3),
                            ctypes.POINTER(ctypes.c_double))
    _light_time.restype = ctypes.c_short
    _light_time.errcheck = _check_c_errors
    _light_time.c_errors = {
        1: (NonConvergentError, "from C function 'light_time': algorithm failed to converge after 10 iterations."),
        11: (ValueError, "from C function 'solarsystem': Invalid value of body or origin."),
        12: (RuntimeError, "from C function 'solarsystem': Error detected by JPL software.")
        }

    pos = (ctypes.c_double*3)()
    tlight = ctypes.c_double

    _light_time(jd_tdb, ctypes.byref(ss_object),
                ctypes.byref((ctypes.c_double*3)(*pos_obs)), estimate,
                accuracy, ctypes.byref(pos), ctypes.byref(tlight))

    return tuple([i for i in pos]), tlight.value


def d_light(pos_star, pos_obs):
    """
    Computes the difference in light-time, for a star, between the
    barycenter of the solar system and the observer (or the geocenter).

    Parameters
    ----------
    pos_star : tuple of floats, of length 3
        Position vector of star, with respect to origin at solar
        system barycenter.
    pos_obs : tuple of floats, of length 3
        Position vector of observer (or the geocenter), with respect
        to origin at solar system barycenter, components in AU.

    Returns
    -------
    diflt : float
        Difference in light time, in the sense star to barycenter
        minus star to earth, in days.

    Notes
    -----
    .. [N1] Alternatively, this function returns the light-time from
        the observer (or the geocenter) to a point on a light ray
        that is closest to a specific solar system body. For this
        purpose, 'pos1' is the position vector toward observed
        object, with respect to origin at observer (or the
        geocenter); 'pos_obs' is the position vector of solar system
        body, with respect to origin at observer (or the geocenter),
        components in AU; and the returned value is the light time
        to point on line defined by 'pos1' that is closest to solar
        system body (positive if light passes body before hitting
        observer, i.e., if 'pos1' is within 90 degrees of
        'pos_obs').

    """

    if len(pos_star) != 3:
        raise ValueError(_vector_len_err.format(name='pos_star'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))

    _d_light = novaslib.d_light
    _d_light.argtypes = (ctypes.POINTER(ctypes.c_double*3), ctypes.POINTER(ctypes.c_double*3))
    _d_light.restype = ctypes.c_double

    diflt = _d_light(ctypes.byref((ctypes.c_double*3)(*pos_star)),
                    ctypes.byref((ctypes.c_double*3)(*pos_obs)))

    return diflt


def grav_def(jd_tdb, pos_obj, pos_obs, location, accuracy=0):
    """
    Compute the total gravitational deflection of light for the observed
    object due to the major gravitating bodies in the solar system. This
    function valid for an observed body within the solar system as well
    as for a star.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date of observation.
    pos_obj : tuple of floats, of length 3
        Position vector of observed object, with respect to origin
        at observer (or the geocenter), referred to ICRS axes,
        components in AU.
    pos_obs : tuple of floats, of length 3
        Position vector of observer (or the geocenter), with respect
        to origin at solar system barycenter, referred to ICRS axes,
        components in AU.
    location : {0, 1}, optional
        Code for location of observer, determining whether the
        gravitational deflection due to the earth itself is applied.
            = 0 ... No earth deflection (normally means observer is
                at geocenter)
            = 1 ... Add in earth deflection (normally means observer
                is on or above surface of earth, including earth
                orbit)
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    positon : tuple of floats, of length 3
        Position vector of observed object, with respect to origin
        at observer (or the geocenter), referred to ICRS axes,
        corrected for gravitational deflection, components in AU.

    Notes
    -----
    .. [N1] If 'accuracy' is set to zero (full accuracy), three
        bodies (Sun, Jupiter, and Saturn) are used in the
        calculation. If the reduced-accuracy option is set, only the
        Sun is used in the calculation. In both cases, if the
        observer is not at the geocenter, the deflection due to the
        Earth is included.
    .. [N2] The number of bodies used at full and reduced accuracy
        can be set by making a change to the code in this function
        as indicated in the comments.

    References
    ----------
    .. [R1] Klioner, S. (2003), Astronomical Journal 125, 1580-1597,
        Section 6.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obj'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))
    if location not in [0, 1]:
        raise ValueError(_option_err.format(name='location', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _grav_def = novaslib.grav_def
    _grav_def.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short,
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3))
    _grav_def.restype = ctypes.c_short
    _grav_def.errcheck = _check_c_errors
    _grav_def.c_errors = {
        1: (ValueError, "from C function 'ephemeris': Invalid value of 'origin'."),
        2: (ValueError, "from C function 'ephemeris': Invalid value of 'type' in 'cel_obj'."),
        3: (MemoryError, "from C function 'ephemeris': Unable to allocate memory."),
        31: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        32: (ValueError, "from C function 'make_object': 'number' out of range"),
        33: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        34: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        35: (ValueError, "from C function 'make_object': 'name' is out of string bounds.")
        }

    pos2 = (ctypes.c_double*3)()

    _grav_def(jd_tdb, location, accuracy,
              ctypes.byref((ctypes.c_double*3)(*pos_obj)),
              ctypes.byref((ctypes.c_double*3)(*pos_obs)), ctypes.byref(pos2))

    return tuple([i for i in pos2])


def grav_vec(pos_obj, pos_obs, pos_body, rmass):
    """
    Correct the position vector for the deflection of light in the
    gravitational field of an arbitrary body. This function valid for an
    observed body within the solar system as well as for a star.

    Parameters
    ----------
    pos_obj : tuple of floats, of length 3
        Position vector of observed object, with respect to origin
        at observer (or the geocenter), components in AU.
    pos_obs : tuple of floats, of length 3
        Position vector of observer (or the geocenter), with respect
        to origin at solar system barycenter, components in AU.
    pos_body : tuple of floats, of length 3
        Position vector of gravitating body, with respect to origin
        at solar system barycenter, components in AU.
    rmass : float
        Reciprocal mass of gravitating body in solar mass units,
        that is, Sun mass / body mass.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector of observed object, with respect to origin
        at observer (or the geocenter), corrected for gravitational
        deflection, components in AU.

    References
    ----------
    .. [R1] Murray, C.A. (1981) Mon. Notices Royal Ast. Society 195,
        639-648.
    .. [R2] See also formulae in Section B of the Astronomical
        Almanac, or Kaplan, G. et al. (1989) Astronomical Journal
        97, 1197-1210, section iii f.

    """

    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obj'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))
    if len(pos_body) != 3:
        raise ValueError(_vector_len_err.format(name='pos_body'))
    if rmass < 0.0:
        raise ValueError(_neg_err.format(name='rmass'))

    _grav_vec = novaslib.grav_vec
    _grav_vec.argtypes = (ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3), ctypes.c_double,
                          ctypes.POINTER(ctypes.c_double*3))
    _grav_vec.restype = None

    pos2 = (ctypes.c_double*3)()

    _grav_vec(ctypes.byref((ctypes.c_double*3)(*pos_obj)),
                           ctypes.byref((ctypes.c_double*3)(*pos_obs)),
                           ctypes.byref((ctypes.c_double*3)(*pos_body)), rmass,
                           ctypes.byref(pos2))

    return tuple([i for i in pos2])


def aberration(position, vel_earth, lighttime):
    """
    Correct the position vector for aberration of light, including
    relativistic terms.

    Parameters
    ----------
    position : tuple of floats, of length 3
        Position vector, referred to origin at center of mass of the
        Earth, components in AU.
    vel_earth : tuple of floats, of length 3
        Velocity vector of center of mass of the Earth, referred to
        origin at solar system barycenter, components in AU/day.
    lighttime : float
        Light time from object to Earth in days.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, referred to origin at center of mass of the
        Earth, corrected for aberration, components in AU.

    Notes
    -----
    .. [N1] If 'lighttime' = 0 on input, this function will compute
        it.

    References
    ----------
    .. [R1] Murray, C. A. (1981) Mon. Notices Royal Ast. Society
        195, 639-648.

    """

    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if len(vel_earth) != 3:
        raise ValueError(_vector_len_err.format(name='vel_earth'))
    if lighttime < 0.0:
        raise ValueError(_neg_err.format(name='lighttime'))

    _aberration = novaslib.aberration
    _aberration.argtypes = (ctypes.POINTER(ctypes.c_double*3),
                            ctypes.POINTER(ctypes.c_double*3), ctypes.c_double,
                            ctypes.POINTER(ctypes.c_double*3))
    _aberration.restype = None

    pos2 = (ctypes.c_double*3)()

    _aberration(ctypes.byref((ctypes.c_double*3)(*position)),
                ctypes.byref((ctypes.c_double*3)(*vel_earth)), lighttime,
                ctypes.byref(pos2))

    return tuple([i for i in pos2])


def rad_vel(cel_object, pos, vel, vel_obs, distance_geo, distance_sun,
            distance_obj_sun):
    """
    Predicts the radial velocity of the observed object as it would be
    measured by spectroscopic means. Radial velocity is here defined as
    the radial velocity measure (z) times the speed of light. For a
    solar system body, it applies to a fictitious emitter at the center
    of the observed object, assumed massless (no gravitational red
    shift), and does not in general apply to reflected light. For stars,
    it includes all effects, such as gravitational red shift, contained
    in the catalog barycentric radial velocity measure, a scalar derived
    from spectroscopy. Nearby stars with a known kinematic velocity
    vector (obtained independently of spectroscopy) can be treated like
    solar system objects.

    Parameters
    ----------
    cel_object : Object
        Specifies the celestial object of interest.
    pos : tuple of floats, of length 3
        Geometric position vector of object with respect to
        observer, corrected for light-time, in AU.
    vel : tuple of floats, of length 3
        Velocity vector of object with respect to solar system
        barycenter, in AU/day.
    vel_obs : tuple of floats, of length 3
        Velocity vector of observer with respect to solar system
        barycenter, in AU/day.
    distance_geo : float
        Distance from observer to geocenter, in AU.
    distance_sun : float
        Distance from observer to Sun, in AU.
    distance_obj_sun : float
        Distance from object to Sun, in AU.

    Returns
    -------
    rv : float
        The observed radial velocity measure times the speed of
        light, in kilometers/second.

    Notes
    -----
    .. [N1] All the input arguments are BCRS quantities, expressed
        with respect to the ICRS axes. 'vel' and 'vel_obs' are
        kinematic velocities - derived from geometry or dynamics,
        not spectroscopy.
    .. [N2] If the object is outside the solar system, the algorithm
        used will be consistent with the IAU definition of stellar
        radial velocity, specifically, the barycentric radial
        velocity measure, which is derived from spectroscopy. In
        that case, the vector 'vel' can be very approximate -- or,
        for distant stars or galaxies, zero -- as it will be used
        only for a small geometric correction that is proportional
        to proper motion.
    .. [N3] Any of the distances (last three input arguments) can be
        set to zero (0.0) if the corresponding general relativistic
        gravitational potential term is not to be evaluated. These
        terms generally are important only at the meter/second
        level. If 'distance_geo' and 'distance_sun' are both zero,
        an average value will be used for the relativistic term for
        the observer, appropriate for an observer on the surface of
        the Earth. 'distance_sun', if given, is used only for solar
        system objects.

    References
    ----------
    .. [R1] Lindegren & Dravins (2003), Astronomy & Astrophysics
        401, 1185-1201.

    """

    if len(pos) != 3:
        raise ValueError(_vector_len_err.format(name='pos'))
    if len(vel) != 3:
        raise ValueError(_vector_len_err.format(name='vel'))
    if len(vel_obs) != 3:
        raise ValueError(_vector_len_err.format(name='vel_obs'))
    if distance_geo < 0.0:
        raise ValueError(_neg_err.format(name='distance_geo'))
    if distance_sun < 0.0:
        raise ValueError(_neg_err.format(name='distance_sun'))
    if distance_obj_sun < 0.0:
        raise ValueError(_neg_err.format(name='distance_obj_sun'))

    _rad_vel = novaslib.rad_vel
    _rad_vel.argtypes = (ctypes.POINTER(Object),
                         ctypes.POINTER(ctypes.c_double*3),
                         ctypes.POINTER(ctypes.c_double*3),
                         ctypes.POINTER(ctypes.c_double*3), ctypes.c_double,
                         ctypes.c_double, ctypes.c_double,
                         ctypes.POINTER(ctypes.c_double))
    _rad_vel.restype = None

    rv = ctypes.c_double()

    _rad_vel(ctypes.byref(cel_object), ctypes.byref((ctypes.c_double*3)(*pos)),
             ctypes.byref((ctypes.c_double*3)(*vel)),
             ctypes.byref((ctypes.c_double*3)(*vel_obs)),
             distance_geo, distance_sun, distance_obj_sun, ctypes.byref(rv))

    return rv.value


def precession(jd_tdb1, position, jd_tdb2):
    """
    Precesses equatorial rectangular coordinates from one epoch to
    another. One of the two epochs must be J2000.0. The coordinates
    are referred to the mean dynamical equator and equinox of the two
    respective epochs.

    Parameters
    ----------
    jd_tdb1 : float
        TDB Julian date of first epoch. See Note [N1]_ below.
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to mean dynamical equator and equinox
        of first epoch.
    jd_tdb2 : float
        TDB Julian date of second epoch. See Note [N1]_ below.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to mean dynamical equator and equinox
        of second epoch.

    Notes
    -----
    .. [N1] Either 'date' or 'newdate' must be 2451545.0 (J2000.0)
        TDB.

    References
    ----------
    .. [R1] Explanatory Supplement To The Astronomical Almanac, pp.
        103-104.
    .. [R2] Capitaine, N. et al. (2003), Astronomy And Astrophysics
        412, pp. 567-586.
    .. [R3] Hilton, J. L. et al. (2006), IAU WG report, Celest.
        Mech., 94, pp. 351-367.

    """

    if jd_tdb1 < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb1'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if jd_tdb2 < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb2'))

    _precession = novaslib.precession
    _precession.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double*3), ctypes.c_double,
                            ctypes.POINTER(ctypes.c_double*3))
    _precession.restype = ctypes.c_short
    _precession.errcheck = _check_c_errors
    _precession.c_errors = {
        1: (ValueError, "from C function 'precession': Precession not to or from J2000.0; 'jd_tdb1' or 'jd_tdb2' not 2451545.0.")
        }

    pos2 = (ctypes.c_double*3)()

    _precession(jd_tdb1, ctypes.byref((ctypes.c_double*3)(*position)),
                jd_tdb2, ctypes.byref(pos2))

    return tuple([i for i in pos2])


def nutation(jd_tdb, position, direction=0, accuracy=0):
    """
    Nutates equatorial rectangular coordinates from mean equator and
    equinox of epoch to true equator and equinox of epoch. Inverse
    transformation may be applied by setting flag 'direction'.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date of epoch.
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to mean equator and equinox of epoch.
    direction : {0, 1}, optional
            = 0 ... transform from mean to true (default)
            = 1 ... transform from true to mean
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular
        coordinates, referred to true equator and equinox of epoch.

    References
    ----------
    .. [R1] Explanatory Supplement To The Astronomical Almanac, pp.
        114-115.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if direction not in [0, 1]:
        raise ValueError(_option_err.format(name='direction', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _nutation = novaslib.nutation
    _nutation.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short,
                          ctypes.POINTER(ctypes.c_double*3),
                          ctypes.POINTER(ctypes.c_double*3))
    _nutation.restype = None

    pos2 = (ctypes.c_double*3)()

    _nutation(jd_tdb, direction, accuracy,
              ctypes.byref((ctypes.c_double*3)(*position)), ctypes.byref(pos2))

    return tuple([i for i in pos2])


def nutation_angles(time, accuracy=0):
    """
    This function returns the values for nutation in longitude and
    nutation in obliquity for a given TDB Julian date. The nutation
    model selected depends upon the input value of 'accuracy'. See
    notes below for important details.

    Parameters
    ----------
    time : float
        TDB time in Julian centuries since J2000.0
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        nutation.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (dpsi, deps) : tuple of floats
        Nutation in (longitude, obliquity) in arcseconds.

    Notes
    -----
    .. [N1] This function selects the nutation model depending first
        upon the input value of 'accuracy'. If 'accuracy' = 0 (full
        accuracy), the IAU 2000A nutation model is used. If
        'accuracy' = 1 (reduced accuracy, the model used depends
        upon the value of local variable 'low_acc_choice', which is
        set in 'novas.c'.
    .. [N2] If local variable 'low_acc_choice' = 1 (the default), a
        specially truncated version of IAU 2000A, called 'NU2000K'
        is used. If 'low_acc_choice' = 2, the IAU 2000B nutation
        model is used.
    .. [N3] See the docstrings of the nutation functions in
        novas.nutation for details concerning the models.

    References
    ----------
    .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _nutation_angles = novaslib.nutation_angles
    _nutation_angles.argtypes = (ctypes.c_double, ctypes.c_short,
                                 ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double))
    _nutation_angles.restype = None

    dpsi = ctypes.c_double()
    deps = ctypes.c_double()

    _nutation_angles(time, accuracy, ctypes.byref(dpsi), ctypes.byref(deps))

    return (dpsi.value, deps.value)


def fund_args(time):
    """
    To compute the fundamental arguments (mean elements) of the Sun and
    Moon.

    Parameters
    ----------
    time : float
        TDB time in Julian centuries since J2000.0

    Returns
    -------
    (l, l', F, D, omega) : tuple
        l = mean anomaly of the Moon
        l' = mean anomaly of the Sun
        F = mean argument of the latitude of the Moon
        D = mean elongation of the Moon from the Sun
        omega = mean longitude of the Moon's ascending node;
            from Simon section 3.4(b.3),
            precession = 5028.8200 arcsec/cy

    References
    ----------
    .. [R1] Simon et al. (1994) Astronomy and Astrophysics 282, 663-683,
        esp. Sections 3.4-3.5.

    """

    _fund_args = novaslib.fund_args
    _fund_args.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double*5))
    _fund_args.restype = None

    a = (ctypes.c_double*5)()

    _fund_args(time, ctypes.byref(a))

    return tuple([i for i in a])


def mean_obliq(jd_tdb):
    """
    Return the mean obliquity of the ecliptic.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date.

    Returns
    -------
    epsilon : float
        Mean obliquity of the ecliptic in arcseconds.

    References
    ----------
    .. [R1] Capitaine et al. (2003), Astronomy and Astrophysics 412,
        567-586.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))

    _mean_obliq = novaslib.mean_obliq
    _mean_obliq.argtypes = (ctypes.c_double,)
    _mean_obliq.restype = ctypes.c_double

    obliq = _mean_obliq(jd_tdb)

    return obliq


def vector2radec(position):
    """
    Converts an vector in equatorial rectangular coordinates to
    equatorial spherical coordinates.

    Parameters
    ----------
    position : tuple of floats, of length 3
        Position vector, equatorial rectangular coordinates.

    Returns
    -------
    (rightascension, declination) : tuple of floats
        (Right ascension in hours, declination in degrees)

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
        C3.1', C38-C39.

    """

    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))

    _vector2radec = novaslib.vector2radec
    _vector2radec.argtypes = (ctypes.POINTER(ctypes.c_double*3),
                              ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_double))
    _vector2radec.restype = ctypes.c_short
    _vector2radec.errcheck = _check_c_errors
    _vector2radec.c_errors = {
        1: (IndeterminateError, "from C function 'vector2radec': All vector components are zero; 'ra' and 'dec' are indeterminate."),
        2: (IndeterminateError, "from C function 'vector2radec': Both pos[0] and pos[1] are zero, but pos[2] is nonzero; 'ra' is indeterminate.")
        }

    ra = ctypes.c_double()
    dec = ctypes.c_double()

    _vector2radec(ctypes.byref((ctypes.c_double*3)(*position)),
                  ctypes.byref(ra), ctypes.byref(dec))

    return (ra.value, dec.value)


def radec2vector(ra, dec, distance):
    """
    Converts equatorial spherical coordinates to a vector (equatorial
    rectangular coordinates).

    Parameters
    ----------
    ra : float
        Right ascension in hours.
    dec : float
        Declination in degrees.
    distance : float
        Distance in AU.

    Returns
    -------
    position : tuple of floats, of length 3
        Position vector, equatorial rectangular coordinates in AU.

    """

    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err.format(name='rai'))
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err.format(name='deci'))
    if distance < 0.0:
        raise ValueError(_neg_err.format(name='distance'))

    _radec2vector = novaslib.radec2vector
    _radec2vector.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double,
                              ctypes.POINTER(ctypes.c_double*3))
    _radec2vector.restype = None

    vector = (ctypes.c_double*3)()

    _radec2vector(ra, dec, distance, ctypes.byref(vector))

    return tuple([i for i in vector])


def starvectors(star):
    """
    Converts angular quantities for stars to vectors.

    Parameters
    ----------
    star : CatEntry
        Instance of CatEntry type object containing ICRS catalog
        data.

    Returns
    -------
    (position, velocity) : tuple of tuple of floats, of length 3
        Position and velocity vectors in equatorial rectangular
        coordinates. Components are in AU and AU/day.

    References
    ----------
    .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
        Version C3.1', C16, C113.

    """

    _starvectors = novaslib.starvectors
    _starvectors.argtypes = (ctypes.POINTER(CatEntry),
                             ctypes.POINTER(ctypes.c_double*3),
                             ctypes.POINTER(ctypes.c_double*3))
    _starvectors.restype = None

    pos = (ctypes.c_double*3)()
    vel = (ctypes.c_double*3)()

    _starvectors(ctypes.byref(star), ctypes.byref(pos), ctypes.byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])


def tdb2tt(tdb_jd):
    """
    Computes the Terrestrial Time (TT) or Terrestrial Dynamical Time
    (TDT) Julian date corresponding to a Barycentric Dynamical Time
    (TDB) Julian date.

    Parameters
    ----------
    tdb_jd : float
        TDB Julian date.

    Returns
    -------
    (tt_jd, secdiff) : tuple of float
        tt_jd = TT Julian date, secdiff = difference 'tdb_jd' -
        'tt_jd', in seconds.

    Notes
    -----
    .. [N1] Expression used in this function is a truncated form of
        a longer and more precise series given in the first
        reference [R1]_. The result is good to about 10
        microseconds.

    References
    ----------
    .. [R1] Fairhead, L. & Bretagnon, P. (1990) Astron. & Astrophys.
        229, 240.
    .. [R2] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if tdb_jd < 0.0:
        raise ValueError(_neg_err.format(name='tdb_jd'))

    _tdb2tt = novaslib.tdb2tt
    _tdb2tt.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double))
    _tdb2tt.restype = None

    tt_jd = ctypes.c_double()
    secdiff = ctypes.c_double()

    _tdb2tt(tdb_jd, ctypes.byref(tt_jd), ctypes.byref(secdiff))

    return tt_jd.value, secdiff.value


def cio_ra(jd_tt, accuracy=0):
    """
    Returns the true right ascension of the celestial intermediate
    origin (CIO) at a given TT Julian date. This is -(equation of the
    origins).

    Parameters
    ----------
    jd_tt : float
        TT Julian day.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    ra_cio : float
        Right ascension of the CIO, with respect to the true equinox
        of date, in hours (+ or -).

    References
    ----------
    .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if jd_tt < 0.0:
        raise ValueError(_neg_err.format(name='jd_tt'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _cio_ra = novaslib.cio_ra
    _cio_ra.argtypes = (ctypes.c_double, ctypes.c_short,
                        ctypes.POINTER(ctypes.c_double))
    _cio_ra.restype = ctypes.c_short
    _cio_ra.errcheck = _check_c_errors
    _cio_ra.c_errors = {
        1: (ValueError, "from C function 'cio_ra': invalid value of 'accuracy'."),
        11: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        21: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        22: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        23: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        24: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        25: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        26: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure."),
#        21: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    ra_cio = ctypes.c_double()

    _cio_ra(jd_tt, accuracy, ctypes.byref(ra_cio))

    return ra_cio.value


def cio_location(jd_tdb, accuracy=0):
    """
    Returns the location of the celestial intermediate origin (CIO) for
    a given Julian date, as a right ascension with respect to either the
    GCRS (geocentric ICRS) origin or the true equinox of date.  The CIO
    is always located on the true equator (= intermediate equator) of
    date.

    Parameters
    ----------
    jd_tdb : float
         TDB Julian date.
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    ra_cio : float
        Right ascension of the CIO, in hours.
    ref_sys : {1, 2}
        Reference system in which right ascension is given
            = 1 ... GCRS
            = 2 ... True equator and equinox of date.

    Notes
    -----
    .. [N1] If an external file of CIO right ascensions is available, it
        will be used and 'ref_sys' will be set to 1. Otherwise an
        internal computation will be used and 'ref_sys' will be set to 2.
    .. [N2] The external binary file of CIO right ascensions is assumed
        to be named 'cio_ra.bin'.  Utility program 'cio_file.c',
        provided with the NOVAS-C package, creates this file from a text
        file also provided with NOVAS-C.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _cio_location = novaslib.cio_location
    _cio_location.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                              ctypes.POINTER(ctypes.c_short))
    _cio_location.restype = ctypes.c_short
    _cio_location.errcheck = _check_c_errors
    _cio_location.c_errors = {
        1: (MemoryError, "from C function 'cio_location': unable to allocate memory for the 'cio' array."),
        11: (IOError, "from C function 'cio_array': error opening the 'cio_ra.bin' file."),
        12: (ValueError, "from C function 'cio_array': 'jd_tdb' not in the range of the CIO file."),
        13: (ValueError, "from C function 'cio_array': 'n_pts' out of range."),
        14: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 't' array."),
        15: (MemoryError, "from C function 'cio_array': unable to allocate memory for the internal 'ra' array."),
        16: (ValueError, "from C function 'cio_array': 'jd_tdb' is too close to either end of the CIO file; unable to put 'n_pts' data points into the output ctypes.Structure.")
        }

    ra_cio = ctypes.c_double()
    ref_sys = ctypes.c_short()

    _cio_location(jd_tdb, accuracy, ctypes.byref(ra_cio),
                  ctypes.byref(ref_sys))

    return ra_cio.value, ref_sys.value


def cio_basis(jd_tdb, ra_cio, system, accuracy=0):
    """
    Compute the orthonormal basis vectors, with respect to the GCRS
    (geocentric ICRS), of the celestial intermediate system defined by
    the celestial intermediate pole (CIP) (in the z direction) and the
    celestial intermediate origin (CIO) (in the x direction). A TDB
    Julian date and the right ascension of the CIO at that date is
    required as input. The right ascension of the CIO can be with
    respect to either the GCRS origin or the true equinox of date --
    different algorithms are used in the two cases.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian date of epoch.
    ra_cio : float
        Right ascension of the CIO at epoch in hours.
    system : {1, 2}, optional
        Reference system in which right ascension is given (output
        from function 'cio_location')
            = 1 ... GCRS
            = 2 ... True equator and equinox of date
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (x, y, z) : tuple of floats
        Unit vector toward (the CIO, the y-direction, north
        celestial pole (CIP)), equatorial rectangular coordinates,
        referred to the GCRS

    Notes
    -----
    .. [N1] This function effectively constructs the matrix C in eq.
        (3) of the reference.

    References
    ----------
    .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if not 0.0 <= ra_cio < 24.0:
        raise ValueError(_hour_range_err.format(name='ra_cio'))
    if system not in [1, 2]:
        raise ValueError(_option_err.format(name='system', allowed=[1, 2]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _cio_basis = novaslib.cio_basis
    _cio_basis.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_short,
                           ctypes.c_short, ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double))
    _cio_basis.restype = ctypes.c_short
    _cio_basis.errcheck = _check_c_errors
    _cio_basis.c_errors = {
        1: (ValueError, "from C function 'cio_basis': invalid value of input variable 'ref_sys'.")
        }

    x = ctypes.c_double()
    y = ctypes.c_double()
    z = ctypes.c_double()

    _cio_basis(jd_tdb, ra_cio, system, accuracy, ctypes.byref(x),
               ctypes.byref(y), ctypes.byref(z))

    return x.value, y.value, z.value


def ira_equinox(jd_tdb, equinox, accuracy=0):
    """
    To compute the intermediate right ascension of the equinox at
    the input Julian date, using an analytical expression for the
    accumulated precession in right ascension. For the true equinox,
    the result is the equation of the origins.

    Parameters
    ----------
    jd_tdb : float
        TDB Julian day.
    equinox : {0, 1}, optional
        = 0 ... mean equinox
        = 1 ... true equinox
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output
        position.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    ira_eq : float
        Intermediate right ascension of the equinox, in hours (+ or
        -). If 'equinox' = 1 (i.e. true equinox), then the returned
        value is the equation of the origins.

    References
    ----------
    .. [R1] Capitaine, N. et al. (2003), Astronomy and Astrophysics
        412, 567-586, eq. (42).

    """

    if jd_tdb < 0.0:
        raise ValueError(_neg_err.format(name='jd_tdb'))
    if equinox not in [0, 1]:
        raise ValueError(_option_err.format(name='equinox', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _ira_equinox = novaslib.ira_equinox
    _ira_equinox.argtypes = (ctypes.c_double, ctypes.c_short, ctypes.c_short)
    _ira_equinox.restype = ctypes.c_double

    ira_eq = _ira_equinox(jd_tdb, equinox, accuracy)

    return ira_eq


def ephemeris(jd, ss_body, origin=1, accuracy=0):
    """

    Retrieves the position and velocity of a solar system body from
    a fundamental ephemeris.

    Parameters
    ----------
    jd : tuple of floats, of length 2
        TDB Julian date split into two parts, where the sum jd[0] +
        jd[1] is the TDB Julian date.
    ss_body : Object
        Instance of Object type object containing the designation of
        the body of interest.
    origin : {0, 1}, optional
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun (default)
    accuracy : {0, 1}, optional
        Code specifying the relative accuracy of the output position
        and velocity.
            = 0 ... full accuracy (default)
            = 1 ... reduced accuracy

    Returns
    -------
    (pos, vel) : tuple of tuples of three floats
        pos is position vector of the body at 'jd_tdb'; equatorial
        rectangular coordinates in AU referred to the ICRS.
        vel is the velocity vector of the body at 'jd_tdb';
        equatorial rectangular system referred the ICRS, in AU/Day.

    """

    if jd[0] + jd[1] < 0.0:
        raise ValueError(_neg_err.format(name='jd[0]+jd[1]'))
    if origin not in [0, 1]:
        raise ValueError(_option_err.format(name='origin', allowed=[0, 1]))
    if accuracy not in [0, 1]:
        raise ValueError(_option_err.format(name='accuracy', allowed=[0, 1]))

    _ephemeris = novaslib.ephemeris
    _ephemeris.argtypes = (ctypes.c_double*2, ctypes.POINTER(Object),
                           ctypes.c_short, ctypes.c_short,
                           ctypes.POINTER(ctypes.c_double*3),
                           ctypes.POINTER(ctypes.c_double*3))
    _ephemeris.restype = ctypes.c_short
    _ephemeris.errcheck = _check_c_errors
    _ephemeris.c_errors = {
        1: (ValueError, "from C function 'ephemeris': Invalid value of 'origin'."),
        2: (ValueError, "from C function 'ephemeris': Invalid value of 'type' in 'cel_obj'."),
        3: (MemoryError, "from C function 'ephemeris': Unable to allocate memory."),
        11: (ValueError, "from C function 'solarsystem': Invalid value of body or origin."),
        12: (RuntimeError, "from C function 'solarsystem': Error detected by JPL software."),
        21: (MemoryError, "from C function 'readeph': memory allocation error"),
        22: (ValueError, "from C function 'readeph': mismatch between asteroid name and number"),
        23: (ValueError, "from C function 'readeph': Julian date out of bounds"),
        24: (IOError, "from C function 'readeph': can not find Chebyshev polynomial file"),
        29: (RuntimeError, "from C function 'readeph': dummy 'readeph' function called (see note 1 in readeph0.c)")
        }

    pos = (ctypes.c_double*3)()
    vel = (ctypes.c_double*3)()

    _ephemeris((ctypes.c_double*2)(*jd), ctypes.byref(ss_body),
               origin, accuracy, ctypes.byref(pos), ctypes.byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])


def transform_hip(hipparcos):
    """
    To convert Hipparcos catalog data at epoch J1991.25 to epoch
    J2000.0, for use within NOVAS. To be used only for Hipparcos or
    Tycho stars with linear space motion. Both input and output data
    is in the ICRS.

    Parameters
    ----------
    hipparcos : CatEntry
        Instance of CatEntry type object containing an entry from
        the Hipparcos catalog, at epoch J1991.25, with all members
        having Hipparcos catalog units. See Note [N1]_ below.

    Returns
    -------
    hip_2000 : CatEntry
        Instance of CatEntry type object containing the transformed
        input entry, at epoch J2000.0. See Note [N2]_ below.

    Notes
    -----
    .. [N1] Input (Hipparcos catalog) epoch and units:
        Epoch: J1991.25
        Right ascension (RA): degrees
        Declination (Dec): degrees
        Proper motion in RA: milliarcseconds per year
        Proper motion in Dec: milliarcseconds per year
        Parallax: milliarcseconds
        Radial velocity: kilometers per second (not in catalog)
    .. [N2] Output (modified Hipparcos) epoch and units:
        Epoch: J2000.0
        Right ascension: hours
        Declination: degrees
        Proper motion in RA: milliarcseconds per year
        Proper motion in Dec: milliarcseconds per year
        Parallax: milliarcseconds
        Radial velocity: kilometers per second

    """

    _transform_hip = novaslib.transform_hip
    _transform_hip.argtypes = (ctypes.POINTER(CatEntry),
                               ctypes.POINTER(CatEntry))
    _transform_hip.restype = None

    hip_2000 = CatEntry()

    _transform_hip(ctypes.byref(hipparcos), ctypes.byref(hip_2000))

    return hip_2000


def transform_cat(option, date_incat, incat, date_newcat, newcat_id):
    """
    To transform a star's catalog quantities for a change of epoch
    and/or equator and equinox. Also used to rotate catalog
    quantities on the dynamical equator and equinox of J2000.0 to the
    ICRS or vice versa.

    Parameters
    ----------
    option : {1, 2, 3, 4, 5}
        Transformation option
            = 1 ... change epoch; same equator and equinox
            = 2 ... change equator and equinox; same epoch
            = 3 ... change equator and equinox and epoch
            = 4 ... change equator and equinox J2000.0 to ICRS
            = 5 ... change ICRS to equator and equinox of J2000.0
    date_incat : float
        TT Julian date, or year, of input catalog data.
    incat : CatEntry
        Instance of CatEntry type object containing an entry from
        the input catalog, with units as given in the type
        definition.
    date_newcat : float
        TT Julian date, or year, of transformed catalog data.
    newcat_id : string
        Abbreviated name of the transformed catalog. Max length
        determined by SIZE_OF_CAT_NAME defined in ``novas.h``;
        defaults to 3.

    Returns
    -------
    newcat : CatEntry
        Instance of CatEntry type object containing the transformed
        catalog entry, with units as given in the type definition.

    Notes
    -----
    .. [N1] 'date_incat' and 'date_newcat' may be specified either
        as a Julian date (e.g., 2433282.5) or a Julian year and
        fraction (e.g., 1950.0). Values less than 10000 are assumed
        to be years. For 'option' = 2 or 'option' = 3, either
        'date_incat' or 'date_newcat' must be 2451545.0 or 2000.0
        (J2000.0). For 'option' = 4 and 'option' = 5, 'date_incat'
        and 'date_newcat' are ignored.
    .. [N2] 'option' = 1 updates the star's data to account for the
        star's space motion between the first and second dates,
        within a fixed reference frame.
            'option' = 2 applies a rotation of the reference frame
        corresponding to precession between the first and second
        dates, but leaves the star fixed in space.
            'option' = 3 provides both transformations.
            'option' = 4 and 'option' = 5 provide a a fixed rotation
        about very small angles (<0.1 arcsecond) to take data from
        the dynamical system of J2000.0 to the ICRS ('option' = 4)
        or vice versa ('option' = 5).
    .. [N3] For 'option' = 1, input data can be in any fixed
        reference system. for 'option' = 2 or 'option' = 3, this
        function assumes the input data is in the dynamical system
        and produces output in the dynamical system. for 'option' =
        4, the input data must be on the dynamical equator and
        equinox of J2000.0. for 'option' = 5, the input data must be
        in the ICRS.
    .. [N4] This function cannot be properly used to bring data from
        old star catalogs into the modern system, because old
        catalogs were compiled using a set of constants that are
        incompatible with modern values. In particular, it should
        not be used for catalogs whose positions and proper motions
        were derived by assuming a precession constant significantly
        different from the value implicit in function 'precession'.

    """

    if option not in [1, 2, 3, 4, 5]:
        raise ValueError(_option_err.format(name='option',
                                            allowed=[1, 2, 3, 4, 5]))
    if date_incat < 0.0:
        raise ValueError(_neg_err.format(name='date_incat'))
    if date_newcat < 0.0:
        raise ValueError(_neg_err.format(name='date_newcat'))

    _transform_cat = novaslib.transform_cat
    _transform_cat.argtypes = (ctypes.c_short, ctypes.c_double,
                               ctypes.POINTER(CatEntry), ctypes.c_double,
                               ctypes.c_char*4, ctypes.POINTER(CatEntry))
    _transform_cat.restype = ctypes.c_short
    _transform_cat.errcheck = _check_c_errors
    _transform_cat.c_errors = {
        1: (ValueError, "from C function 'transform_cat': Invalid value of an input date for option 2 or 3 (see Note 1)."),
        2: (ValueError, "from C function 'transform_cat': length of 'newcat_id' out of bounds.")
        }

    newcat = CatEntry()

    _transform_cat(option, date_incat, ctypes.byref(incat),
                   date_newcat,
                   ctypes.create_string_buffer(newcat_id.encode('ascii'), 4),
                   ctypes.byref(newcat))

    return newcat


def limb_angle(pos_obj, pos_obs):
    """
    Compute the angle of an object above or below the Earth's limb
    (horizon). The geometric limb is computed, assuming the Earth to be
    an airless sphere (no refraction or oblateness is included). The
    observer can be on or above the Earth. For an observer on the
    surface of the Earth, this function returns the approximate
    unrefracted altitude.

    Parameters
    ----------
    pos_obj : tuple of floats, of length 3
        Position vector of observed object, with respect to origin
        at geocenter, components in AU.
    pos_obs : tuple of floats, of length 3
        Position vector of observer, with respect to origin at
        geocenter, components in AU.

    Returns
    -------
    (limb_angle, nadir_angle) : tuple of floats
        Angle of observed object above (+) or below (-) limb in
        degrees and nadir angle of observed object as a fraction of
        apparent radius of limb where: < 1.0, below the limb; = 1.0,
        on the limb; > 1.0, above the limb

    """

    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obj'))
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err.format(name='pos_obs'))

    _limb_angle = novaslib.limb_angle
    _limb_angle.argtypes = (ctypes.c_double*3, ctypes.c_double*3,
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double))
    _limb_angle.restype = None

    limb_ang = ctypes.c_double()
    nadir_ang = ctypes.c_double()

    _limb_angle((ctypes.c_double*3)(*pos_obj), (ctypes.c_double*3)(*pos_obs),
                ctypes.byref(limb_ang), ctypes.byref(nadir_ang))

    return limb_ang.value, nadir_ang.value


def refract(location, zd_obs, atmosphere=1):
    """

    This function computes atmospheric refraction in zenith
    distance. This version computes approximate refraction for
    optical wavelengths.

    Parameters
    ----------
    location : OnSurface
        Instance of OnSurface type object  containing observer's
        location. This ctypes.Structure also contains weather data
        (optional) for the observer's location.
    zenith : float
        Observed zenith distance, in degrees.
    atmosphere : {1, 2}, optional
            = 1 ... use 'standard' atmospheric conditions (default)
            = 2 ... use atmospheric parameters input in the
                'location' object.

    Returns
    -------
    refraction : float
        Atmospheric refraction, in degrees.

    Notes
    -----
    .. [N1] This function can be used for planning observations or
        telescope pointing, but should not be used for the reduction
        of precise observations.

    References
    ----------
    .. [R1] Explanatory Supplement to the Astronomical Almanac,
        p. 144.
    .. [R2] Bennett, G. (1982), Journal of Navigation (Royal
        Institute) 35, pp. 255-259.

    """

    if atmosphere not in [1, 2]:
        raise ValueError(_option_err.format(name='atmosphere', allowed=[1, 2]))

    _refract = novaslib.refract
    _refract.argtypes = (ctypes.POINTER(OnSurface), ctypes.c_short,
                         ctypes.c_double)
    _refract.restype = ctypes.c_double

    refraction = _refract(ctypes.byref(location), atmosphere, zd_obs)

    return refraction


def julian_date(year, month, day, hour=0.0):
    """
    This function will compute the Julian date for a given calendar
    date (year, month, day, hour).

    Parameters
    ----------
    year : integer
        Year.
    month : integer
        Month number.
    day : integer
        Day-of-month.
    hour : float, optional
        Hour-of-day.

    Returns
    -------
    jd : float
        Julian day.

    Notes
    -----
    .. [N1] This function makes no checks for a valid input calendar
        date.
    .. [N2] Input calendar date must be Gregorian.
    .. [N3] Input time value can be based on any UT-like time scale
        (UTC, UT1, TT, etc.) -- output Julian date will have the
        same basis.

    References
    ----------
    .. [R1] Fliegel, H. & Van Flandern, T. Comm. of the ACM,
        Vol. 11, No. 10, October 1968, p. 657.

    """

    if year < -4712:
        raise ValueError(_jd_year_err.format(name='year'))
    if not 1 <= month <= 12:
        raise ValueError(_month_range_err.format(name='month'))
    if not 1 <= day <= _days_in_month[month]:
        raise ValueError(_day_range_err.format(name='day', daysinmonth=_days_in_month[month]))
    if not 0.0 <= hour < 24.0:
        raise ValueError(_hour_range_err.format(name='hour'))

    _julian_date = novaslib.julian_date
    _julian_date.argtypes = (ctypes.c_short, ctypes.c_short, ctypes.c_short,
                             ctypes.c_double)
    _julian_date.restype = ctypes.c_double

    jd = _julian_date(year, month, day, hour)

    return jd


def cal_date(jd):
    """
    Return the Gregorian date for a given Julian day.

    Parameters
    ----------
    jd : float
        Julian day.

    Returns
    -------
    date : tuple of three ints and a float
        The elements are year, month, day and hour.

    Notes
    -----
    .. [N1] This routine valid for any 'jd' greater than zero.
    .. [N2] Input Julian date can be based on any UT-like time scale
        (UTC, UT1, TT, etc.) -- output time value will have same
        basis.

    References
    ----------
    .. [R1] Fliegel, H. & Van Flandern, T. Comm. of the ACM,
        Vol. 11, No. 10, October 1968, p. 657.

    """

    if jd < 0.0:
        raise ValueError(_neg_err.format(name='jd'))

    _cal_date = novaslib.cal_date
    _cal_date.argtypes = (ctypes.c_double, ctypes.POINTER(ctypes.c_short),
                          ctypes.POINTER(ctypes.c_short),
                          ctypes.POINTER(ctypes.c_short),
                          ctypes.POINTER(ctypes.c_double))
    _cal_date.restype = None

    year = ctypes.c_short()
    month = ctypes.c_short()
    day = ctypes.c_short()
    hour = ctypes.c_double()

    _cal_date(jd, ctypes.byref(year), ctypes.byref(month), ctypes.byref(day),
              ctypes.byref(hour))

    return year.value, month.value, day.value, hour.value


def norm_ang(angle):
    """
    Normalize angle into the range 0 <= angle < (2 * pi).

    Parameters
    ----------
    angle : float
        Input angle in radians.

    Returns
    -------
    norm_angle : float
        The input angle, normalized as described above, in radians.

    """

    _norm_ang = novaslib.norm_ang
    _norm_ang.argtypes = (ctypes.c_double,)
    _norm_ang.restype = ctypes.c_double

    norm_angle = _norm_ang(angle)

    return norm_angle


def make_cat_entry(star_name, catalog, star_num, ra, dec, pm_ra, pm_dec,
                   parallax, rad_vel):
    """
    Create an instance of CatEntry containing catalog data for a star or
    "star-like" object.

    Parameters
    ----------
    star_name : string
        Object name. Max length determined by SIZE_OF_OBJ_NAME
        defined in ``novas.h``; defaults to 50.
    catalog : string
        Catalog identifier (e.g. HIP = Hipparcos, TY2 = Tycho-2).
        Max length determined by SIZE_OF_CAT_NAME defined in
        ``novas.h``; defaults to 3.
    star_num : integer
        Object number in the catalog.
    ra : float
        Right ascension of the object (hours).
    dec : float
        Declination of the object (degrees).
    pm_ra : float (optional)
        Proper motion in right ascension (milliarcseconds/year).
    pm_dec : float (optional)
        Proper motion in declination (milliarcseconds/year).
    parallax : float (optional)
        Parallax (milliarcseconds).
    rad_vel : float (optional)
        Radial velocity (kilometers/second).

    Returns
    -------
    star : CatEntry
        Instance of CatEntry type object containing the input data.

    Notes
    -----
    .. [N1] This function is equivalent to calling the object with
        arguments,e.g., CatEntry(star_name, catalog, star_num, ...);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """

    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err.format(name='ra'))
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err.format(name='dec'))
    if parallax < 0.0:
        raise ValueError(_neg_err.format(name='parallax'))

    _make_cat_entry = novaslib.make_cat_entry
    _make_cat_entry.argtypes = (ctypes.c_char*51, ctypes.c_char*4,
                                ctypes.c_long, ctypes.c_double,
                                ctypes.c_double, ctypes.c_double,
                                ctypes.c_double, ctypes.c_double,
                                ctypes.c_double, ctypes.POINTER(CatEntry))
    _make_cat_entry.restype = ctypes.c_short
    _make_cat_entry.errcheck = _check_c_errors
    _make_cat_entry.c_errors = {
        1: (ValueError, "from C function 'make_cat_entry': length of 'star_name' out of bounds."),
        2: (ValueError, "from C function 'make_cat_entry': length of 'catalog' out of bounds.")
        }

    star = CatEntry()

    _make_cat_entry(ctypes.create_string_buffer(star_name.encode(), 51),
                    ctypes.create_string_buffer(catalog.encode(), 4), star_num,
                    ra, dec, pm_ra, pm_dec, parallax, rad_vel,
                    ctypes.byref(star))

    return star


def make_object(type, number, name, star_data):
    """
    Makes an instance of Object--specifying a celestial object--based on
    the input parameters.

    Parameters
    ----------
    type : integer
        Type of object
            = 0 ... major planet, Sun, or Moon
            = 1 ... minor planet
            = 2 ... object located outside the solar system
                (i.e. star, galaxy, nebula, etc.)
    number : integer
        For 'type' = 0: Mercury = 1, ..., Pluto = 9, Sun = 10,
            Moon = 11
        For 'type' = 1: minor planet number
        For 'type' = 2: set to 0 (zero)
    name : string
        Name of the object. Max length determined by
        SIZE_OF_OBJ_NAME defined in ``novas.h``; defaults to 50.
    star_data : CatEntry (or ``None`` if 'type' = 0 or 'type' = 1;
        see [N1]_)
        Instance of CatEntry type object containing basic
        astrometric data for any celestial object located outside
        the solar system; the catalog data for a star.

    Returns
    -------
    object : Object
        Instance of Object type object containing the object
        definition.

    Notes
    -----
    .. [N1] If 'type' = 0 or 'type' = 1, ``None`` may be given for
        'star_data'; in this case, a dummy star CatEntry will be
        constructed for 'star_data' automatically.
    .. [N2] This function is equivalent to calling the object with
        arguments, e.g., Object(type, number, name, star_data);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """



    if type not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='type',
                                            allowed=[0, 1, 2]))
    if number not in list(range(0, 12)):
        raise ValueError(_option_err.format(name='number',
                                            allowed=list(range(0, 12))))

    _make_object = novaslib.make_object
    _make_object.argtypes = (ctypes.c_short, ctypes.c_short, ctypes.c_char*51,
                             ctypes.POINTER(CatEntry), ctypes.POINTER(Object))
    _make_object.restype = ctypes.c_short
    _make_object.errcheck = _check_c_errors
    _make_object.c_errors = {
        1: (ValueError, "from C function 'make_object': invalid value of 'type'"),
        2: (ValueError, "from C function 'make_object': 'number' out of range"),
        3: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (object name)."),
        4: (InitializationError, "from C function 'make_object': Initialization of 'cel_obj' failed (catalog name)."),
        5: (ValueError, "from C function 'make_object': 'name' is out of string bounds.")
        }

    if not star_data:
        star_data = CatEntry()

    cel_obj = Object()

    _make_object(type, number, ctypes.create_string_buffer(name.encode(), 51),
                 ctypes.byref(star_data), ctypes.byref(cel_obj))

    return cel_obj


def make_observer(location, obs_surface, obs_space):
    """
    Makes an instance of Observer, specifying the location of the
    observer.

    Parameters
    ----------
    location : {0, 1, 2}
        = 0 ... observer at geocenter
        = 1 ... observer on surface of Earth
        = 2 ... observer on near-earth spacecraft
    obs_surface : OnSurface (or ``None`` if 'location' != 1;
            see [N1]_)
        Instance of OnSurface type object containing data for an
        observer's location on the surface of the Earth; used when
        'location' = 1.
    obs_space : InSpace (or ``None`` if 'location' != 2; see [N1]_)
        Instance of InSpace type object containing data for an
        observer's location on a near-earth spacecraft; used when
        'location' = 2.

    Returns
    -------
    observer : Observer
        Instance of Observer type object specifying the location of
        the observer.

    Notes
    -----
    .. [N1] If 'location' = 0 or 'location' = 2, ``None`` may be
        given for 'obs_surface'; likewise, if 'location' = 0 or
        'location' = 1, ``None`` may be given for 'in_space'. This
        also means that if 'location' = 0, ``None`` may be given for
        both 'obs_surface' and 'in_space'.
    .. [N2] This function is equivalent to calling the object with
        arguments, e.g., Observer(obs_surface, in_space, location);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """

    if location not in [0, 1, 2]:
        raise ValueError(_option_err.format(name='location',
                                            allowed=[0, 1, 2]))

    _make_observer = novaslib.make_observer
    _make_observer.argtypes = (ctypes.c_short, ctypes.POINTER(OnSurface),
                               ctypes.POINTER(InSpace),
                               ctypes.POINTER(Observer))
    _make_observer.restype = ctypes.c_short
    _make_observer.errcheck = _check_c_errors
    _make_observer.c_errors = {
        1: (ValueError, "from C function 'make_observer': input value of 'where' is out-of-range.")
        }

    obs = Observer()

    _make_observer(location, ctypes.byref(obs_surface), ctypes.byref(obs_space),
                   ctypes.byref(obs))

    return obs


def make_observer_at_geocenter():
    """
    Makes an instance of Observer, specifying an observer at the
    geocenter.

    Parameters
    ----------
    None

    Returns
    -------
    observer : Observer
        Instance of Observer type object specifying the location of
        the observer at the geocenter.

    Notes
    -----
    .. [N1] This function is equivalent to calling the object with
        arguments, e.g., Observer(0, None, None); this function
        exists for the purpose of direct compatibility with NOVAS C.

    """

    _make_observer_at_geocenter = novaslib.make_observer_at_geocenter
    _make_observer_at_geocenter.argtypes = (ctypes.POINTER(Observer),)
    _make_observer_at_geocenter.restype = None

    obs_at_geocenter = Observer()

    _make_observer_at_geocenter(ctypes.byref(obs_at_geocenter))

    return obs_at_geocenter


def make_observer_on_surface(latitude, longitude, height, temperature,
                             pressure):
    """
    Makes an instance of Observer, specifying the location of and
    weather for an observer on the surface of the Earth.

    Parameters
    ----------
    latitude : float
        Geodetic (ITRS) latitude in degrees; north positive.
    longitude : float
        Geodetic (ITRS) longitude in degrees; east positive.
    height : float
        Height of the observer (meters).
    temperature : float
        Temperature (degrees Celcius).
    pressure : float
        Atmospheric pressure (millibars).

    Returns
    -------
    observer : Observer
        Instance of Observer type object containing the location of
        and weather for an observer on the surface of the Earth.

    Notes
    -----
    .. [N1] This function is equivalent to calling the object with
        arguments, e.g., Observer(1, obs_surface, None);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """

    if not -90.0 <= latitude <= 90.0:
        raise ValueError(_elev_range_err.format(name='latitude'))
    if not -180.0 <= longitude <= 180.0:
        raise ValueError(_az180_range_err.format(name='longitude'))

    _make_observer_on_surface = novaslib.make_observer_on_surface
    _make_observer_on_surface.argtypes = (ctypes.c_double, ctypes.c_double,
                                          ctypes.c_double, ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.POINTER(Observer))
    _make_observer_on_surface.restype = None

    obs_on_surface = Observer()

    _make_observer_on_surface(latitude, longitude, height, temperature,
                              pressure, ctypes.byref(obs_on_surface))

    return obs_on_surface


def make_observer_in_space(position, velocity):
    """
    Makes an instance of Observer, specifying the position and velocity
    of observer situated on a near-Earth spacecraft.

    Parameters
    ----------
    position : tuple of floats, of length 3
        Geocentric position vector (x, y, z) in km.
    velocity : tuple of floats, of length 3
        Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.

    Returns
    -------
    observer : Observer
        Instance of Observer type object containing the position and
        velocity of an observer situated on a near-Earth spacecraft.

    Notes
    -----
    .. [N1] Both input vector tuples are with respect to true
        equator and equinox of date.
    .. [N2] This function is equivalent to calling the object with
        arguments, e.g., Observer(2, None, obs_space);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """

    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if len(velocity) != 3:
        raise ValueError(_vector_len_err.format(name='velocity'))

    _make_observer_in_space = novaslib.make_observer_in_space
    _make_observer_in_space.argtypes = (ctypes.c_double*3, ctypes.c_double*3,
                                        ctypes.POINTER(Observer))
    _make_observer_in_space.restype = None

    obs_in_space = Observer()

    _make_observer_in_space((ctypes.c_double*3)(*position),
                            (ctypes.c_double*3)(*velocity),
                            ctypes.byref(obs_in_space))

    return obs_in_space


def make_on_surface(latitude, longitude, height, temperature, pressure):
    """
    Makes an instance of OnSurface, specifying the location of and
    weather for an observer on the surface of the Earth.

    Parameters
    ----------
    latitude : float
        Geodetic (ITRS) latitude in degrees; north positive.
    longitude : float
        Geodetic (ITRS) longitude in degrees; east positive.
    height : float
        Height of the observer (meters).
    temperature : float
        Temperature (degrees Celcius).
    pressure : float
        Atmospheric pressure (millibars).

    Returns
    -------
    location : OnSurface
        Instance of OnSurface type object containing the location of
        and weather for an observer on the surface of the Earth.

    Notes
    -----
    .. [N1] This function is equivalent to calling the object with
        arguments, e.g., OnSurface(latitude, longitude, height,
        temperature, pressure); this function exists for the purpose
        of direct compatibility with NOVAS C.

    """

    if not -90.0 <= latitude <= 90.0:
        raise ValueError(_elev_range_err.format(name='latitude'))
    if not -180.0 <= longitude <= 180.0:
        raise ValueError(_az180_range_err.format(name='longitude'))

    _make_on_surface = novaslib.make_on_surface
    _make_on_surface.argtypes = (ctypes.c_double, ctypes.c_double,
                                 ctypes.c_double, ctypes.c_double,
                                 ctypes.c_double, ctypes.POINTER(OnSurface))
    _make_on_surface.restype = None

    obs_surface = OnSurface()

    _make_on_surface(latitude, longitude, height, temperature, pressure,
                     ctypes.byref(obs_surface))

    return obs_surface


def make_in_space(position, velocity):
    """
    Makes an instance of InSpace, specifying the position and velocity
    of an observer situated on a near-Earth spacecraft.

    Parameters
    ----------
    position : tuple of floats, of length 3
        Geocentric position vector (x, y, z) in km.
    velocity : tuple of floats, of length 3
        Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.

    Returns
    -------
    object : InSpace
        Instance of InSpace type object containing the position and
        velocity of an observer situated on a near-Earth spacecraft.

    Notes
    -----
    .. [N1] Both input vector tuples are with respect to true
        equator and equinox of date.
    .. [N2] This function is equivalent to calling the object with
        arguments, e.g., InSpace(sc_pos, sc_vel);
        this function exists for the purpose of direct compatibility
        with NOVAS C.

    """

    if len(position) != 3:
        raise ValueError(_vector_len_err.format(name='position'))
    if len(velocity) != 3:
        raise ValueError(_vector_len_err.format(name='velocity'))

    _make_in_space = novaslib.make_in_space
    _make_in_space.argtypes = (ctypes.c_double*3, ctypes.c_double*3,
                               ctypes.POINTER(InSpace))
    _make_in_space.restype = None

    obs_space = InSpace()

    _make_in_space((ctypes.c_double*3)(*position),
                   (ctypes.c_double*3)(*velocity),
                   ctypes.byref(obs_space))

    return obs_space
