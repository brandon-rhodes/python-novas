# -*- coding: utf-8 -*-

from ctypes import *
from novas import novaslib

class CatEntry(Structure):
    _fields_ = [
        ('starname', c_char*51),
        ('catalog', c_char*4),
        ('starnumber', c_long),
        ('ra', c_double),
        ('dec', c_double),
        ('promora', c_double),
        ('promodec', c_double),
        ('parallax', c_double),
        ('radialvelocity', c_double)
    ]

class Object(Structure):
    _fields_ = [
        ('type', c_short),
        ('number', c_short),
        ('name', c_char*51),
        ('star', CatEntry)
    ]

class OnSurface(Structure):
    _fields_ = [
        ('latitude', c_double),
        ('longitude', c_double),
        ('height', c_double),
        ('temperature', c_double),
        ('pressure', c_double)
    ]

class InSpace(Structure):
    _fields_ = [
        ('_sc_pos_array', c_double*3),
        ('_sc_vel_array', c_double*3)
    ]

    def _get_sc_pos(self):
        return tuple([i for i in self._sc_pos_array])

    def _set_sc_pos(self, sc_pos_tuple):
        self._sc_pos_array = (c_double*3)(*sc_pos_tuple)

    def _get_sc_vel(self):
        return tuple([i for i in self._sc_vel_array])

    def _set_sc_vel(self, sc_vel_tuple):
        self._sc_vel_array = (c_double*3)(*sc_vel_tuple)

    sc_pos = property(_get_sc_pos, _set_sc_pos)
    sc_vel = property(_get_sc_vel, _set_sc_vel)

class Observer(Structure):
    _fields_ = [
        ('where', c_short),
        ('on_surf', OnSurface),
        ('near_earth', InSpace)
    ]

class SkyPos(Structure):
    _fields_ = [
        ('_r_hat_array', c_double*3),
        ('ra', c_double),
        ('dec', c_double),
        ('dis', c_double),
        ('rv', c_double)
    ]

    def _get_r_hat(self):
        return tuple([i for i in self._r_hat_array])

    def _set_r_hat(self, r_hat_tuple):
        self._r_hat_array = (c_double*3)(*r_hat_tuple)

    r_hat = property(_get_r_hat, _set_r_hat)

class RAofCIO(Structure):
    _fields_ = [
        ('jd_tdb', c_double),
        ('ra_cio', c_double)
    ]

_neg_err = "%(name)s must be >= 0.0"
_hour_range_err = "%(name)s must be in the range 0.0 <= %(name)s < 24.0"
_elev_range_err = "%(name)s must be in the range -90.0 <= %(name)s <= 90.0"
_az180_range_err = "%(name)s must be in the range -180.0 <= %(name)s <= 180.0"
_az360_range_err = "%(name)s must be in the range 0.0 <= %(name)s <= 360.0"
_vector_len_err = "%(name)s must be a sequence of length 3"
_option_err = "%(name)s must be in %(allowed)s"
_jd_year_err = "%(name)s must be >= -4712"
_month_range_err = "%(name)s must be 1 <= %(name)s <= 12"
_day_range_err = "%(name)s must be 1 <= %(name)s <= %(daysinmonth)i"

_days_in_month = (None, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

def app_star(jd_tt, star, accuracy=0):
    """
    Computes the apparent place of a star at 'date', given its catalog
    mean place, proper motion, parallax, and radial velocity.
    
    *Parameters*
        jd_tt : float
            TT Julian date for apparent place.
        star : CatEntry
            Instance of CatEntry type object containing catalog data for
            the object in the ICRS.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            position.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (ra, dec) : tuple of floats
            Apparent (right ascension in hours, declination in degrees),
            referred to true equator and equinox of date 'jd_tt'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C61.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _app_star = novaslib.app_star
    _app_star.argtypes = (c_double, POINTER(CatEntry), c_short,
                          POINTER(c_double), POINTER(c_double))
    _app_star.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _app_star(jd_tt, byref(star), accuracy,
                             byref(ra), byref(dec))

    return (ra.value, dec.value)

def virtual_star(jd_tt, star, accuracy=0):
    """
    Computes the virtual place of a star at 'date', given its catalog
    mean place, proper motion, parallax, and radial velocity.
    
    *Parameters*
        jd_tt : float
            TT Julian date for virtual place.
        star : CatEntry
            Instance of CatEntry type object containing catalog data for
            the object in the ICRS.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            position.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (ra, dec) : tuple of floats
            Virtual (right ascension in hours, declination in degrees),
            referred to the GCRS.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C64.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _virtual_star = novaslib.virtual_star
    _virtual_star.argtypes = (c_double, POINTER(CatEntry), c_short,
                              POINTER(c_double), POINTER(c_double))
    _virtual_star.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _virtual_star(jd_tt, byref(star), accuracy,
                                 byref(ra), byref(dec))

    return (ra.value, dec.value)

def astro_star(jd_tt, star, accuracy=0):
    """
    Computes the astrometric place of a star at 'date', given its 
    catalog mean place, proper motion, parallax, and radial velocity.
    
    *Parameters*
        jd_tt : float
            TT Julian date for astrometric place.
        star : CatEntry
            Instance of CatEntry type object containing catalog data for
            the object in the ICRS.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            position.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (ra, dec) : tuple of floats
            Astrometric (right ascension in hours, declination in
            degrees), referred to the ICRS without light deflection or
            aberration.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C66.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _astro_star = novaslib.astro_star
    _astro_star.argtypes = (c_double, POINTER(CatEntry), c_short,
                            POINTER(c_double), POINTER(c_double))
    _astro_star.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _astro_star(jd_tt, byref(star), c_short,
                               byref(ra), byref(dec))

    return (ra.value, dec.value)

def app_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the apparent place of a planet or other solar system body.
    
    *Parameters*
        jd_tt : float
            TT Julian date for apparent place.
        ss_body : Object
            Instance of Object type object containing the body
            designation for the solar system body.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            position.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (ra, dec, dis) : tuple of floats
            Apparent (right ascension in hours, declination in degrees,
            ...), referred to true equator and equinox of date, and true
            (..., ..., distance in AU) from Earth to solar system body
            at 'jd_tt'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C67.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _app_planet = novaslib.app_planet
    _app_planet.argtypes = (c_double, POINTER(Object), c_short,
                            POINTER(c_double), POINTER(c_double),
                            POINTER(c_double))
    _app_planet.restype = c_short

    ra = c_double()
    dec = c_double()
    dis = c_double()

    return_value = _app_planet(jd_tt, byref(ss_body), accuracy,
                               byref(ra), byref(dec), byref(dis))

    return (ra.value, dec.value, dis.value)

def virtual_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the virtual place of a planet or other solar system body.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec, dis) : tuple of floats
            Virtual (right ascension in hours, declination in degrees,
            ...), referred to the GCRS, and true (..., ..., distance in
            AU) from Earth to solar system body.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C70.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _virtual_planet = novaslib.virtual_planet
    _virtual_planet.argtypes = (c_double, POINTER(Object), c_short,
                                POINTER(c_double), POINTER(c_double),
                                POINTER(c_double))
    _virtual_planet.restype = c_short

    ra = c_double()
    dec = c_double()
    dis = c_double()

    return_value = _virtual_planet(jd_tt, byref(ss_body), accuracy,
                                   byref(ra), byref(dec), byref(dis))

    return (ra.value, dec.value, dis.value)

def astro_planet(jd_tt, ss_body, accuracy=0):
    """
    Compute the astrometric place of a planet or other solar system
    body.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec, dis) : tuple of floats
            Astrometric (right ascension in hours, declination in
            degrees, ...), referred to the ICRS without light deflection
            or aberration, and true (..., ..., distance in AU) from
            Earth to solar system body.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C73.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _astro_planet = novaslib.astro_planet
    _astro_planet.argtypes = (c_double, POINTER(Object), c_short,
                              POINTER(c_double), POINTER(c_double),
                              POINTER(c_double))
    _astro_planet.restype = c_short

    ra = c_double()
    dec = c_double()
    dis = c_double()

    return_value = _astro_planet(jd_tt, byref(ss_body), accuracy,
                                 byref(ra), byref(dec), byref(dis))

    return (ra.value, dec.value, dis.value)

def topo_star(jd_tt, delta_t, star, position, accuracy=0):
    """
    Computes the topocentric place of a star at 'date', given its
    catalog mean place, proper motion, parallax, and radial velocity.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec) : tuple of floats
            Topocentric (right ascension in hours, declination in
            degrees), referred to true equator and equinox of date
            'jd_tt'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C62-C63.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                          'allowed': [0, 1]})

    _topo_star = novaslib.topo_star
    _topo_star.argtypes = (c_double, c_double, POINTER(CatEntry),
                           POINTER(OnSurface), c_short, POINTER(c_double),
                           POINTER(c_double))
    _topo_star.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _topo_star(jd_tt, delta_t, byref(star), byref(position),
                              accuracy, byref(ra), byref(dec))

    return (ra.value, dec.value)

def local_star(jd_tt, delta_t, star, position, accuracy=0):
    """
    Computes the local place of a star at date 'date', given its
    catalog mean place, proper motion, parallax, and radial velocity.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec) : tuple of floats
            Local (right ascension in hours, declination in degrees),
            referred to the 'local GCRS'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C65.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _local_star = novaslib.local_star
    _local_star.argtypes = (c_double, c_double, POINTER(CatEntry),
                            POINTER(OnSurface), c_short, POINTER(c_double),
                            POINTER(c_double))
    _local_star.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _local_star(jd_tt, delta_t, byref(star), byref(position),
                               accuracy, byref(ra), byref(dec))

    return (ra.value, dec.value)

def topo_planet(jd_tt, delta_t, ss_body, position, accuracy=0):
    """
    Computes the topocentric place of a solar system body.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec, dis) : tuple of floats
            Topocentric (right ascension in hours, declination in
            degrees, ...), referred to true equator and equinox of date,
            and true (..., ..., distance in AU) from Earth to solar
            system body at 'jd_tt'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C68-C69.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _topo_planet = novaslib.topo_planet
    _topo_planet.argtypes = (c_double, POINTER(Object), c_double,
                             POINTER(OnSurface), c_short, POINTER(c_double),
                             POINTER(c_double), POINTER(c_double))
    _topo_planet.restype = c_short

    ra = c_double()
    dec = c_double()
    dis = c_double()

    return_value = _topo_planet(jd_tt, byref(ss_body), delta_t,
                                byref(position), accuracy, byref(ra),
                                byref(dec), byref(dis))

    return (ra.value, dec.value, dis.value)

def local_planet(jd_tt, delta_t, ss_body, position, accuracy=0):
    """
    Computes the local place of a solar system body.
    
    *Parameters*
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
    
    *Returns*
        (ra, dec, dis) : tuple of floats
            Local (right ascension in hours, declination in degrees,
            ...), referred to the 'local GCRS', and true (..., ...,
            distance in AU) from Earth to solar system body at 'jd_tt'.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C72-C72.
        .. [R2] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _local_planet = novaslib.local_planet
    _local_planet.argtypes = (c_double, POINTER(Object), c_double,
                              POINTER(OnSurface), c_short, POINTER(c_double),
                              POINTER(c_double), POINTER(c_double))
    _local_planet.restype = c_short

    ra = c_double()
    dec = c_double()
    dis = c_double()

    return_value = _local_planet(jd_tt, byref(ss_body), delta_t,
                                 byref(position), accuracy, byref(ra),
                                 byref(dec), byref(dis))

    return (ra.value, dec.value, dis.value)

def mean_star(jd_tt, ra, dec, accuracy=0):
    """
    Computes the ICRS position of a star, given its apparent place
    at 'date'. Proper motion, parallax, and radial velocity are assumed
    to be zero.
    
    *Parameters*
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
    
    *Returns*
        (ira, idec) : tuple of floats
            ICRS (right ascension in hours, declination in degrees).
    
    *References*
        .. [R1] Explanatory Supplement to the Astronomical Almanac
            (1992), Chapter 3.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err % {'name': 'ra'})
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'dec'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _mean_star = novaslib.mean_star
    _mean_star.argtypes = (c_double, c_double, c_double, c_short,
                           POINTER(c_double), POINTER(c_double))
    _mean_star.restype = c_short

    ira = c_double()
    idec = c_double()

    return_value = _mean_star(jd_tt, ra, dec, accuracy,
                              byref(ira), byref(idec))

    return (ira.value, idec.value)

def place(jd_tt, delta_t, cel_object, location, coord_sys, accuracy=0):
    """
    This function computes the apparent direction of a star or solar
    system body at a specified time and in a specified coordinate
    system.
    
    *Parameters*
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
    
    *Returns*
        output : SkyPos
            Instance of SkyPos type object specifying object's place on
            the sky at time 'jd_tt', with respect to the specified
            output coordinate system.
    
    *Notes*
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
    
    *References*
        .. [R1] Kaplan, G. et al. (1989), Astronomical Journal 97,
            1197-1210.
        .. [R2] Klioner, S. (2003), Astronomical Journal 125, 1580-1597.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if coord_sys not in [0, 1, 2, 3]:
        raise ValueError(_option_err % {'name': 'coord_sys',
                                        'allowed': [0, 1, 2, 3]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _place = novaslib.place
    _place.argtypes = (c_double, POINTER(Object), POINTER(Observer), c_double,
                       c_short, c_short, POINTER(SkyPos))
    _place.restype = c_short

    output = SkyPos()

    return_value = _place(jd_tt, byref(cel_object), byref(location), delta_t,
                          coord_sys, accuracy, byref(output))

    return output

def equ2gal(rai, deci):
    """
    Convert ICRS equatorial position to galactic position
    
    *Parameters*
        ra : float
            ICRS right ascension in hours.
        dec : float
            ICRS declination in degrees.
    
    *Returns*
        (glon, glat) : tuple of floats
            Galactic (longitude, latitude) in degrees.
    
    *References*
        .. [R1] Hipparcos and Tycho Catalogues, Vol. 1, Section 1.5.3.
    
    """

    if not 0.0 <= rai < 24.0:
        raise ValueError(_hour_range_err % {'name': 'rai'})
    if not -90.0 <= deci <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'deci'})

    _equ2gal = novaslib.equ2gal
    _equ2gal.argtypes = (c_double, c_double, POINTER(c_double),
                         POINTER(c_double))
    _equ2gal.restype = None

    glon = c_double()
    glat = c_double()

    _equ2gal(rai, deci, byref(glon), byref(glat))

    return (glon.value, glat.value)

def equ2ecl(jd_tt, ra, dec, coord_sys=2, accuracy=0):
    """
    Convert equatorial position to ecliptic position
    
    *Parameters*
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
    
    *Returns*
        (elon, elat) : tuple of floats
            Ecliptic (longitude, latitude) in degrees, referred to
            specified ecliptic and equinox of date.
    
    *Notes*
        .. [N1] To convert ICRS RA and dec to ecliptic coordinates (mean
            ecliptic and equinox of J2000.0), set 'system' = 2; the
            value of 'jd_tt' can be set to anything, since J2000.0 is
            assumed. Except for the input to this case, all input
            coordinates are dynamical.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err % {'name': 'ra'})
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'dec'})
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'coord_sys',
                                        'allowed': [0, 1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _equ2ecl = novaslib.equ2ecl
    _equ2ecl.argtypes = (c_double, c_short, c_short, c_double, c_double,
                         POINTER(c_double), POINTER(c_double))
    _equ2ecl.restype = c_short

    elon = c_double()
    elat = c_double()

    return_value = _equ2ecl(jd_tt, coord_sys, accuracy, ra, dec,
                            byref(elon), byref(elat))

    return elon.value, elat.value

def equ2ecl_vec(jd_tt, position, coord_sys=2, accuracy=0):
    """
    Convert an equatorial position vector to an ecliptic position
    vector.
    
    *Parameters*
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
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, referred to specified ecliptic and equinox
            of date.
    
    *Notes*
        .. [N1] To convert an ICRS vector to an ecliptic vector (mean
            ecliptic and equinox of J2000.0 only), set 'system' = 2; the
            value of 'jd_tt' can be set to anything, since J2000.0 is
            assumed. Except for the input to this case, all vectors are
            assumed to be with respect to a dynamical system.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'coord_sys',
                                        'allowed': [0, 1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _equ2ecl_vec = novaslib.equ2ecl_vec
    _equ2ecl_vec.argtypes = (c_double, c_short, c_short, POINTER(c_double*3),
                             POINTER(c_double*3))
    _equ2ecl_vec.restype = c_short

    pos2 = (c_double*3)()

    return_value = _equ2ecl_vec(jd_tt, coord_sys, accuracy,
                                byref((c_double*3)(*position)), byref(pos2))

    return tuple([i for i in pos2])

def ecl2equ_vec(jd_tt, position, coord_sys=2, accuracy=0):
    """
    Converts an ecliptic position vector to an equatorial position
    vector.
    
    *Parameters*
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
    
    *Returns*
        position2 : tuple of floats, of length 3
            Position vector, referred to specified ecliptic and equinox
            of date.
    
    *Notes*
        .. [N1] To convert an ecliptic vector (mean ecliptic and equinox
            of J2000.0 only) to an ICRS vector, set 'system' = 2; the
            value of 'jd_tt' can be set to anything, since J2000.0 is
            assumed. Except for the output from this case, all vectors
            are assumed to be with respect to a dynamical system.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'coord_sys',
                                        'allowed': [0, 1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _ecl2equ_vec = novaslib.ecl2equ_vec
    _ecl2equ_vec.argtypes = (c_double, c_short, c_short, POINTER(c_double*3),
                             POINTER(c_double*3))
    _ecl2equ_vec.restype = c_short

    pos2 = (c_double*3)()

    return_value = _ecl2equ_vec(jd_tt, coord_sys, accuracy,
                                byref((c_double*3)(*position)), byref(pos2))

    return tuple([i for i in pos2])

def equ2hor(jd_ut1, delta_t, xp, yp, location, ra, dec, ref_option=0,
            accuracy=0):
    """
    This function transforms topocentric right ascension and declination
    to zenith distance and azimuth. It uses a method that properly
    accounts for polar motion, which is significant at the sub-arcsecond
    level. This function can also adjust coordinates for atmospheric
    refraction.
    
    *Parameters*
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
            Topocentric declination of object of interest, in hours,
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
    
    *Returns*
        (zd, az) : tuple of floats
            Topocentric (zenith distance, azimuth) in degrees. Zenith
            distance is affected by refraction if 'ref_option', is
            non-zero. Azimuth is measured east from north.
        (rar, decr) : tuple of floats
            Topocentric (right ascension in hours, declination in
            degrees) of object of interest, referred to true equator and
            equinox of date, affected by refraction if 'ref_option' is
            non-zero.
    
    *Notes*
        .. [N1] 'xp' and 'yp' can be set to zero if sub-arcsecond
            accuracy is not needed. 'ra' and 'dec' can be obtained from
            functions 'topo_star' or 'topo_planet'.
        .. [N2] The directions 'zd'= 0 (zenith) and 'az'= 0 (North) are
            here considered fixed in the terrestrial system.
            Specifically, the zenith is along the geodetic normal, and
            North is toward the ITRS pole.
        .. [N3] If 'ref_option'= 0, then 'rar'='ra' and 'decr'='dec'.
    
    *References*
        .. [R1] Kaplan, G. (2008). USNO/AA Technical Note of 28 Apr
            2008, "Refraction as a Vector."
    
    """

    if jd_ut1 < 0.0: raise ValueError(_neg_err % {'name': 'jd_ut1'})
    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err % {'name': 'ra'})
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'dec'})
    if ref_option not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'ref_option',
                                        'allowed': [0, 1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _equ2hor = novaslib.equ2hor
    _equ2hor.argtypes = (c_double, c_double, c_short, c_double, c_double,
                         POINTER(OnSurface), c_double, c_double, c_short,
                         POINTER(c_double), POINTER(c_double),
                         POINTER(c_double), POINTER(c_double))
    _equ2hor.restype = None

    zd = c_double()
    az = c_double()
    rar = c_double()
    decr = c_double()

    _equ2hor(jd_ut1, delta_t, accuracy, xp, yp, byref(location), ra, dec,
             ref_option, byref(zd), byref(az), byref(rar), byref(decr))

    return (zd.value, az.value), (rar.value, decr.value)

def gcrs2equ(jd_tt, rag, decg, coord_sys=1, accuracy=0):
    """
    Convert GCRS right ascension and declination to coordinates with
    respect to the equator of date (mean or true). For coordinates with
    respect to the true equator of date, the origin of right ascension
    can be either the true equinox or the celestial intermediate origin
    (CIO).
    
    *Parameters*
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
    
    *Returns*
        (ra, dec) : tuple of floats
            (Right ascension in hours, declination in degrees), referred
            to specified equator (and right ascension origin, for 'ra')
            of date.
    
    *Notes*
        .. [N1] Set input value of 'accuracy' equal to any short int if
            'coord_sys' equals 0 or 1. It is not used in these cases.
        .. [N2] This function only supports the CIO-based method.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if not 0.0 <= rag < 24.0:
        raise ValueError(_hour_range_err % {'name': 'rag'})
    if not -90.0 <= decg <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'decg'})
    if coord_sys not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'coord_sys',
                                        'allowed': [0, 1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _gcrs2equ = novaslib.gcrs2equ
    _gcrs2equ.argtypes = (c_double, c_short, c_short, c_double, c_double,
                          POINTER(c_double), POINTER(c_double))
    _gcrs2equ.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _gcrs2equ(jd_tt, coord_sys, accuracy, rag, decg,
                             byref(ra), byref(dec))

    return (ra.value, dec.value)

def sidereal_time(jd_high, jd_low, delta_t, gst_type=1, method=0, accuracy=0):
    """
    Computes the Greenwich sidereal time, either mean or apparent, at
    Julian date 'jd_high' + 'jd_low'.
    
    *Parameters*
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
    
    *Returns*
        gst : float
            Greenwich (mean or apparent) sidereal time, in hours.
    
    *Notes*
        .. [N1] The Julian date may be split at any point, but for
            highest precision, set 'jd_high' to be the integral part of
            the Julian date, and set 'jd_low' to be the fractional part.
    
    *References*
        .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if jd_high + jd_low < 0.0:
        raise ValueError(_neg_err % {'name': 'jd_high+jd_low'})
    if gst_type not in [0, 1]:
        raise ValueError(_option_err % {'name': 'gst_type',
                                        'allowed': [0, 1]})
    if method not in [0, 1]:
        raise ValueError(_option_err % {'name': 'method',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _sidereal_time = novaslib.sidereal_time
    _sidereal_time.argtypes = (c_double, c_double, c_double, c_short, c_short,
                               c_short, POINTER(c_double))
    _sidereal_time.restype = c_short

    gst = c_double()

    return_value = _sidereal_time(jd_high, jd_low, delta_t, gst_type, method,
                                  accuracy, byref(gst))

    return gst.value

def era(jd_ut1_high, jd_ut1_low=0.0):
    """
    Compute the Earth Rotation Angle (theta) for a given UT1 Julian
    date. The expression used is taken from the note to IAU Resolution
    B1.8 of 2000[R1]_.
    
    *Parameters*
        jd_ut1_high : float
            High-order part of UT1 Julian date.
        jd_ut1_low : float (optional)
            Low-order part of UT1 Julian date.
    
    *Returns*
        theta : float
            The Earth Rotation Angle in degrees.
    
    *Notes*
        .. [N1] The algorithm used here is equivalent to the canonical
            theta = 0.7790572732640 + 1.00273781191135448 * t,
            where t is the time in days from J2000 (t = jd_ut1_high +
            jd_ut1_low - T0), but it avoids many two-PI 'wraps' that
            decrease precision (adopted from SOFA Fortran routine
            iau_era00; see also expression at top of page 35 of IERS
            Conventions (1996)).
    
    *References*
        .. [R1] IAU Resolution B1.8, adopted at the 2000 IAU General
            Assembly, Manchester, UK.
        .. [R2] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err % {'name': 'jd_ut1_high+jd_ut1_low'})

    _era = novaslib.era
    _era.argtypes = (c_double, c_double)
    _era.restype = c_double

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
    
    *Parameters*
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
    
    *Returns*
        vec2 : tuple of floats, of length 3
            Position vector, geocentric equatorial rectangular
            coordinates, referred to GCRS axes (celestial system) or
            with respect to the equator and equinox of date, depending
            on 'option'.
    
    *Notes*
        .. [N1] 'xp' = 'yp' = 0 means no polar motion transformation.
        .. [N2] The 'option' flag only works for the equinox-based
            method.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C53-C54.
        .. [R2] Kaplan, G. H. (2003), 'Another Look at Non-Rotating
            Origins', Proceedings of IAU XXV Joint Discussion 16.
    
    """

    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err % {'name': 'jd_ut1_high+jd_ut1_low'})
    if len(vec1) != 3: raise ValueError(_vector_len_err % {'name': 'vec1'})
    if method not in [0, 1]:
        raise ValueError(_option_err % {'name': 'method',
                                        'allowed': [0, 1]})
    if option not in [0, 1]:
        raise ValueError(_option_err % {'name': 'option',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})


    _ter2cel = novaslib.ter2cel
    _ter2cel.argtypes = (c_double, c_double, c_double, c_short, c_short,
                         c_short, c_double, c_double, POINTER(c_double*3),
                         POINTER(c_double*3))
    _ter2cel.restype = c_short

    vec2 = (c_double*3)()

    return_value = _ter2cel(jd_ut1_high, jd_ut1_low, delta_t, method, accuracy,
                            option, xp, yp, byref((c_double*3)(*vec1)),
                            byref(vec2))

    return tuple([i for i in vec2])

def cel2ter(jd_ut1_high, jd_ut1_low, delta_t, xp, yp, vec1, method=0, option=0,
            accuracy=0):
    """
    Rotates a vector from the celestial to the terrestrial system.
    Specifically, transforms a vector in the GCRS (a local space-fixed
    system) to the ITRS (a rotating earth-fixed system) by applying
    rotations for the GCRS-to-dynamical frame tie, precession, nutation,
    Earth rotation, and polar motion.
    
    *Parameters*
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
    
    *Returns*
        vec2 : tuple of floats, of length 3
        Position vector, geocentric equatorial rectangular coordinates,
        referred to ITRS axes (terrestrial system) in the normal case
        where 'option' = 0.
    
    *Notes*
        .. [N1] 'xp' = 'yp' = 0 means no polar motion transformation.
    .. [N2] The 'option' flag only works for the equinox-based method.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C54.
        .. [R2] Kaplan, G. H. (2003), 'Another Look at Non-Rotating
            Origins', Proceedings of IAU XXV Joint Discussion 16.
    
    """

    if jd_ut1_high + jd_ut1_low < 0.0:
        raise ValueError(_neg_err % {'name': 'jd_ut1_high+jd_ut1_low'})
    if len(vec1) != 3: raise ValueError(_vector_len_err % {'name': 'vec1'})
    if method not in [0, 1]:
        raise ValueError(_option_err % {'name': 'method',
                                        'allowed': [0, 1]})
    if option not in [0, 1]:
        raise ValueError(_option_err % {'name': 'option',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _cel2ter = novaslib.cel2ter
    _cel2ter.argtypes = (c_double, c_double, c_double, c_short, c_short,
                         c_short, c_double, c_double, POINTER(c_double*3),
                         POINTER(c_double*3))
    _cel2ter.restype = c_short

    vec2 = (c_double*3)()

    return_value = _cel2ter(jd_ut1_high, jd_ut1_low, delta_t, method, accuracy,
                            option, xp, yp, byref((c_double*3)(*vec1)),
                            byref(vec2))

    return tuple([i for i in vec2])

def spin(angle, pos1):
    """
    This function transforms a vector from one coordinate system to
    another with same origin and axes rotated about the z-axis.
    
    *Parameters*
        angle : float
            Angle of coordinate system rotation, positive
            counterclockwise when viewed from +z, in degrees.
        position : tuple of floats, of length 3
            Position vector.
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector expressed in new coordinate system rotated
            about z by 'angle'.
    
    """

    if not 0.0 <= angle < 360.0:
        raise ValueError(_az360_range_err % {'name': 'angle'})
    if len(pos1) != 3: raise ValueError(_vector_len_err % {'name': 'pos1'})

    _spin = novaslib.spin
    _spin.argtypes = (c_double, POINTER(c_double*3), POINTER(c_double*3))
    _spin.restype = None

    pos2 = (c_double*3)()

    _spin(angle, byref((c_double*3)(*pos1)), byref(pos2))

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
    
    *Parameters*
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
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, geocentric equatorial rectangular
            coordinates, referred to true equator and TIO.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C23, C110, C111, C113, C114.
        .. [R2] Lambert & Bizouard (2002), Astronomy and Astrophysics
            394, 317-321.
    
    """

    if tjd < 0.0: raise ValueError(_neg_err % {'name': 'tjd'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if direction not in [0, 1]:
        raise ValueError(_option_err % {'name': 'direction',
                                        'allowed': [0, 1]})

    _wobble = novaslib.wobble
    _wobble.argtypes = (c_double, c_short, c_double, c_double,
                        POINTER(c_double*3), POINTER(c_double*3))
    _wobble.restype = None

    pos2 = (c_double*3)()

    _wobble(tjd, direction, x, y, byref((c_double*3)(*position)), byref(pos2))

    return tuple([i for i in pos2])

def terra(location, time):
    """
    Computes the position and velocity vectors of a terrestrial observer
    with respect to the center of the Earth.
    
    *Parameters*
        location : OnSurface
            Instance of OnSurface type object containing observer's
            location.
        time : float
            Local apparent sidereal time at reference meridian in hours.
    
    *Returns*
        (position, velocity) : tuple of tuple of floats, of length 3
            Position and velocity vector of observer with respect to
            center of Earth in equatorial rectangular coordinates,
            referred to true equator and equinox of date. Components in
            AU and AU/day.
    
    *Notes*
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
        raise ValueError(_hour_range_err % {'name': 'time'})

    _terra = novaslib.terra
    _terra.argtypes = (POINTER(OnSurface), c_double,
                       POINTER(c_double*3), POINTER(c_double*3))
    _terra.restype = None

    pos = (c_double*3)()
    vel = (c_double*3)()

    _terra(byref(location), time, byref(pos), byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])

def e_tilt(jd_tdb, accuracy=0):
    """
    Computes quantities related to the orientation of the Earth's
    rotation axis at the given Julian day.
    
    *Parameters*
        jd_tdb : float
            TDB Julian date.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            quantities.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (mobl, tobl, ee, dpsi, deps) : tuple of five floats
            mobl = mean obliquity of the ecliptic in degrees.
            tobl = true obliquity of the ecliptic in degrees.
            ee   = equation of the equinoxes in seconds of time.
            dpsi = nutation in longitude in arcseconds.
            deps = nutation in obliquity in arcseconds.
    
    *Notes*
        .. [N1] Values of the celestial pole offsets 'PSI_COR' and
            'EPS_COR' are set using function 'cel_pole', if desired. See
            the docstring for 'cel_pole' for details.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _e_tilt = novaslib.e_tilt
    _e_tilt.argtypes = (c_double, c_short, POINTER(c_double),
                        POINTER(c_double), POINTER(c_double),
                        POINTER(c_double), POINTER(c_double))
    _e_tilt.restype = None

    mobl = c_double()
    tobl = c_double()
    ee = c_double()
    dpsi = c_double()
    deps = c_double()

    _e_tilt(jd_tdb, accuracy, byref(mobl), byref(tobl), byref(ee),
            byref(dpsi), byref(deps))

    return mobl.value, tobl.value, ee.value, dpsi.value, deps.value

def cel_pole(tjd, type, dpole1, dpole2):
    """
    This function allows for the specification of celestial pole
    offsets for high-precision applications. Each set of offsets is
    a correction to the modeled position of the pole for a specific
    date, derived from observations and published by the IERS.
    
    *Parameters*
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
    
    *Notes*
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
    
    *References*
        .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
        .. [R2] Kaplan, G. (2003), USNO/AA Technical Note 2003-03.
    
    """

    if tjd < 0.0: raise ValueError(_neg_err % {'name': 'tjd'})
    if type not in [0, 1]:
        raise ValueError(_option_err % {'name': 'type',
                                        'allowed': [0, 1]})

    _cel_pole = novaslib.celpole
    _cel_pole.argtypes = (c_double, c_short, c_double, c_double)
    _cel_pole.restype = c_short

    return_value = _cel_pole(tjd, type, dpole1, dpole2)

    return None

def ee_ct(jd_tt_high, jd_tt_low, accuracy=0):
    """
    To compute the \"complementary terms\" of the equation of the
    equinoxes.
    
    *Parameters*
        jd_tt_high : float
            High-order (integer) part of TT Julian day.
        jd_tt_low : float
            Low-order (fractional) part of TT Julian day.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the terms.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        comp_terms : float
            Complementary terms, in radians.
    
    *Notes*
        .. [N1] The series used in this function was derived from the
            first reference.  This same series was also adopted for use
            in the IAU's Standards of Fundamental Astronomy (SOFA)
            software (i.e., subroutine eect00.for and function
            eect00.c).
        .. [N2] The low-accuracy series used in this function is a
            simple implementation derived from the first reference, in
            which terms smaller than 2 microarcseconds have been
            omitted.
    
    *References*
        .. [R1] Capitaine, N., Wallace, P.T., and McCarthy, D.D. (2003).
            Astron. & Astrophys. 406, p. 1135-1149. Table 3.
        .. [R2] IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
            (Table 5.2e presented in the printed publication is a
            truncated series. The full series, which is used in NOVAS,
            is available on the IERS Conventions Center website in file
            tab5.2e.txt.) ftp://tai.bipm.org/iers/conv2010/chapter5/
    
    """

    if jd_tt_high + jd_tt_low < 0.0:
        raise ValueError(_neg_err % {'name': 'jd_tt_high+jd_tt_low'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _ee_ct = novaslib.ee_ct
    _ee_ct.argtypes = (c_double, c_double, c_short)
    _ee_ct.restype = c_double

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
    
    *Parameters*
        position : tuple of floats, of length 3
            Position vector in equatorial rectangular coordinates.
        direction : {0, -1}, optional
            = 0 ... ICRS to dynamical (default)
            = -1 ... dynamical to ICRS
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector in equatorial rectangular coordinates.
    
    *Notes*
        .. [N1] For geocentric coordinates, the same transformation is
            used between the dynamical reference system and the GCRS.
    
    *References*
        .. [R1] Hilton, J. and Hohenkerk, C. (2004), Astronomy and
            Astrophysics 413, 765-770, eq. (6) and (8).
        .. [R2] IERS (2003) Conventions, Chapter 5.
    
    """

    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if direction not in [0, -1]:
        raise ValueError(_option_err % {'name': 'direction',
                                        'allowed': [0, -1]})

    _frame_tie = novaslib.frame_tie
    _frame_tie.argtypes = (POINTER(c_double*3), c_short, POINTER(c_double*3))
    _frame_tie.restype = None

    pos2 = (c_double*3)()

    _frame_tie(byref((c_double*3)(*position)), direction, byref(pos2))

    return tuple([i for i in pos2])

def proper_motion(jd_tdb1, position, velocity, jd_tdb2):
    """
    Applies proper motion, including foreshortening effects, to a star's
    position.
    
    *Parameters*
        jd_tdb1 : float
            TDB Julian date of first epoch.
        position : tuple of floats, of length 3
            Position vector at first epoch.
        velocity : tuple of floats, of length 3
            Velocity vector at first epoch.
        jd_tdb2 : float
            TDB Julian date of second epoch.
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector at second epoch.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C16.
    
    """

    if jd_tdb1 < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb1'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if len(velocity) != 3:
        raise ValueError(_vector_len_err % {'name': 'velocity'})
    if jd_tdb2 < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb2'})

    _proper_motion = novaslib.proper_motion
    _proper_motion.argtypes = (c_double, POINTER(c_double*3),
                               POINTER(c_double*3), c_double,
                               POINTER(c_double*3))
    _proper_motion.restype = None

    pos2 = (c_double*3)()

    _proper_motion(jd_tdb1, byref((c_double*3)(*position)),
                   byref((c_double*3)(*velocity)), jd_tdb2, byref(pos2))

    return tuple([i for i in pos2])

def bary2obs(pos, pos_obs):
    """
    Transform the origin from the solar system barycenter to the
    observer (or the geocenter); i.e., correct for parallax
    (annual+geocentric or just annual).
    
    *Parameters*
        pos : tuple of floats, of length 3
            Position vector, referred to origin at solar system
            barycenter, components in AU.
        pos_obs : tuple of floats, of length 3
            Position vector of observer (or the geocenter), with respect
            to origin at solar system barycenter, components in AU.
    
    *Returns*
        (position, lighttime) : tuple of tuple and float
            Position vector, referred to origin at center of mass of the
            Earth, components in AU, and lighttime of object from Earth
            in days.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C17, C109.
    
    """

    if len(pos) != 3: raise ValueError(_vector_len_err % {'name': 'pos'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})

    _bary2obs = novaslib.bary2obs
    _bary2obs.argtypes = (POINTER(c_double*3), POINTER(c_double*3),
                          POINTER(c_double*3), POINTER(c_double))
    _bary2obs.restype = None

    pos2 = (c_double*3)()
    lighttime = c_double()

    _bary2obs(byref((c_double*3)(*pos)), byref((c_double*3)(*pos_obs)),
              byref(pos2), byref(lighttime))

    return tuple([i for i in pos2]), lighttime.value

def geo_posvel(jd_tt, delta_t, observer, accuracy=0):
    """
    Compute the geocentric position and velocity of an observer on the
    surface of the earth or on a near-earth spacecraft. The final
    vectors are expressed in the GCRS.
    
    *Parameters*
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
    
    *Returns*
        (position, velocity) : tuple of tuple of floats, of length 3
            Position vector of observer, with respect to origin at
            geocenter, referred to GCRS axes, components in AU and
            AU/day, respectively.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _geo_posvel = novaslib.geo_posvel
    _geo_posvel.argtypes = (c_double, c_double, c_short, POINTER(Observer),
                            POINTER(c_double*3), POINTER(c_double*3))
    _geo_posvel.restype = c_short

    pos = (c_double*3)()
    vel = (c_double*3)()

    return_value = _geo_posvel(jd_tt, delta_t, accuracy, byref(observer),
                               byref(pos), byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])

def light_time(jd_tdb, ss_object, pos_obs, estimate=0.0, accuracy=0):
    """
    Compute the geocentric position of a solar system body, as antedated
    for light-time.
    
    *Parameters*
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
    
    *Returns*
        tuple of tuple of floats, of length 3, and float
            Position vector of body, with respect to origin at observer
            (or the geocenter), referred to ICRS axes, components in AU
            and final light-time, in days.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})
    if estimate < 0.0: raise ValueError(_neg_err % {'name': 'estimate'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _light_time = novaslib.light_time
    _light_time.argtypes = (c_double, POINTER(Object), POINTER(c_double*3),
                            c_double, c_short, POINTER(c_double*3),
                            POINTER(c_double))
    _light_time.restype = c_short

    pos = (c_double*3)()
    tlight = c_double

    return_value = _light_time(jd_tdb, byref(ss_object),
                               byref((c_double*3)(*pos_obs)), estimate,
                               accuracy, byref(pos), byref(tlight))

    return tuple([i for i in pos]), tlight.value

def d_light(pos_star, pos_obs):
    """
    Computes the difference in light-time, for a star, between the
    barycenter of the solar system and the observer (or the geocenter).
    
    *Parameters*
        pos_star : tuple of floats, of length 3
            Position vector of star, with respect to origin at solar
            system barycenter.
        pos_obs : tuple of floats, of length 3
            Position vector of observer (or the geocenter), with respect
            to origin at solar system barycenter, components in AU.
    
    *Returns*
        diflt : float
            Difference in light time, in the sense star to barycenter
            minus star to earth, in days.
    
    *Notes*
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
        raise ValueError(_vector_len_err % {'name': 'pos_star'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})

    _d_light = novaslib.d_light
    _d_light.argtypes = (POINTER(c_double*3), POINTER(c_double*3))
    _d_light.restype = c_double

    diflt = _d_light(byref((c_double*3)(*pos_star)),
                    byref((c_double*3)(*pos_obs)))

    return diflt

def grav_def(jd_tdb, pos_obj, pos_obs, location, accuracy=0):
    """
    Compute the total gravitational deflection of light for the observed
    object due to the major gravitating bodies in the solar system. This
    function valid for an observed body within the solar system as well
    as for a star.
    
    *Parameters*
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
    
    *Returns*
        positon : tuple of floats, of length 3
            Position vector of observed object, with respect to origin
            at observer (or the geocenter), referred to ICRS axes,
            corrected for gravitational deflection, components in AU.
    
    *Notes*
        .. [N1] If 'accuracy' is set to zero (full accuracy), three
            bodies (Sun, Jupiter, and Saturn) are used in the
            calculation. If the reduced-accuracy option is set, only the
            Sun is used in the calculation. In both cases, if the
            observer is not at the geocenter, the deflection due to the
            Earth is included.
        .. [N2] The number of bodies used at full and reduced accuracy
            can be set by making a change to the code in this function
            as indicated in the comments.
    
    *References*
        .. [R1] Klioner, S. (2003), Astronomical Journal 125, 1580-1597,
            Section 6.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obj'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})
    if location not in [0, 1]:
        raise ValueError(_option_err % {'name': 'location',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _grav_def = novaslib.grav_def
    _grav_def.argtypes = (c_double, c_short, c_short, POINTER(c_double*3),
                          POINTER(c_double*3), POINTER(c_double*3))
    _grav_def.restype = c_short

    pos2 = (c_double*3)()

    return_value = _grav_def(jd_tdb, location, accuracy,
                             byref((c_double*3)(*pos_obj)),
                             byref((c_double*3)(*pos_obs)), byref(pos2))

    return tuple([i for i in pos2])

def grav_vec(pos_obj, pos_obs, pos_body, rmass):
    """
    Correct the position vector for the deflection of light in the
    gravitational field of an arbitrary body. This function valid for an
    observed body within the solar system as well as for a star.
    
    *Parameters*
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
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector of observed object, with respect to origin
            at observer (or the geocenter), corrected for gravitational
            deflection, components in AU.
    
    *References*
        .. [R1] Murray, C.A. (1981) Mon. Notices Royal Ast. Society 195,
            639-648.
        .. [R2] See also formulae in Section B of the Astronomical
            Almanac, or Kaplan, G. et al. (1989) Astronomical Journal
            97, 1197-1210, section iii f.
    
    """

    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obj'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})
    if len(pos_body) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_body'})
    if rmass < 0.0: raise ValueError(_neg_err % {'name': 'rmass'})

    _grav_vec = novaslib.grav_vec
    _grav_vec.argtypes = (POINTER(c_double*3), POINTER(c_double*3),
                          POINTER(c_double*3), c_double, POINTER(c_double*3))
    _grav_vec.restype = None

    pos2 = (c_double*3)()

    _grav_vec(byref((c_double*3)(*pos_obj)), byref((c_double*3)(*pos_obs)),
              byref((c_double*3)(*pos_body)), rmass, byref(pos2))

    return tuple([i for i in pos2])

def aberration(position, vel_earth, lighttime):
    """
    Correct the position vector for aberration of light, including
    relativistic terms.
    
    *Parameters*
        position : tuple of floats, of length 3
            Position vector, referred to origin at center of mass of the
            Earth, components in AU.
        vel_earth : tuple of floats, of length 3
            Velocity vector of center of mass of the Earth, referred to
            origin at solar system barycenter, components in AU/day.
        lighttime : float
            Light time from object to Earth in days.
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, referred to origin at center of mass of the
            Earth, corrected for aberration, components in AU.
    
    *Notes*
        .. [N1] If 'lighttime' = 0 on input, this function will compute
            it.
    
    *References*
        .. [R1] Murray, C. A. (1981) Mon. Notices Royal Ast. Society
            195, 639-648.
    
    """

    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if len(vel_earth) != 3:
        raise ValueError(_vector_len_err % {'name': 'vel_earth'})
    if lighttime < 0.0: raise ValueError(_neg_err % {'name': 'lighttime'})

    _aberration = novaslib.aberration
    _aberration.argtypes = (POINTER(c_double*3), POINTER(c_double*3), c_double,
                            POINTER(c_double*3))
    _aberration.restype = None

    pos2 = (c_double*3)()

    _aberration(byref((c_double*3)(*position)),
                byref((c_double*3)(*vel_earth)), lighttime, byref(pos2))

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
    
    *Parameters*
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
    
    *Returns*
        rv : float
            The observed radial velocity measure times the speed of
            light, in kilometers/second.
    
    *Notes*
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
    
    *References*
        .. [R1] Lindegren & Dravins (2003), Astronomy & Astrophysics
            401, 1185-1201.
    
    """

    if len(pos) != 3: raise ValueError(_vector_len_err % {'name': 'pos'})
    if len(vel) != 3: raise ValueError(_vector_len_err % {'name': 'vel'})
    if len(vel_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'vel_obs'})
    if distance_geo < 0.0:
        raise ValueError(_neg_err % {'name': 'distance_geo'})
    if distance_sun < 0.0:
        raise ValueError(_neg_err % {'name': 'distance_sun'})
    if distance_obj_sun < 0.0:
        raise ValueError(_neg_err % {'name': 'distance_obj_sun'})

    _rad_vel = novaslib.rad_vel
    _rad_vel.argtypes = (POINTER(Object), POINTER(c_double*3),
                         POINTER(c_double*3), POINTER(c_double*3), c_double,
                         c_double, c_double, POINTER(c_double))
    _rad_vel.restype = None

    rv = c_double()

    _rad_vel(byref(cel_object), byref((c_double*3)(*pos)),
             byref((c_double*3)(*vel)), byref((c_double*3)(*vel_obs)),
             distance_geo, distance_sun, distance_obj_sun, byref(rv))

    return rv.value

def precession(jd_tdb1, position, jd_tdb2):
    """
    Precesses equatorial rectangular coordinates from one epoch to
    another. One of the two epochs must be J2000.0. The coordinates
    are referred to the mean dynamical equator and equinox of the two
    respective epochs.
    
    *Parameters*
        jd_tdb1 : float
            TDB Julian date of first epoch. See Note [N1]_ below.
        position : tuple of floats, of length 3
            Position vector, geocentric equatorial rectangular
            coordinates, referred to mean dynamical equator and equinox
            of first epoch.
        jd_tdb2 : float
            TDB Julian date of second epoch. See Note [N1]_ below.
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, geocentric equatorial rectangular
            coordinates, referred to mean dynamical equator and equinox
            of second epoch.
    
    *Notes*
        .. [N1] Either 'date' or 'newdate' must be 2451545.0 (J2000.0)
            TDB.
    
    *References*
        .. [R1] Explanatory Supplement To The Astronomical Almanac, pp.
            103-104.
        .. [R2] Capitaine, N. et al. (2003), Astronomy And Astrophysics
            412, pp. 567-586.
        .. [R3] Hilton, J. L. et al. (2006), IAU WG report, Celest.
            Mech., 94, pp. 351-367.
    
    """

    if jd_tdb1 < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb1'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if jd_tdb2 < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb2'})

    _precession = novaslib.precession
    _precession.argtypes = (c_double, POINTER(c_double*3), c_double,
                            POINTER(c_double*3))
    _precession.restype = c_short

    pos2 = (c_double*3)()

    return_value = _precession(jd_tdb1, byref((c_double*3)(*position)),
                               jd_tdb2, byref(pos2))

    return tuple([i for i in pos2])

def nutation(jd_tdb, position, direction=0, accuracy=0):
    """
    Nutates equatorial rectangular coordinates from mean equator and
    equinox of epoch to true equator and equinox of epoch. Inverse
    transformation may be applied by setting flag 'direction'.
    
    *Parameters*
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
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, geocentric equatorial rectangular
            coordinates, referred to true equator and equinox of epoch.
    
    *References*
        .. [R1] Explanatory Supplement To The Astronomical Almanac, pp.
            114-115.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if direction not in [0, 1]:
        raise ValueError(_option_err % {'name': 'direction',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _nutation = novaslib.nutation
    _nutation.argtypes = (c_double, c_short, c_short, POINTER(c_double*3),
                          POINTER(c_double*3))
    _nutation.restype = None

    pos2 = (c_double*3)()

    _nutation(jd_tdb, direction, accuracy, byref((c_double*3)(*position)),
              byref(pos2))

    return tuple([i for i in pos2])

def nutation_angles(time, accuracy=0):
    """
    This function returns the values for nutation in longitude and
    nutation in obliquity for a given TDB Julian date. The nutation
    model selected depends upon the input value of 'accuracy'. See
    notes below for important details.
    
    *Parameters*
        time : float
            TDB time in Julian centuries since J2000.0
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            nutation.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        (dpsi, deps) : tuple of floats
            Nutation in (longitude, obliquity) in arcseconds.
    
    *Notes*
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
    
    *References*
        .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _nutation_angles = novaslib.nutation_angles
    _nutation_angles.argtypes = (c_double, c_short, POINTER(c_double),
                                 POINTER(c_double))
    _nutation_angles.restype = None

    dpsi = c_double()
    deps = c_double()

    _nutation_angles(time, accuracy, byref(dpsi), byref(deps))

    return (dpsi.value, deps.value)

def fund_args(time):
    """
    To compute the fundamental arguments (mean elements) of the Sun and
    Moon.
    
    *Parameters*
        time : float
            TDB time in Julian centuries since J2000.0
    
    *Returns*
        (l, l', F, D, omega) : tuple
            l = mean anomaly of the Moon
            l' = mean anomaly of the Sun
            F = mean argument of the latitude of the Moon
            D = mean elongation of the Moon from the Sun
            omega = mean longitude of the Moon's ascending node;
                from Simon section 3.4(b.3),
                precession = 5028.8200 arcsec/cy
    
    *References*
        .. [R1] Simon et al. (1994) Astronomy and Astrophysics 282, 663-683,
            esp. Sections 3.4-3.5.
    
    """

    _fund_args = novaslib.fund_args
    _fund_args.argtypes = (c_double, POINTER(c_double*5))
    _fund_args.restype = None

    a = (c_double*5)()

    _fund_args(time, byref(a))

    return tuple([i for i in a])

def mean_obliq(jd_tdb):
    """
    Return the mean obliquity of the ecliptic.
    
    *Parameters*
        jd_tdb : float
            TDB Julian date.
    
    *Returns*
        epsilon : float
            Mean obliquity of the ecliptic in arcseconds.
    
    *References*
        .. [R1] Capitaine et al. (2003), Astronomy and Astrophysics 412,
            567-586.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})

    _mean_obliq = novaslib.mean_obliq
    _mean_obliq.argtypes = (c_double,)
    _mean_obliq.restype = c_double

    obliq = _mean_obliq(jd_tdb)

    return obliq

def vector2radec(position):
    """
    Converts an vector in equatorial rectangular coordinates to
    equatorial spherical coordinates.
    
    *Parameters*
        position : tuple of floats, of length 3
            Position vector, equatorial rectangular coordinates.
    
    *Returns*
        (rightascension, declination) : tuple of floats
            (Right ascension in hours, declination in degrees)
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS Version
            C3.1', C38-C39.
    
    """

    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})

    _vector2radec = novaslib.vector2radec
    _vector2radec.argtypes = (POINTER(c_double*3), POINTER(c_double),
                              POINTER(c_double))
    _vector2radec.restype = c_short

    ra = c_double()
    dec = c_double()

    return_value = _vector2radec(byref((c_double*3)(*position)),
                                 byref(ra), byref(dec))

    return (ra.value, dec.value)

def radec2vector(ra, dec, distance):
    """
    Converts equatorial spherical coordinates to a vector (equatorial
    rectangular coordinates).
    
    *Parameters*
        ra : float
            Right ascension in hours.
        dec : float
            Declination in degrees.
        distance : float
            Distance in AU.
    
    *Returns*
        position : tuple of floats, of length 3
            Position vector, equatorial rectangular coordinates in AU.
    
    """

    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err % {'name': 'rai'})
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'deci'})
    if distance < 0.0: raise ValueError(_neg_err % {'name': 'distance'})

    _radec2vector = novaslib.radec2vector
    _radec2vector.argtypes = (c_double, c_double, c_double,
                              POINTER(c_double*3))
    _radec2vector.restype = None

    vector = (c_double*3)()

    _radec2vector(ra, dec, distance, byref(vector))

    return tuple([i for i in vector])

def starvectors(star):
    """
    Converts angular quantities for stars to vectors.
    
    *Parameters*
        star : CatEntry
            Instance of CatEntry type object containing ICRS catalog
            data.
    
    *Returns*
        (position, velocity) : tuple of tuple of floats, of length 3
            Position and velocity vectors in equatorial rectangular
            coordinates. Components are in AU and AU/day.
    
    *References*
        .. [R1] Bangert, J. et. al. (2011), 'User's Guide to NOVAS
            Version C3.1', C16, C113.
    
    """

    _starvectors = novaslib.starvectors
    _starvectors.argtypes = (POINTER(CatEntry), POINTER(c_double*3),
                             POINTER(c_double*3))
    _starvectors.restype = None

    pos = (c_double*3)()
    vel = (c_double*3)()

    _starvectors(byref(star), byref(pos), byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])

def tdb2tt(tdb_jd):
    """
    Computes the Terrestrial Time (TT) or Terrestrial Dynamical Time
    (TDT) Julian date corresponding to a Barycentric Dynamical Time
    (TDB) Julian date.
    
    *Parameters*
        tdb_jd : float
            TDB Julian date.
    
    *Returns*
        (tt_jd, secdiff) : tuple of float
            tt_jd = TT Julian date, secdiff = difference 'tdb_jd' -
            'tt_jd', in seconds.
    
    *Notes*
        .. [N1] Expression used in this function is a truncated form of
            a longer and more precise series given in the first
            reference [R1]_. The result is good to about 10
            microseconds.
    
    *References*
        .. [R1] Fairhead, L. & Bretagnon, P. (1990) Astron. & Astrophys.
            229, 240.
        .. [R2] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if tdb_jd < 0.0: raise ValueError(_neg_err % {'name': 'tdb_jd'})

    _tdb2tt = novaslib.tdb2tt
    _tdb2tt.argtypes = (c_double, POINTER(c_double), POINTER(c_double))
    _tdb2tt.restype = None

    tt_jd = c_double()
    secdiff = c_double()

    _tdb2tt(tdb_jd, byref(tt_jd), byref(secdiff))

    return tt_jd.value, secdiff.value

def cio_ra(jd_tt, accuracy=0):
    """
    Returns the true right ascension of the celestial intermediate
    origin (CIO) at a given TT Julian date. This is -(equation of the
    origins).
    
    *Parameters*
        jd_tt : float
            TT Julian day.
        accuracy : {0, 1}, optional
            Code specifying the relative accuracy of the output
            position.
                = 0 ... full accuracy (default)
                = 1 ... reduced accuracy
    
    *Returns*
        ra_cio : float
            Right ascension of the CIO, with respect to the true equinox
            of date, in hours (+ or -).
    
    *References*
        .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if jd_tt < 0.0: raise ValueError(_neg_err % {'name': 'jd_tt'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _cio_ra = novaslib.cio_ra
    _cio_ra.argtypes = (c_double, c_short, POINTER(c_double))
    _cio_ra.restype = c_short

    ra_cio = c_double()

    return_value = _cio_ra(jd_tt, accuracy, byref(ra_cio))

    return ra_cio.value

def cio_location(jd_tdb, accuracy=0):
    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _cio_location = novaslib.cio_location
    _cio_location.argtypes = (c_double, c_short, POINTER(c_double),
                              POINTER(c_short))
    _cio_location.restype = c_short

    ra_cio = c_double()
    ref_sys - c_short()

    return_value = _cio_location(jd_tdb, accuracy, byref(ra_cio),
                                 byref(ref_sys))

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
    
    *Parameters*
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
    
    *Returns*
        (x, y, z) : tuple of floats
            Unit vector toward (the CIO, the y-direction, north
            celestial pole (CIP)), equatorial rectangular coordinates,
            referred to the GCRS
    
    *Notes*
        .. [N1] This function effectively constructs the matrix C in eq.
            (3) of the reference.
    
    *References*
        .. [R1] Kaplan, G. (2005), US Naval Observatory Circular 179.
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if not 0.0 <= ra_cio < 24.0:
        raise ValueError(_hour_range_err % {'name': 'ra_cio'})
    if system not in [1, 2]:
        raise ValueError(_option_err % {'name': 'system',
                                        'allowed': [1, 2]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _cio_basis = novaslib.cio_basis
    _cio_basis.argtypes = (c_double, c_double, c_short, c_short,
                           POINTER(c_double), POINTER(c_double),
                           POINTER(c_double))
    _cio_basis.restype = c_short

    x = c_double()
    y = c_double()
    z = c_double()

    return_value = _cio_basis(jd_tdb, ra_cio, system, accuracy, byref(x),
                              byref(y), byref(z))

    return x.value, y.value, z.value

def ira_equinox(jd_tdb, equinox, accuracy=0):
    """
    To compute the intermediate right ascension of the equinox at
    the input Julian date, using an analytical expression for the
    accumulated precession in right ascension. For the true equinox,
    the result is the equation of the origins.
    
    *Parameters*
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
    
    *Returns*
        ira_eq : float
            Intermediate right ascension of the equinox, in hours (+ or
            -). If 'equinox' = 1 (i.e. true equinox), then the returned
            value is the equation of the origins.
    
    *References*
        .. [R1] Capitaine, N. et al. (2003), Astronomy and Astrophysics
            412, 567-586, eq. (42).
    
    """

    if jd_tdb < 0.0: raise ValueError(_neg_err % {'name': 'jd_tdb'})
    if equinox not in [0, 1]:
        raise ValueError(_option_err % {'name': 'equinox',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _ira_equinox = novaslib.ira_equinox
    _ira_equinox.argtypes = (c_double, c_short, c_short)
    _ira_equinox.restype = c_double

    ira_eq = _ira_equinox(jd_tdb, equinox, accuracy)

    return ira_eq

def ephemeris(jd, ss_body, origin=1, accuracy=0):
    """
    
    Retrieves the position and velocity of a solar system body from
    a fundamental ephemeris.
    
    *Parameters*
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
    
    *Returns*
        (pos, vel) : tuple of tuples of three floats
            pos is position vector of the body at 'jd_tdb'; equatorial
            rectangular coordinates in AU referred to the ICRS.
            vel is the velocity vector of the body at 'jd_tdb';
            equatorial rectangular system referred the ICRS, in AU/Day.
    
    """

    if jd[0] + jd[1] < 0.0:
        raise ValueError(_neg_err % {'name': 'jd[0]+jd[1]'})
    if origin not in [0, 1]:
        raise ValueError(_option_err % {'name': 'origin',
                                        'allowed': [0, 1]})
    if accuracy not in [0, 1]:
        raise ValueError(_option_err % {'name': 'accuracy',
                                        'allowed': [0, 1]})

    _ephemeris = novaslib.ephemeris
    _ephemeris.argtypes = (c_double*2, POINTER(Object), c_short, c_short,
                           POINTER(c_double*3), POINTER(c_double*3))
    _ephemeris.restype = c_short

    pos = (c_double*3)()
    vel = (c_double*3)()

    return_value = _ephemeris((c_double*2)(*jd), byref(ss_body), origin,
                              accuracy, byref(pos), byref(vel))

    return tuple([i for i in pos]), tuple([j for j in vel])

def transform_hip(hipparcos):
    """
    To convert Hipparcos catalog data at epoch J1991.25 to epoch
    J2000.0, for use within NOVAS. To be used only for Hipparcos or
    Tycho stars with linear space motion. Both input and output data
    is in the ICRS.
    
    *Parameters*
        hipparcos : CatEntry
            Instance of CatEntry type object containing an entry from
            the Hipparcos catalog, at epoch J1991.25, with all members
            having Hipparcos catalog units. See Note [N1]_ below.
    
    *Returns*
        hip_2000 : CatEntry
            Instance of CatEntry type object containing the transformed
            input entry, at epoch J2000.0. See Note [N2]_ below.
    
    *Notes*
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
    _transform_hip.argtypes = (POINTER(CatEntry), POINTER(CatEntry))
    _transform_hip.restype = None

    hip_2000 = CatEntry()

    _transform_hip(byref(hipparcos), byref(hip_2000))

    return hip_2000

def transform_cat(option, date_incat, incat, date_newcat, newcat_id):
    """
    To transform a star's catalog quantities for a change of epoch
    and/or equator and equinox. Also used to rotate catalog
    quantities on the dynamical equator and equinox of J2000.0 to the
    ICRS or vice versa.
    
    *Parameters*
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
    
    *Returns*
        newcat : CatEntry
            Instance of CatEntry type object containing the transformed
            catalog entry, with units as given in the type definition.
    
    *Notes*
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
        raise ValueError(_option_err % {'name': 'option',
                                        'allowed': [1, 2, 3, 4, 5]})
    if date_incat < 0.0: raise ValueError(_neg_err % {'name': 'date_incat'})
    if date_newcat < 0.0:
        raise ValueError(_neg_err % {'name': 'date_newcat'})

    _transform_cat = novaslib.transform_cat
    _transform_cat.argtypes = (c_short, c_double, POINTER(CatEntry), c_double,
                               c_char*4, POINTER(CatEntry))
    _transform_cat.restype = c_short

    newcat = CatEntry()

    return_value = _transform_cat(option, date_incat, byref(incat),
                                  date_newcat,
                                  create_string_buffer(newcat_id, 51),
                                  byref(newcat))

    return newcat

def limb_angle(pos_obj, pos_obs):
    """
    Compute the angle of an object above or below the Earth's limb
    (horizon). The geometric limb is computed, assuming the Earth to be
    an airless sphere (no refraction or oblateness is included). The
    observer can be on or above the Earth. For an observer on the
    surface of the Earth, this function returns the approximate
    unrefracted altitude.
    
    *Parameters*
        pos_obj : tuple of floats, of length 3
            Position vector of observed object, with respect to origin
            at geocenter, components in AU.
        pos_obs : tuple of floats, of length 3
            Position vector of observer, with respect to origin at
            geocenter, components in AU.
    
    *Returns*
        (limb_angle, nadir_angle) : tuple of floats
            Angle of observed object above (+) or below (-) limb in
            degrees and nadir angle of observed object as a fraction of
            apparent radius of limb where: < 1.0, below the limb; = 1.0,
            on the limb; > 1.0, above the limb
    
    """

    if len(pos_obj) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obj'})
    if len(pos_obs) != 3:
        raise ValueError(_vector_len_err % {'name': 'pos_obs'})

    _limb_angle = novaslib.limb_angle
    _limb_angle.argtypes = (c_double*3, c_double*3, POINTER(c_double),
                            POINTER(c_double))
    _limb_angle.restype = None

    limb_ang = c_double()
    nadir_ang = c_double()

    _limb_angle = ((c_double*3)(*pos_obj), (c_double*3)(*pos_obs),
                   byref(limb_ang), byref(nadir_ang))

    return limb_ang.value, nadir_ang.value

def refract(location, zd_obs, atmosphere=1):
    """
    
    This function computes atmospheric refraction in zenith
    distance. This version computes approximate refraction for
    optical wavelengths.
    
    *Parameters*
        location : OnSurface
            Instance of OnSurface type object  containing observer's
            location. This structure also contains weather data
            (optional) for the observer's location.
        zenith : float
            Observed zenith distance, in degrees.
        atmosphere : {1, 2}, optional
                = 1 ... use 'standard' atmospheric conditions (default)
                = 2 ... use atmospheric parameters input in the
                    'location' object.
    
    *Returns*
        refraction : float
            Atmospheric refraction, in degrees.
    
    *Notes*
        .. [N1] This function can be used for planning observations or
            telescope pointing, but should not be used for the reduction
            of precise observations.
    
    *References*
        .. [R1] Explanatory Supplement to the Astronomical Almanac,
            p. 144.
        .. [R2] Bennett, G. (1982), Journal of Navigation (Royal
            Institute) 35, pp. 255-259.
    
    """

    if atmosphere not in [1, 2]:
        raise ValueError(_option_err % {'name': 'atmosphere',
                                        'allowed': [1, 2]})

    _refract = novaslib.refract
    _refract.argtypes = (POINTER(OnSurface), c_short, c_double)
    _refract.restype = c_double

    refraction = _refract(byref(location), atmosphere, zd_obs)

    return refraction

def julian_date(year, month, day, hour=0.0):
    """
    This function will compute the Julian date for a given calendar
    date (year, month, day, hour).
    
    *Parameters*
        year : integer
            Year.
        month : integer
            Month number.
        day : integer
            Day-of-month.
        hour : float, optional
            Hour-of-day.
    
    *Returns*
        jd : float
            Julian day.
    
    *Notes*
        .. [N1] This function makes no checks for a valid input calendar
            date.
        .. [N2] Input calendar date must be Gregorian.
        .. [N3] Input time value can be based on any UT-like time scale
            (UTC, UT1, TT, etc.) -- output Julian date will have the
            same basis.
    
    *References*
        .. [R1] Fliegel, H. & Van Flandern, T. Comm. of the ACM,
            Vol. 11, No. 10, October 1968, p. 657.
    
    """

    if year < -4712: raise ValueError(_jd_year_err % {'name': 'year'})
    if not 1 <= month <= 12:
        raise ValueError(_month_range_err % {'name': 'month'})
    if not 1 <= day <= _days_in_month[month]:
        raise ValueError(_day_range_err %
                         {'name': 'day', 'daysinmonth': _days_in_month[month]})
    if not 0.0 <= hour < 24.0:
        raise ValueError(_hour_range_err % {'name': 'hour'})

    _julian_date = novaslib.julian_date
    _julian_date.argtypes = (c_short, c_short, c_short, c_double)
    _julian_date.restype = c_double

    jd = _julian_date(year, month, day, hour)

    return jd

def cal_date(day):
    """
    Return the Gregorian date for a given Julian day.
    
    *Parameters*
        day : float
            Julian day.
    
    *Returns*
        date : tuple of three ints and a float
            The elements are year, month, day and hour.
    
    *Notes*
        .. [N1] This routine valid for any 'jd' greater than zero.
        .. [N2] Input Julian date can be based on any UT-like time scale
            (UTC, UT1, TT, etc.) -- output time value will have same
            basis.
    
    *References*
        .. [R1] Fliegel, H. & Van Flandern, T. Comm. of the ACM,
            Vol. 11, No. 10, October 1968, p. 657.
    
    """

    if day < 0.0: raise ValueError(_neg_err % {'name': 'day'})

    _cal_date = novaslib.cal_date
    _cal_date.argtypes = (c_double, POINTER(c_short), POINTER(c_short),
                          POINTER(c_short), POINTER(c_double))
    _cal_date.restype = None

    year = c_short()
    month = c_short()
    day = c_short()
    hour = c_double()

    _cal_date(day, byref(year), byref(month), byref(day), byref(hour))

    return year.value, month.value, day.value, hour.value

def norm_ang(angle):
    """
    Normalize angle into the range 0 <= angle < (2 * pi).
    
    *Parameters*
        angle : float
            Input angle in radians.
    
    *Returns*
        norm_angle : float
            The input angle, normalized as described above, in radians.
    
    """

    _norm_ang = novaslib.norm_ang
    _norm_ang.argtypes = (c_double,)
    _norm_ang.restype = c_double

    norm_angle = _norm_ang(angle)

    return norm_angle

def make_cat_entry(star_name, catalog, star_num, ra, dec, pm_ra, pm_dec,
                   parallax, rad_vel):
    """
    Create an instance of CatEntry containing catalog data for a star or
    "star-like" object.
    
    *Parameters*
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
    
    *Returns*
        star : CatEntry
            Instance of CatEntry type object containing the input data.
    
    *Notes*
        .. [N1] This function is equivalent to calling the object with
            arguments,e.g., CatEntry(star_name, catalog, star_num, ...);
            this function exists for the purpose of direct compatibility
            with NOVAS C.
    
    """

    if not 0.0 <= ra < 24.0:
        raise ValueError(_hour_range_err % {'name': 'ra'})
    if not -90.0 <= dec <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'dec'})
    if parallax < 0.0: raise ValueError(_neg_err % {'name': 'parallax'})

    _make_cat_entry = novaslib.make_cat_entry
    _make_cat_entry.argtypes = (c_char*51, c_char*4, c_long, c_double,
                                c_double, c_double, c_double, c_double,
                                c_double, POINTER(CatEntry))
    _make_cat_entry.restype = c_short

    star = CatEntry()

    return_value = _make_cat_entry(create_string_buffer(star_name, 51),
                                   create_string_buffer(catalog, 4), star_num,
                                   ra, dec, pm_ra, pm_dec, parallax, rad_vel,
                                   byref(star))

    return star

def make_object(type, number, name, star_data):
    """
    Makes an instance of Object--specifying a celestial object--based on
    the input parameters.
    
    *Parameters*
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
    
    *Returns*
        object : Object
            Instance of Object type object containing the object
            definition.
    
    *Notes*
        .. [N1] If 'type' = 0 or 'type' = 1, ``None`` may be given for
            'star_data'; in this case, a dummy star CatEntry will be
            constructed for 'star_data' automatically.
        .. [N2] This function is equivalent to calling the object with
            arguments, e.g., Object(type, number, name, star_data);
            this function exists for the purpose of direct compatibility
            with NOVAS C.
    
    """

    if type not in [0, 1, 2]:
        raise ValueError(_option_err % {'name': 'type',
                                        'allowed': [0, 1, 2]})

    _make_object = novaslib.make_object
    _make_object.argtypes = (c_short, c_short, c_char*51, POINTER(CatEntry),
                             POINTER(Object))
    _make_object.restype = c_short

    if not star_data:
        star_data = CatEntry()

    cel_obj = Object()

    return_value = _make_object(type, number, create_string_buffer(name, 51),
                                byref(star_data), byref(cel_obj))

    return cel_obj

def make_observer(location, obs_surface, obs_space):
    """
    Makes an instance of Observer, specifying the location of the
    observer.
    
    *Parameters*
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
    
    *Returns*
        observer : Observer
            Instance of Observer type object specifying the location of
            the observer.
    
    *Notes*
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
        raise ValueError(_option_err % {'name': 'location',
                                        'allowed': [0, 1, 2]})

    _make_observer = novaslib.make_observer
    _make_observer.argtypes = (c_short, POINTER(OnSurface), POINTER(InSpace),
                               POINTER(Observer))
    _make_observer.restype = c_short

    obs = Observer()

    return_value = _make_observer(where, byref(obs_surface), byref(obs_space),
                                  byref(obs))

    return obs

def make_observer_at_geocenter():
    """
    Makes an instance of Observer, specifying an observer at the
    geocenter.
    
    *Parameters*
        None
    
    *Returns*
        observer : Observer
            Instance of Observer type object specifying the location of
            the observer at the geocenter.
    
    *Notes*
        .. [N1] This function is equivalent to calling the object with
            arguments, e.g., Observer(0, None, None); this function
            exists for the purpose of direct compatibility with NOVAS C.
    
    """

    _make_observer_at_geocenter = novaslib.make_observer_at_geocenter
    _make_observer_at_geocenter.argtypes = (POINTER(Observer),)
    _make_observer_at_geocenter.restype = None

    obs_at_geocenter = Observer()

    _make_observer_at_geocenter(byref(obs_at_geocenter))

    return obs_at_geocenter

def make_observer_on_surface(latitude, longitude, height, temperature,
                             pressure):
    """
    Makes an instance of Observer, specifying the location of and
    weather for an observer on the surface of the Earth.
    
    *Parameters*
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
    
    *Returns*
        observer : Observer
            Instance of Observer type object containing the location of
            and weather for an observer on the surface of the Earth.
    
    *Notes*
        .. [N1] This function is equivalent to calling the object with
            arguments, e.g., Observer(1, obs_surface, None);
            this function exists for the purpose of direct compatibility
            with NOVAS C.
    
    """

    if not -90.0 <= latitude <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'latitude'})
    if not -180.0 <= longitude <= 180.0:
        raise ValueError(_az180_range_err % {'name': 'longitude'})

    _make_observer_on_surface = novaslib.make_observer_on_surface
    _make_observer_on_surface.argtypes = (c_double, c_double, c_double,
                                          c_double, c_double,
                                          POINTER(Observer))
    _make_observer_on_surface.restype = None

    obs_on_surface = Observer()

    _make_observer_on_surface(latitude, longitude, height, temperature,
                              pressure, byref(obs_on_surface))

    return obs_on_surface

def make_observer_in_space(position, velocity):
    """
    Makes an instance of Observer, specifying the position and velocity
    of observer situated on a near-Earth spacecraft.
    
    *Parameters*
        position : tuple of floats, of length 3
        Geocentric position vector (x, y, z) in km.
    velocity : tuple of floats, of length 3
        Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.
    
    *Returns*
        observer : Observer
        Instance of Observer type object containing the position and
        velocity of an observer situated on a near-Earth spacecraft.
    
    *Notes*
        .. [N1] Both input vector tuples are with respect to true
            equator and equinox of date.
        .. [N2] This function is equivalent to calling the object with
            arguments, e.g., Observer(2, None, obs_space);
            this function exists for the purpose of direct compatibility
            with NOVAS C.
    
    """

    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if len(velocity) != 3:
        raise ValueError(_vector_len_err % {'name': 'velocity'})

    _make_observer_in_space = novaslib.make_observer_in_space
    _make_observer_in_space.argtypes = (c_double*3, c_double*3,
                                        POINTER(Observer))
    _make_observer_in_space.restype = None

    obs_in_space = Observer()

    _make_observer_in_space((c_double*3)(*position), (c_double*3)(*velocity),
                            byref(obs_in_space))

    return obs_in_space

def make_on_surface(latitude, longitude, height, temperature, pressure):
    """
    Makes an instance of OnSurface, specifying the location of and
    weather for an observer on the surface of the Earth.
    
    *Parameters*
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
    
    *Returns*
        location : OnSurface
            Instance of OnSurface type object containing the location of
            and weather for an observer on the surface of the Earth.
    
    *Notes*
        .. [N1] This function is equivalent to calling the object with
            arguments, e.g., OnSurface(latitude, longitude, height,
            temperature, pressure); this function exists for the purpose
            of direct compatibility with NOVAS C.
    
    """

    if not -90.0 <= latitude <= 90.0:
        raise ValueError(_elev_range_err % {'name': 'latitude'})
    if not -180.0 <= longitude <= 180.0:
        raise ValueError(_az180_range_err % {'name': 'longitude'})

    _make_on_surface = novaslib.make_on_surface
    _make_on_surface.argtypes = (c_double, c_double, c_double, c_double,
                                 c_double, POINTER(OnSurface))
    _make_on_surface.restype = None

    obs_surface = OnSurface()

    _make_on_surface(latitude, longitude, height, temperature, pressure,
                     byref(obs_surface))

    return obs_surface

def make_in_space(position, velocity):
    """
    Makes an instance of InSpace, specifying the position and velocity
    of an observer situated on a near-Earth spacecraft.
    
    *Parameters*
        position : tuple of floats, of length 3
            Geocentric position vector (x, y, z) in km.
        velocity : tuple of floats, of length 3
            Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.
    
    *Returns*
        object : InSpace
            Instance of InSpace type object containing the position and
            velocity of an observer situated on a near-Earth spacecraft.
    
    *Notes*
        .. [N1] Both input vector tuples are with respect to true
            equator and equinox of date.
        .. [N2] This function is equivalent to calling the object with
            arguments, e.g., InSpace(sc_pos, sc_vel);
            this function exists for the purpose of direct compatibility
            with NOVAS C.
    
    """

    if len(position) != 3:
        raise ValueError(_vector_len_err % {'name': 'position'})
    if len(velocity) != 3:
        raise ValueError(_vector_len_err % {'name': 'velocity'})

    _make_in_space = novaslib.make_in_space
    _make_in_space.argtypes = (c_double*3, c_double*3, POINTER(InSpace))
    _make_in_space.restype = None

    obs_space = InSpace()

    _make_in_space((c_double*3)(*position), (c_double*3)(*velocity),
                   byref(obs_space))

    return obs_space