# -*- coding: utf-8 -*-

import os
import sys
from ctypes import *
from novas import novaslib

_neg_err = "%(name)s must be >= 0.0"
_option_err = "%(name)s must be in %(allowed)s"

def ephem_open(ephem_name=None):
    """
    This function opens a JPL planetary ephemeris file and sets initial
    values. This function must be called prior to calls to the other JPL
    ephemeris functions.

    *Parameters*
        ephem_name : string (optional)
            Name (including path if necessary) of direct-access JPL
            ephemeris file. This argument is optional if the full path
            to the ephemeris file is set in an environment variable
            called EPHEMERIS_FILE.  If EPHEMERIS_FILE is set but
            ``ephem_name`` is given when the function is called, the
            function will use the value of ``ephem_name`` and ignore
            EPHEMERIS_FILE.

    *Returns*
        jd_begin : float
            Beginning Julian date of the ephemeris file.
        jd_end : float
            Ending Julian date of the ephemeris file.
        de_number : integer
            DE number of the ephemeris file opened.

    *References*
        .. [R1] Standish, E.M. and Newhall, X X (1988). "The JPL Export
            Planetary Ephemeris"; JPL document dated 17 June 1988.

    """

    _ephem_open = novaslib.ephem_open
    _ephem_open.argtypes = (c_char_p, POINTER(c_double), POINTER(c_double),
                            POINTER(c_short))
    _ephem_open.restype = c_short

    jd_begin = c_double()
    jd_end = c_double()
    de_number = c_short()

    # Begin customizations by Brandon Rhodes for release on PyPI
    if ephem_name is None and 'EPHEMERIS_FILE' not in os.environ:
        try:
            import novas_de405
        except ImportError:
            pass
        else:
            directory = os.path.dirname(novas_de405.__file__)
            ephem_name = os.path.join(directory, 'DE405.bin')
    # End customizations by Brandon Rhodes for release on PyPI

    if ephem_name is None:
        try:
            ephem_name = os.environ['EPHEMERIS_FILE']
        except KeyError:
            message = """ephem_open needs the full path to a binary ephemeris
                         file. Either pass the path to ephem_open or set an
                         environment variable called EPHEMERIS_FILE to the
                         path. If you choose the later you may continue to call
                         ephem_open without any arguments."""
            raise IOError(message)

    return_value = _ephem_open(c_char_p(ephem_name), byref(jd_begin),
                               byref(jd_end), byref(de_number))

    if return_value == 1:
        raise IOError("Ephemeris file %s not found" % ephem_name)
    elif return_value in range(2, 11):
        raise IOError("Error reading from ephemeris file header")
    elif return_value == 11:
        raise IOError("Unable to set record length; ephemeris not supported")

    return jd_begin.value, jd_end.value, de_number.value

def planet_ephemeris(tjd, target, center):
    """
    This function accesses the JPL planetary ephemeris to give the
    position and velocity of the target object with respect to the
    center object.

    *Parameters*
        tjd : tuple of floats, of length 2
            Julian date split into two parts. 'tjd' may be split any way
            (although the first element is usually the "integer" part,
            and the second element is the "fractional" part).  Julian
            date is in the TDB or "T_eph" time scale.
        target : {0, ..., 13}
            Number of 'target' point.
        center : {0, ..., 13}
            Number of 'center' (origin) point.

            The numbering convention for 'target' and'center' is:
                = 0 ... Mercury
                = 1 ... Venus
                = 2 ... Earth
                = 3 ... Mars
                = 4 ... Jupiter
                = 5 ... Saturn
                = 6 ... Uranus
                = 7 ... Neptune
                = 8 ... Pluto
                = 9 ... Moon
                = 10 ... Sun
                = 11 ... Solar system bary.
                = 12 ... Earth-Moon bary.
                = 13 ... Nutations (long int. and obliq.)
                (If nutations are desired, set 'target' = 13;
                 'center' will be ignored on that call.)

    *Returns*
        position : tuple of floats, of length 3
            Position vector array of target relative to center, measured
            in AU.
        velocity : tuple of floats, of length 3
            Velocity vector array of target relative to center, measured
            in AU/day.

    *References*
        .. [R1] Standish, E.M. and Newhall, X X (1988). "The JPL Export
            Planetary Ephemeris"; JPL document dated 17 June 1988.

    """

    if tjd[0] + tjd[1] < 0:
        raise ValueError(_neg_err % {'name': 'tjd[0]+tjd[1]'})
    if target not in range(0, 14):
        raise ValueError(_option_err % {'name': 'target',
                                        'allowed': range(0, 14)})
    if center not in range(0, 14):
        raise ValueError(_option_err % {'name': 'center',
                                        'allowed': range(0, 14)})

    _planet_ephemeris = novaslib.planet_ephemeris
    _planet_ephemeris.argtypes = (c_double*2, c_short, c_short,
                                  POINTER(c_double*3), POINTER(c_double*3))
    _planet_ephemeris.restype = c_short

    position = (c_double*3)()
    velocity = (c_double*3)()

    return_value = _planet_ephemeris((c_double*2)(*tjd), target, center,
                                     byref(position), byref(velocity))

    if return_value == 1:
        raise IOError("Error reading from ephemeris file")
    elif return_value == 2:
        raise ValueError("tjd[0]+tjd[1] does not fall within the time span of \
                         the open ephemeris file")

    return tuple([i for i in position]), tuple([j for j in velocity])

def state(jed, target):
    """
    This function reads and interpolates the JPL planetary ephemeris
    file.

    *Parameters*
        jed : tuple of floats, of length 2
            Julian date (TDB) at which interpolation is wanted. Any
            combination of jed[0]+jed[1] which falls within the time
            span on the file is a permissible epoch. See Note [N1]_
            below.
        target : {0, ..., 10}
            The requested body to get data for from the ephemeris file.
            The designation of the astronomical bodies is:
                = 0 ... Mercury
                = 1 ... Venus
                = 2 ... Earth-Moon barycenter
                = 3 ... Mars
                = 4 ... Jupiter
                = 5 ... Saturn
                = 6 ... Uranus
                = 7 ... Neptune
                = 8 ... Pluto
                = 9 ... geocentric Moon
                = 10 ... Sun

    *Returns*
        target_pos : tuple of floats, of length 3
            The barycentric position vector array of the requested
            object, in AU.
        target_vel : tuple of floats, of length 3
            The barycentric velocity vector array of the requested
            object, in AU/Day.

        Both vectors are referenced to the Earth mean equator and
        equinox of epoch.

    *Notes*
        .. [N1] For ease in programming, the user may put the entire
            epoch in jed[0] and set jed[1] = 0. For maximum
            interpolation accuracy set jed[0] = the most recent midnight
            at or before interpolation epoch, and set jed[1] =
            fractional part of a day elapsed between jed[0] and epoch.
            As an alternative, it may prove convenient to set jed[0] =
            some fixed epoch, such as start of the integration and
            jed[1] = elapsed interval between then and epoch.

    *References*
        .. [R1] Standish, E.M. and Newhall, X X (1988). "The JPL Export
            Planetary Ephemeris"; JPL document dated 17 June 1988.

    """

    if tjd[0] + tjd[1] < 0:
        raise ValueError(_neg_err % {'name': 'tjd[0]+tjd[1]'})
    if target not in range(0, 11):
        raise ValueError(_option_err % {'name': 'target',
                                        'allowed': range(0, 11)})

    _state = novaslib.state
    _state.argtypes = (POINTER(c_double*2), c_short, POINTER(c_double*3),
                       POINTER(c_double*3))
    _state.restype = c_short

    target_pos = (c_double*3)()
    target_vel = (c_double*3)()

    return_value = _state(byref((c_double*2)(*jed)), target, byref(target_pos),
                          byref(target_vel))

    if return_value == 1:
        raise IOError("Error reading from ephemeris file")
    elif return_value == 2:
        raise ValueError("jed[0]+jed[1] does not fall within the time span of \
                         the open ephemeris file")

    return tuple([i for i in target_pos]), tuple([j for j in target_vel])

def split(tt):
    """
    This function breaks up a float number into a double integer part
    and a fractional part.

    *Parameters*
        tt : float
            Input number.

    *Returns*
        fr : tuple of floats, of length 2
            Output tuple;
                fr[0] contains integer part,
                fr[1] contains fractional part.
            For negative input numbers,
                fr[0] contains the next more negative integer;
                fr[1] contains a positive fraction.

    *References*
        .. [R1] Standish, E.M. and Newhall, X X (1988). "The JPL Export
            Planetary Ephemeris"; JPL document dated 17 June 1988.

    """

    _split = novaslib.split
    _split.argtypes = (c_double, POINTER(c_double*2))
    _split.restype = None

    fr = (c_double*2)()

    _split(tt, byref(fr))

    return tuple([i for i in fr])