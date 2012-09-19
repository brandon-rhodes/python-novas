# -*- coding: utf-8 -*-

from ctypes import *
from novas import novaslib

_neg_err = "%(name)s must be >= 0.0"
_option_err = "%(name)s must be in %(allowed)s"

def solarsystem(tjd, body, origin):
	"""
	Provides an interface between the JPL direct-access solar system
    ephemerides and NOVAS-C.

    *Parameters*
        tjd : float
        	Julian date of the desired time, on the TDB time scale.
        body : integer
        	Body identification number for the solar system object of
        	interest; Mercury = 1, ..., Pluto = 9, Sun = 10, Moon = 11.
        origin : integer
        	Origin code
        		= 0 ... solar system barycenter
        		= 1 ... center of mass of the Sun
        		= 2 ... center of Earth

    *Returns*
        position : tuple of floats, of length 3
        	Position vector of 'body' at 'tjd'; equatorial rectangular
        	coordinates in AU referred to the ICRS.
        velocity : tuple of floats, of length 3
        	Velocity vector of 'body' at 'tjd'; equatorial rectangular
        	system referred to the ICRS, in AU/day.

    *Notes*
        .. [N1] This function and function 'planet_ephemeris' were
            designed to work with the 1997 version of the JPL
            ephemerides, as noted in the references.
        .. [N2] The user must have a JPL binary ephemeris file (see
            package README) and open the file using function
            'ephem_open' prior to calling this function.
        .. [N3] This function places the entire Julian date in the first
            element of the input time to 'planet_ephemeris'. This is
            adequate for all but the highest precision applications. For
            highest precision, use function 'solarsystem_hp'.

    *References*
        .. [R1] JPL. 2007, "JPL Planetary and Lunar Ephemerides: Export
        	Information," (Pasadena, CA: JPL)
        	http://ssd.jpl.nasa.gov/?planet_eph_export.
        .. [R2] Kaplan, G. H. "NOVAS: Naval Observatory Vector
            Astrometry Subroutines"; USNO internal document dated 20 Oct
            1988; revised 15 Mar 1990.

	"""

	if tjd < 0.0: raise ValueError(_neg_err % {'name': 'tjd'})
	if body not in range(1, 12):
		raise ValueError(_option_err % {'name': 'body',
										'allowed': range(1, 12)})
	if origin not in [0, 1, 2]:
		raise ValueError(_option_err % {'name': 'origin',
										'allowed': [0, 1, 2]})

	_solarsystem = novaslib.solarsystem
	_solarsystem.argtypes = (c_double, c_short, c_short, POINTER(c_double*3),
							 POINTER(c_double*3))
	_solarsystem.restype = c_short

	position = (c_double*3)()
	velocity = (c_double*3)()

	return_value = _solarsystem(tjd, body, origin, byref(position),
						        byref(velocity))

	return tuple([i for i in position]), tuple([j for j in velocity])

def solarsystem_hp(tjd, body, origin):
	"""
	Provides an interface between the JPL direct-access solar system
    ephemerides and NOVAS-C for highest precision applications.

    *Parameters*
        tjd : tuple of floats, of length 2
        	Two-element tuple containing the Julian date, which may be
            split any way (although the first element is usually the
            "integer" part, and the second element is the "fractional"
            part). Julian date is on the TDB or "T_eph" time scale.
        body : integer
        	Body identification number for the solar system object of
        	interest; Mercury = 1, ..., Pluto = 9, Sun = 10, Moon = 11.
        origin : integer
        	Origin code
        		= 0 ... solar system barycenter
        		= 1 ... center of mass of the Sun
        		= 2 ... center of Earth

    *Returns*
        position : tuple of floats, of length 3
        	Position vector of 'body' at 'tjd'; equatorial rectangular
        	coordinates in AU referred to the ICRS.
        velocity : tuple of floats, of length 3
        	Velocity vector of 'body' at 'tjd'; equatorial rectangular
        	system referred to the ICRS, in AU/day.

    *Notes*
        .. [N1] This function and function 'planet_ephemeris' were
            designed to work with the 1997 version of the JPL
            ephemerides, as noted in the references.
        .. [N2] The user must have a JPL binary ephemeris file (see
            package README) and open the file using function
            'ephem_open' prior to calling this function.
        .. [N3] This function supports the "split" Julian date feature
            of function 'planet_ephemeris' for highest precision. For
            usual applications, use function 'solarsystem'.

    *References*
        .. [R1] JPL. 2007, "JPL Planetary and Lunar Ephemerides: Export
        	Information," (Pasadena, CA: JPL)
        	http://ssd.jpl.nasa.gov/?planet_eph_export.
        .. [R2] Kaplan, G. H. "NOVAS: Naval Observatory Vector
            Astrometry Subroutines"; USNO internal document dated 20 Oct
            1988; revised 15 Mar 1990.

	"""

	if tjd < 0.0: raise ValueError(_neg_err % {'name': 'tjd'})
	if body not in range(1, 12):
		raise ValueError(_option_err % {'name': 'body',
										'allowed': range(1, 12)})
	if origin not in [0, 1, 2]:
		raise ValueError(_option_err % {'name': 'origin',
										'allowed': [0, 1, 2]})

	_solarsystem_hp = novaslib.solarsystem_hp
	_solarsystem_hp.argtypes = (c_double*2, c_short, c_short,
								POINTER(c_double*3), POINTER(c_double*3))
	_solarsystem_hp.restype = c_short

	position = (c_double*3)()
	velocity = (c_double*3)()

	return_value = _solarsystem_hp((c_double*2)(*tjd), body, origin,
								   byref(position), byref(velocity))

	return tuple([i for i in position]), tuple([j for j in velocity])