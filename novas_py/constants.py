# -*- coding: utf-8 -*-

# TDB Julian date of epoch J2000.0.
T0 = 2451545.00000000

# Speed of light in meters/second is a defining physical constant.
C = 299792458.0

# Light-time for one astronomical unit (AU) in seconds, from DE-405.
AU_SEC = 499.0047838061

# Speed of light in AU/day. Value is 86400 / AU_SEC.
C_AUDAY = 173.1446326846693

# Astronomical unit in meters. Value is AU_SEC * C.
AU = 1.4959787069098932e+11

# Astronomical Unit in kilometers.
AU_KM = 1.4959787069098932e+8

# Heliocentric gravitational constant in meters^3 / second^2, from DE-405.
GS = 1.32712440017987e+20

# Geocentric gravitational constant in meters^3 / second^2, from DE-405.
GE = 3.98600433e+14

# Radius of Earth in meters from IERS Conventions (2003).
ERAD = 6378136.6

# Earth ellipsoid flattening from IERS Conventions (2003).
# Value is 1 / 298.25642.
F = 0.003352819697896

# Rotational angular velocity of Earth in radians/sec from IERS
# Conventions (2003).
ANGVEL = 7.2921150e-5

# Reciprocal masses of solar system bodies, from DE-405
# (Sun mass / body mass).
RMASS = {
    'Earth-Moon Barycenter': 328900.561400,
    'Mercury': 6023600.0,
    'Venus': 408523.71,
    'Earth': 332946.050895,
    'Mars': 3098708.0,
    'Jupiter': 1047.3486,
    'Saturn': 3497.898,
    'Uranus': 22902.98,
    'Neptune': 19412.24,
    'Pluto': 135200000.0,
    'Sun': 1.0,
    'Moon': 27068700.387534
}

# Value of 2 * pi in radians.
TWOPI = 6.283185307179586476925287

# Number of arcseconds in 360 degrees.
ASEC360 = 1296000.0

# Angle conversion constants.
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321
