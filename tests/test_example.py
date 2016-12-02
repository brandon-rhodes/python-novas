# -*- coding: utf-8 -*-

import unittest

from math import sqrt, sin, cos
from novas.compat import *
from novas.compat.eph_manager import ephem_open
from novas.constants import T0, DEG2RAD

try:
    import novas_de405
except ImportError:
    raise Exception('please run `pip install novas_de405` before running'
                    ' these tests')

jd_begin, jd_end, de_number = ephem_open()

year = 2008
month = 4
day = 24
leap_secs = 33.0
accuracy = 0
error = 0

hour = 10.605
ut1_utc = -0.387845

latitude = 42.0
longitude = -70.0
height = 0.0
temperature = 10.0
pressure = 1010.0

x_pole = -0.002
y_pole = +0.529

#jd_utc = 2454580.941875
#jd_tt = jd_utc + (leap_secs + 32.184) / 86400.0
#jd_ut1 = jd_utc + ut1_utc / 86400.0
#delta_t = 32.184 + leap_secs - ut1_utc

#jd_tdb = jd_tt # approximation good to 0.0017 seconds


class TestCalendarFunctions(unittest.TestCase):
    def setUp(self):
        self.jd_utc = julian_date(year, month, day, hour)
        self.year_utc, self.month_utc, self.day_utc, self.hour_utc = cal_date(2454580.941875)
        self.jd_tt = self.jd_utc + (leap_secs + 32.184) / 86400.0
        self.jd_ut1 = self.jd_utc + ut1_utc / 86400.0
        self.delta_t = 32.184 + leap_secs - ut1_utc
        self.geo_loc = make_on_surface(latitude, longitude, height, temperature, pressure)

    def test_julian_date_result(self):
        self.assertAlmostEqual(self.jd_utc, 2454580.941875, 6)
        self.assertAlmostEqual(self.jd_tt, 2454580.942629, 6)
        self.assertAlmostEqual(self.jd_ut1, 2454580.941871, 6)

    def test_cal_date_result(self):
        self.assertEqual(self.year_utc, 2008)
        self.assertEqual(self.month_utc, 4)
        self.assertEqual(self.day_utc, 24)
        self.assertAlmostEqual(self.hour_utc, 10.605, 3)

    def test_sidereal_time_result(self):
        self.gast = sidereal_time(self.jd_ut1, 0.0, self.delta_t, 1, 1, accuracy)
        self.last = self.gast + self.geo_loc.longitude / 15.0
        if self.last >= 24.0:
            self.last -= 24.0
        if self.last < 0.0:
            self.last += 24.0
        self.assertAlmostEqual(self.gast, 0.79362134148, 10)
        self.assertAlmostEqual(self.last, 20.12695467481, 10)

    def test_era_result(self):
        self.theta = era(self.jd_ut1, 0.0)
        self.assertAlmostEqual(self.theta, 11.7956158462, 10)


class MakeObjectFunctions(unittest.TestCase):
    def test_make_cat_entry(self):
        self.star = make_cat_entry('GMB 1830', 'FK6', 1307, 11.88299133, 37.71867646, 4003.27, -5815.07, 109.21, -98.8)
        self.assertEqual(self.star.starname, 'GMB 1830'.encode())
        self.assertEqual(self.star.catalog, 'FK6'.encode())
        self.assertEqual(self.star.starnumber, 1307)
        self.assertEqual(self.star.ra, 11.88299133)
        self.assertEqual(self.star.dec, 37.71867646)
        self.assertEqual(self.star.promora, 4003.27)
        self.assertEqual(self.star.promodec, -5815.07)
        self.assertEqual(self.star.parallax, 109.21)
        self.assertEqual(self.star.radialvelocity, -98.8)

    def test_make_cat_entry_exceptions(self):
        self.assertRaises(ValueError, make_cat_entry, 'Lorem ipsum dolor sit amet, consectetur adipiscing elit amet', 'xxx', 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.assertRaises(ValueError, make_cat_entry, 'DUMMY', 'Lorem ipsum', 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

class TestStarFunctions(unittest.TestCase):
    def setUp(self):
        self.jd_utc = julian_date(year, month, day, hour)
        self.jd_tt = self.jd_utc + (leap_secs + 32.184) / 86400.0
        self.delta_t = 32.184 + leap_secs - ut1_utc
        self.star = make_cat_entry('GMB 1830', 'FK6', 1307, 11.88299133, 37.71867646, 4003.27, -5815.07, 109.21, -98.8)
        self.geo_loc = make_on_surface(latitude, longitude, height, temperature, pressure)

    def test_app_star_result(self):
        self.ra, self.dec = app_star(self.jd_tt, self.star, accuracy)
        self.assertAlmostEqual(self.ra, 11.8915509892, 10)
        self.assertAlmostEqual(self.dec, 37.6586357955, 10)

    def test_topo_star_result(self):
        self.rat, self.dect = topo_star(self.jd_tt, self.delta_t, self.star, self.geo_loc, accuracy)
        self.assertAlmostEqual(self.rat, 11.8915479153, 10)
        self.assertAlmostEqual(self.dect, 37.6586695456, 10)


class TestPlanetFunctions(unittest.TestCase):
    def setUp(self):
        self.jd_utc = julian_date(year, month, day, hour)
        self.jd_tt = self.jd_utc + (leap_secs + 32.184) / 86400.0
        self.jd_tdb = self.jd_tt
        self.delta_t = 32.184 + leap_secs - ut1_utc
        self.dummy_star = make_cat_entry('DUMMY', 'xxx', 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.moon = make_object(0, 11, 'Moon', self.dummy_star)
        self.mars = make_object(0, 4, 'Mars', self.dummy_star)
        self.geo_loc = make_on_surface(latitude, longitude, height, temperature, pressure)
        self.obs_loc = make_observer_on_surface(latitude, longitude, height, temperature, pressure)

    def test_app_planet_result(self):
        self.ra, self.dec, self.dis = app_planet(self.jd_tt, self.moon, accuracy)
        self.assertAlmostEqual(self.ra, 17.1390774264, 10)
        self.assertAlmostEqual(self.dec, -27.5374448869, 10)
        self.assertAlmostEqual(self.dis, 2.710296515e-03, 10)

    def test_topo_planet_result(self):
        self.rat, self.dect, self.dist = topo_planet(self.jd_tt, self.delta_t, self.moon, self.geo_loc, accuracy)
        self.assertAlmostEqual(self.rat, 17.1031967646, 10)
        self.assertAlmostEqual(self.dect, -28.2902502967, 10)
        self.assertAlmostEqual(self.dist, 2.703785126e-03, 10)

    def test_place_result(self):
        self.t_place = place(self.jd_tt, self.delta_t, self.moon, self.obs_loc, 1, accuracy)
        self.assertAlmostEqual(self.t_place.ra, 17.1031967646, 10)
        self.assertAlmostEqual(self.t_place.dec, -28.2902502967, 10)
        self.assertAlmostEqual(self.t_place.dis, 2.703785126e-03, 10)

    def test_ephemeris_and_transform_result(self):
        self.pos, self.vel = ephemeris((self.jd_tdb, 0.0), self.mars, 1, accuracy)
        self.pose = equ2ecl_vec(T0, self.pos, 2, accuracy,)
        self.elon, self.elat = vector2radec(self.pose)
        self.elon *= 15.0
        self.r = sqrt(self.pose[0]*self.pose[0]+self.pose[1]*self.pose[1]+self.pose[2]*self.pose[2])
        self.assertAlmostEqual(self.elon, 148.0032235906, 10)
        self.assertAlmostEqual(self.elat, 1.8288284075, 10)
        self.assertAlmostEqual(self.r, 1.664218258879, 10)


class TestCoordinateTransformFunctions(unittest.TestCase):
    def setUp(self):
        self.jd_utc = julian_date(year, month, day, hour)
        self.jd_tt = self.jd_utc + (leap_secs + 32.184) / 86400.0
        self.jd_ut1 = self.jd_utc + ut1_utc / 86400.0
        self.delta_t = 32.184 + leap_secs - ut1_utc
        self.dummy_star = make_cat_entry('DUMMY', 'xxx', 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.moon = make_object(0, 11, 'Moon', self.dummy_star)
        self.geo_loc = make_on_surface(latitude, longitude, height, temperature, pressure)

    def test_equ2hor_result(self):
        self.rat, self.dect, self.dist = topo_planet(self.jd_tt, self.delta_t, self.moon, self.geo_loc, accuracy)
        (self.zd, self.az), (self.rar, self.decr) = equ2hor(self.jd_ut1, self.delta_t, 0.0, 0.0, self.geo_loc, self.rat, self.dect, 1, accuracy)
        self.assertAlmostEqual(self.zd, 81.6891016502, 10)
        self.assertAlmostEqual(self.az, 219.2708903405, 10)

    def test_vector_to_gcrs_result(self):
        self.lon_rad = self.geo_loc.longitude * DEG2RAD
        self.lat_rad = self.geo_loc.latitude * DEG2RAD
        self.sin_lon = sin(self.lon_rad)
        self.cos_lon = cos(self.lon_rad)
        self.sin_lat = sin(self.lat_rad)
        self.cos_lat = cos(self.lat_rad)
        self.vter = (self.cos_lat*self.cos_lon, self.cos_lat*self.sin_lon, self.sin_lat)
        self.vcel = ter2cel(self.jd_ut1, 0.0, self.delta_t, x_pole, y_pole, self.vter, 1, 0, accuracy)
        self.ra, self.dec = vector2radec(self.vcel)
        self.assertAlmostEqual(self.ra, 20.1221838608, 10)
        self.assertAlmostEqual(self.dec, 41.9769823554, 10)

if __name__ == '__main__':
    unittest.main()
