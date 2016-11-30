# -*- coding: utf-8 -*-

import unittest

from novas.compat import *


class TestPlace(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.cel_object = make_object(2, 0, "Polaris", CatEntry("POLARIS".encode(), "HIP".encode(), 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4))
        self.location = make_observer_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        self.delta_t = 67.0
        self.system = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            sky_pos = place(self.jd_tt, self.delta_t, self.cel_object, self.location, self.system)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            sky_pos = place(self.jd_tt, self.delta_t, self.cel_object, self.location, self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            sky_pos = place(self.jd_tt, self.delta_t, self.cel_object, self.location, self.system, self.accuracy)


class TestSiderealTime(unittest.TestCase):
    def setUp(self):
        self.jd_high = 2455519.0
        self.jd_low = 0.5
        self.delta_t = 67.0
        self.gst_type = 0
        self.method = 1
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_high = -2455519.0

        with self.assertRaises(ValueError):
            gst = sidereal_time(self.jd_high, self.jd_low, self.delta_t, self.gst_type, self.method)

    def test_gst_type_value(self):
        self.gst_type = 99

        with self.assertRaises(ValueError):
            gst = sidereal_time(self.jd_high, self.jd_low, self.delta_t, self.gst_type, self.method)

    def test_method_value(self):
        self.method = 99

        with self.assertRaises(ValueError):
            gst = sidereal_time(self.jd_high, self.jd_low, self.delta_t, self.gst_type, self.method)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            gst = sidereal_time(self.jd_high, self.jd_low, self.delta_t, self.gst_type, self.method, self.accuracy)


class TestTer2Cel(unittest.TestCase):
    def setUp(self):
        self.jd_ut_high = 2455519.0
        self.jd_ut_low = 0.5
        self.delta_t = 67.0
        self.xp = 0.0
        self.yp = 0.0
        self.vec1 = (0.0, 0.0, 0.0)
        self.method = 0
        self.option = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_ut_high = -2455519.0

        with self.assertRaises(ValueError):
            vec2 = ter2cel(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_vec1_tuple_length(self):
        self.vec1 = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            vec2 = ter2cel(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_method_value(self):
        self.method = 99

        with self.assertRaises(ValueError):
            vec2 = ter2cel(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_option_value(self):
        self.option = 99

        with self.assertRaises(ValueError):
            vec2 = ter2cel(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            vec2 = ter2cel(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)


class TestCel2Ter(unittest.TestCase):
    def setUp(self):
        self.jd_ut_high = 2455519.0
        self.jd_ut_low = 0.5
        self.delta_t = 67.0
        self.xp = 0.0
        self.yp = 0.0
        self.vec1 = (0.0, 0.0, 0.0)
        self.method = 0
        self.option = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_ut_high = -2455519.0

        with self.assertRaises(ValueError):
            vec2 = cel2ter(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_vec1_tuple_length(self):
        self.vec1 = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            vec2 = cel2ter(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_method_value(self):
        self.method = 99

        with self.assertRaises(ValueError):
            vec2 = cel2ter(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_option_value(self):
        self.option = 99

        with self.assertRaises(ValueError):
            vec2 = cel2ter(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            vec2 = cel2ter(self.jd_ut_high, self.jd_ut_low, self.delta_t, self.xp, self.yp, self.vec1, self.method, self.option, self.accuracy)


class TestEqu2Hor(unittest.TestCase):
    def setUp(self):
        self.jd_ut1 = 2455519.5
        self.delta_t = 67.0
        self.xp = 0.0
        self.yp = 0.0
        self.location = make_on_surface(40.0, -105.0, 1800.0, 20.0, 1020.0)
        self.ra = 0.0
        self.dec = 0.0
        self.ref_option = 0
        self.accuracy= 0

    def test_negative_julian_date(self):
        self.jd_ut1 = -2455519.0

        with self.assertRaises(ValueError):
            (zd, az), (rar, decr) = equ2hor(self.jd_ut1, self.delta_t, self.xp, self.yp, self.location, self.ra, self.dec, self.ref_option, self.accuracy)

    def test_ref_option_value(self):
        self.ref_option = 99

        with self.assertRaises(ValueError):
            (zd, az), (rar, decr) = equ2hor(self.jd_ut1, self.delta_t, self.xp, self.yp, self.location, self.ra, self.dec, self.ref_option, self.accuracy)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            (zd, az), (rar, decr) = equ2hor(self.jd_ut1, self.delta_t, self.xp, self.yp, self.location, self.ra, self.dec, self.ref_option, self.accuracy)


class TestTransformCat(unittest.TestCase):
    def setUp(self):
        pass


class TestTransformHip(unittest.TestCase):
    def setUp(self):
        pass


class TestAppStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.star = make_cat_entry("POLARIS", "HIP", 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = app_star(self.jd_tt, self.star)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = app_star(self.jd_tt, self.star, self.accuracy)


class TestTopoStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.delta_t = 67.0
        self.star = make_cat_entry("POLARIS", "HIP", 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4)
        self.position = make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = topo_star(self.jd_tt, self.delta_t, self.star, self.position)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = topo_star(self.jd_tt, self.delta_t, self.star, self.position, self.accuracy)


class TestVirtualStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.star = make_cat_entry("POLARIS", "HIP", 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = virtual_star(self.jd_tt, self.star)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = virtual_star(self.jd_tt, self.star, self.accuracy)


class TestLocalStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.delta_t = 67.0
        self.star = make_cat_entry("POLARIS", "HIP", 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4)
        self.position = make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = local_star(self.jd_tt, self.delta_t, self.star, self.position)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = local_star(self.jd_tt, self.delta_t, self.star, self.position, self.accuracy)


class TestAstroStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.star = make_cat_entry("POLARIS", "HIP", 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = astro_star(self.jd_tt, self.star)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = astro_star(self.jd_tt, self.star, self.accuracy)


class TestMeanStar(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.ra = 12.0
        self.dec = 0.0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ira, idec = mean_star(self.jd_tt, self.ra, self.dec)

    def test_ra_out_of_range(self):
        self.ra = 25.0

        with self.assertRaises(ValueError):
            ira, idec = mean_star(self.jd_tt, self.ra, self.dec)

    def test_dec_out_of_range(self):
        self.dec = 100.0

        with self.assertRaises(ValueError):
            ira, idec = mean_star(self.jd_tt, self.ra, self.dec)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ira, idec = mean_star(self.jd_tt, self.ra, self.dec, self.accuracy)


class TestAppPlanet(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.ss_body = make_object(0, 11, "Moon", CatEntry())
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec, dis = app_planet(self.jd_tt, self.ss_body)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec, dis = app_planet(self.jd_tt, self.ss_body, self.accuracy)


class TestTopoPlanet(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.delta_t = 67.0
        self.ss_body = make_object(0, 11, "Moon", None)
        self.position = make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec, dis = topo_planet(self.jd_tt, self.delta_t, self.ss_body, self.position)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec, dis = topo_planet(self.jd_tt, self.delta_t, self.ss_body, self.position, self.accuracy)


class TestVirtualPlanet(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.ss_body = make_object(0, 11, "Moon", None)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec, dis = virtual_planet(self.jd_tt, self.ss_body)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec, dis = virtual_planet(self.jd_tt, self.ss_body, self.accuracy)


class TestLocalPlanet(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.delta_t = 67.0
        self.ss_body = make_object(0, 11, "Moon", None)
        self.position = make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec, dis = local_planet(self.jd_tt, self.delta_t, self.ss_body, self.position)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec, dis = local_planet(self.jd_tt, self.delta_t, self.ss_body, self.position, self.accuracy)


class TestAstroPlanet(unittest.TestCase):
    def setUp(self):
        self.jd_tt = 2455519.5
        self.ss_body = make_object(0, 11, "Moon", None)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_tt = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec, dis = astro_planet(self.jd_tt, self.ss_body)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec, dis = astro_planet(self.jd_tt, self.ss_body, self.accuracy)


class TestAberration(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.velocity = (0.0, 0.0, 0.0)
        self.lighttime = 1.0

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = aberration(self.position, self.velocity, self.lighttime)

    def test_velocity_tuple_length(self):
        self.velocity = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = aberration(self.position, self.velocity, self.lighttime)

    def test_negative_lighttime(self):
        self.lighttime = -1.0

        with self.assertRaises(ValueError):
            x, y, z = aberration(self.position, self.velocity, self.lighttime)


class TestBary2Obs(unittest.TestCase):
    def setUp(self):
        self.position1 = (0.0, 0.0, 0.0)
        self.position_obs = (0.0, 0.0, 0.0)

    def test_position1_tuple_length(self):
        self.position1 = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            position, lighttime = bary2obs(self.position1, self.position_obs)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            position, lighttime = bary2obs(self.position1, self.position_obs)


class TestCIOBasis(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.ra_cio = 12.0
        self.system = 1
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = cio_basis(self.date, self.ra_cio, self.system)

    def test_ra_cio_range_high(self):
        self.ra_cio = 25.0

        with self.assertRaises(ValueError):
            x, y, z = cio_basis(self.date, self.ra_cio, self.system)

    def test_ra_cio_range_low(self):
        self.ra_cio = -1.0

        with self.assertRaises(ValueError):
            x, y, z = cio_basis(self.date, self.ra_cio, self.system)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            x, y, z = cio_basis(self.date, self.ra_cio, self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            x, y, z = cio_basis(self.date, self.ra_cio, self.system, self.accuracy)


class TestDLight(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.position_obs = (0.0, 0.0, 0.0)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            diflt = d_light(self.position, self.position_obs)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            diflt = d_light(self.position, self.position_obs)


class TestEcl2EquVec(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.system = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = ecl2equ_vec(self.date, self.position)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = ecl2equ_vec(self.date, self.position)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            x, y, z = ecl2equ_vec(self.date, self.position, coord_sys=self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            x, y, z = ecl2equ_vec(self.date, self.position, accuracy=self.accuracy)


class TestEqu2Ecl(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.ra = 12.0
        self.dec = 0.0
        self.system = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec)

    def test_ra_range_high(self):
        self.ra = 25.0

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec)

    def test_ra_range_low(self):
        self.ra = -1.0

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec)

    def test_dec_range_high(self):
        self.dec = 100.0

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec)

    def test_dec_range_low(self):
        self.dec = -100.0

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec, coord_sys=self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            elon, elat = equ2ecl(self.date, self.ra, self.dec, accuracy=self.accuracy)


class TestEqu2EclVec(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.system = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = equ2ecl_vec(self.date, self.position)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = equ2ecl_vec(self.date, self.position)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            x, y, z = equ2ecl_vec(self.date, self.position, coord_sys=self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            x, y, z = equ2ecl_vec(self.date, self.position, accuracy=self.accuracy)


class TestEqu2Gal(unittest.TestCase):
    def setUp(self):
        self.ra = 12.0
        self.dec = 0.0

    def test_ra_range_high(self):
        self.ra = 25.0

        with self.assertRaises(ValueError):
            glon, glat = equ2gal(self.ra, self.dec)

    def test_ra_range_low(self):
        self.ra = -1.0

        with self.assertRaises(ValueError):
            glon, glat = equ2gal(self.ra, self.dec)

    def test_dec_range_high(self):
        self.dec = 100.0

        with self.assertRaises(ValueError):
            glon, glat = equ2gal(self.ra, self.dec)

    def test_dec_range_low(self):
        self.dec = -100.0

        with self.assertRaises(ValueError):
            glon, glat = equ2gal(self.ra, self.dec)


class TestEra(unittest.TestCase):
    def setUp(self):
        self.date_high = 2455519.0
        self.date_low = 0.0

    def test_negative_julian_date(self):
        self.date_high = -2455519.0

        with self.assertRaises(ValueError):
            theta = era(self.date_high)


class TestFrameTie(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.direction = 0

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = frame_tie(self.position)

    def test_direction_value(self):
        self.direction = 99

        with self.assertRaises(ValueError):
            x, y, z = frame_tie(self.position, self.direction)


class TestGcrs2Equ(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.rag = 12.0
        self.decg = 0.0
        self.system = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg)

    def test_rag_range_high(self):
        self.rag = 25.0

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg)

    def test_rag_range_low(self):
        self.rag = -1.0

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg)

    def test_decg_range_high(self):
        self.decg = 100.0

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg)

    def test_decg_range_low(self):
        self.decg = -100.0

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg)

    def test_system_value(self):
        self.system = 99

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg, coord_sys=self.system)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra, dec = gcrs2equ(self.date, self.rag, self.decg, accuracy=self.accuracy)


class TestGeoPosVel(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.delta_t = 67.0
        self.observer = make_observer_on_surface(40.0, -105.0, 1800.0, 20.0, 1020.0)
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            position, velocity = geo_posvel(self.date, self.delta_t, self.observer)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            position, velocity = geo_posvel(self.date, self.delta_t, self.observer, self.accuracy)


class TestGravDef(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.position_obs = (0.0, 0.0, 0.0)
        self.location = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = grav_def(self.date, self.position, self.position_obs, self.location)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = grav_def(self.date, self.position, self.position_obs, self.location)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = grav_def(self.date, self.position, self.position_obs, self.location)

    def test_location_value(self):
        self.location = 99

        with self.assertRaises(ValueError):
            x, y, z = grav_def(self.date, self.position, self.position_obs, self.location)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            x, y, z = grav_def(self.date, self.position, self.position_obs, self.location, self.accuracy)


class TestGravVec(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.position_obs = (0.0, 0.0, 0.0)
        self.position_body = (0.0, 0.0, 0.0)
        self.rmass = 0.5

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = grav_vec(self.position, self.position_obs, self.position_body, self.rmass)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = grav_vec(self.position, self.position_obs, self.position_body, self.rmass)

    def test_position_body_tuple_length(self):
        self.position_body = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = grav_vec(self.position, self.position_obs, self.position_body, self.rmass)

    def test_negative_rmass(self):
        self.rmass = -0.5

        with self.assertRaises(ValueError):
            x, y, z = grav_vec(self.position, self.position_obs, self.position_body, self.rmass)


class TestLightTime(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.object = make_object(0, 11, 'Moon', None)
        self.position_obs = (0.0, 0.0, 0.0)
        self.estimate = 0.1
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            position, time = light_time(self.date, self.object, self.position_obs, self.estimate)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            position, time = light_time(self.date, self.object, self.position_obs, self.estimate)

    def test_negative_estimate(self):
        self.estimate = -1.0

        with self.assertRaises(ValueError):
            position, time = light_time(self.date, self.object, self.position_obs, self.estimate)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            position, time = light_time(self.date, self.object, self.position_obs, self.estimate, self.accuracy)


class TestLimbAngle(unittest.TestCase):
    def setUp(self):
        self.position_obj = (0.0, 0.0, 0.0)
        self.position_obs = (0.0, 0.0, 0.0)

    def test_position_obj_tuple_length(self):
        self.position_obj = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            limb, nadir = limb_angle(self.position_obj, self.position_obs)

    def test_position_obs_tuple_length(self):
        self.position_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            limb, nadir = limb_angle(self.position_obj, self.position_obs)


class TestMakeCatEntry(unittest.TestCase):
    def setUp(self):
        self.star_name = 'Test'
        self.catalog = 'XXX'
        self.star_num = 9999
        self.ra = 12.0
        self.dec = 0.0
        self.pm_ra = 0.0
        self.pm_dec = 0.0
        self.parallax = 100.0
        self.rad_vel = 0.0

    def test_ra_range_high(self):
        self.ra = 25.0

        with self.assertRaises(ValueError):
            star = make_cat_entry(self.star_name, self.catalog, self.star_num, self.ra, self.dec, self.pm_ra, self.pm_dec, self.parallax, self.rad_vel)

    def test_ra_range_low(self):
        self.ra = -1.0

        with self.assertRaises(ValueError):
            star = make_cat_entry(self.star_name, self.catalog, self.star_num, self.ra, self.dec, self.pm_ra, self.pm_dec, self.parallax, self.rad_vel)

    def test_dec_range_high(self):
        self.dec = 100.0

        with self.assertRaises(ValueError):
            star = make_cat_entry(self.star_name, self.catalog, self.star_num, self.ra, self.dec, self.pm_ra, self.pm_dec, self.parallax, self.rad_vel)

    def test_dec_range_low(self):
        self.dec = -100.0

        with self.assertRaises(ValueError):
            star = make_cat_entry(self.star_name, self.catalog, self.star_num, self.ra, self.dec, self.pm_ra, self.pm_dec, self.parallax, self.rad_vel)

    def test_negative_parallax(self):
        self.parallax = -100.0

        with self.assertRaises(ValueError):
            star = make_cat_entry(self.star_name, self.catalog, self.star_num, self.ra, self.dec, self.pm_ra, self.pm_dec, self.parallax, self.rad_vel)


class TestMakeObject(unittest.TestCase):
    def setUp(self):
        self.type = 0
        self.number = 1
        self.name = 'Test'
        self.star_data = None

    def test_type_value(self):
        self.type = 99

        with self.assertRaises(ValueError):
            object = make_object(self.type, self.number, self.name, self.star_data)


class TestMakeObserver(unittest.TestCase):
    pass


class TestMakeObserverAtGeocenter(unittest.TestCase):
    pass


class TestMakeObserverInSpace(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.velocity = (0.0, 0.0, 0.0)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            observer = make_observer_in_space(self.position, self.velocity)

    def test_velocity_tuple_length(self):
        self.velocity = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            observer = make_observer_in_space(self.position, self.velocity)


class TestMakeObserverOnSurface(unittest.TestCase):
    def setUp(self):
        self.latitude = 40.0
        self.longitude = -105.0
        self.height = 1800.0
        self.temperature = 20.0
        self.pressure = 1020.0

    def test_latitude_range_high(self):
        self.latitude = 100.0

        with self.assertRaises(ValueError):
            observer = make_observer_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_latitude_range_low(self):
        self.latitude = -100.0

        with self.assertRaises(ValueError):
            observer = make_observer_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_longitude_range_high(self):
        self.longitude = 360.0

        with self.assertRaises(ValueError):
            observer = make_observer_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_longitude_range_low(self):
        self.longitude = -360.0

        with self.assertRaises(ValueError):
            observer = make_observer_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)


class TestMakeInSpace(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)
        self.velocity = (0.0, 0.0, 0.0)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            object = make_in_space(self.position, self.velocity)

    def test_velocity_tuple_length(self):
        self.velocity = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            object = make_in_space(self.position, self.velocity)


class TestMakeOnSurface(unittest.TestCase):
    def setUp(self):
        self.latitude = 40.0
        self.longitude = -105.0
        self.height = 1800.0
        self.temperature = 20.0
        self.pressure = 1020.0

    def test_latitude_range_high(self):
        self.latitude = 100.0

        with self.assertRaises(ValueError):
            observer = make_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_latitude_range_low(self):
        self.latitude = -100.0

        with self.assertRaises(ValueError):
            observer = make_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_longitude_range_high(self):
        self.longitude = 360.0

        with self.assertRaises(ValueError):
            observer = make_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)

    def test_longitude_range_low(self):
        self.longitude = -360.0

        with self.assertRaises(ValueError):
            observer = make_on_surface(self.latitude, self.longitude, self.height, self.temperature, self.pressure)


class TestNutation(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.direction = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = nutation(self.date, self.position, self.direction)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = nutation(self.date, self.position, self.direction)

    def test_direction_value(self):
        self.direction = 99

        with self.assertRaises(ValueError):
            x, y, z = nutation(self.date, self.position, self.direction)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            x, y, z = nutation(self.date, self.position, self.direction, self.accuracy)


class TestPrecession(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.newdate = 2455619.5

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = precession(self.date, self.position, self.newdate)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = precession(self.date, self.position, self.newdate)

    def test_negative_julian_newdate(self):
        self.newdate = -2455619.5

        with self.assertRaises(ValueError):
            x, y, z = precession(self.date, self.position, self.newdate)


class TestProperMotion(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.position = (0.0, 0.0, 0.0)
        self.velocity = (0.0, 0.0, 0.0)
        self.newdate = 2455619.5

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = proper_motion(self.date, self.position, self.velocity, self.newdate)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = proper_motion(self.date, self.position, self.velocity, self.newdate)

    def test_velocity_tuple_length(self):
        self.velocity = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = proper_motion(self.date, self.position, self.velocity, self.newdate)

    def test_negative_julian_newdate(self):
        self.newdate = -2455619.5

        with self.assertRaises(ValueError):
            x, y, z = proper_motion(self.date, self.position, self.velocity, self.newdate)


class TestRadVel(unittest.TestCase):
    def setUp(self):
        self.object = make_object(0, 11, 'Moon', None)
        self.position = (0.0, 0.0, 0.0)
        self.velocity = (0.0, 0.0, 0.0)
        self.velocity_obs = (0.0, 0.0, 0.0)
        self.distance_geo = 0.1
        self.distance_sun = 1.0
        self.distance_obj = 1.1

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)

    def test_velocity_tuple_length(self):
        self.velocity = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)

    def test_velocity_obs_tuple_length(self):
        self.velocity_obs = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)

    def test_negative_distance_geo(self):
        self.distance_geo = -0.1

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)

    def test_negative_distance_sun(self):
        self.distance_sun = -0.1

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)

    def test_negative_distance_obj(self):
        self.distance_obj = -0.1

        with self.assertRaises(ValueError):
            vel = rad_vel(self.object, self.position, self.velocity, self.velocity_obs, self.distance_geo, self.distance_sun, self.distance_obj)


class TestRadec2Vector(unittest.TestCase):
    def setUp(self):
        self.ra = 12.0
        self.dec = 0.0
        self.distance = 1.0

    def test_ra_range_high(self):
        self.ra = 25.0

        with self.assertRaises(ValueError):
            x, y, z = radec2vector(self.ra, self.dec, self.distance)

    def test_ra_range_low(self):
        self.ra = -1.0

        with self.assertRaises(ValueError):
            x, y, z = radec2vector(self.ra, self.dec, self.distance)

    def test_dec_range_high(self):
        self.dec = 100.0

        with self.assertRaises(ValueError):
            x, y, z = radec2vector(self.ra, self.dec, self.distance)

    def test_dec_range_low(self):
        self.dec = -100.0

        with self.assertRaises(ValueError):
            x, y, z = radec2vector(self.ra, self.dec, self.distance)

    def test_negative_distance(self):
        self.distance = -1.0

        with self.assertRaises(ValueError):
            x, y, z = radec2vector(self.ra, self.dec, self.distance)


class TestSpin(unittest.TestCase):
    def setUp(self):
        self.angle = 180.0
        self.position = (0.0, 0.0, 0.0)

    def test_angle_range_high(self):
        self.angle = 450.0

        with self.assertRaises(ValueError):
            x, y, z = spin(self.angle, self.position)

    def test_angle_range_low(self):
        self.angle = -90.0

        with self.assertRaises(ValueError):
            x, y, z = spin(self.angle, self.position)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = spin(self.angle, self.position)


class TestStarvectors(unittest.TestCase):
    pass


class TestTerra(unittest.TestCase):
    def setUp(self):
        self.location = make_on_surface(40.0, -105.0, 1800.0, 20.0, 1020.0)
        self.time = 12.0

    def test_time_range_high(self):
        self.time = 25.0

        with self.assertRaises(ValueError):
            position, velocity = terra(self.location, self.time)

    def test_time_range_low(self):
        self.time = -1.0

        with self.assertRaises(ValueError):
            position, velocity = terra(self.location, self.time)


class TestVector2Radec(unittest.TestCase):
    def setUp(self):
        self.position = (0.0, 0.0, 0.0)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            rightascension, declination = vector2radec(self.position)


class TestWobble(unittest.TestCase):
    def setUp(self):
        self.date = 2455519.5
        self.xp = 0.0
        self.yp = 0.0
        self.position = (0.0, 0.0, 0.0)

    def test_negative_julian_date(self):
        self.date = -2455519.5

        with self.assertRaises(ValueError):
            x, y, z = wobble(self.date, self.xp, self.yp, self.position)

    def test_position_tuple_length(self):
        self.position = (0.0, 0.0, 0.0, 0.0)

        with self.assertRaises(ValueError):
            x, y, z = wobble(self.date, self.xp, self.yp, self.position)


class TestCalDate(unittest.TestCase):
    def setUp(self):
        self.jd = 2455519.5

    def test_negative_julian_date(self):
        self.jd = -2455519.5

        with self.assertRaises(ValueError):
            date = cal_date(self.jd)


class TestCelPole(unittest.TestCase):
    def setUp(self):
        self.jd = 2455519.5
        self.type = 1
        self.dpole1 = 0.0
        self.dpole2 = 0.0

    def test_negative_julian_date(self):
        self.jd = -2455519.5

        with self.assertRaises(ValueError):
            cel_pole(self.jd, self.type, self.dpole1, self.dpole2)

    def test_type_value(self):
        self.type = 99

        with self.assertRaises(ValueError):
            cel_pole(self.jd, self.type, self.dpole1, self.dpole2)


class TestCIOLocation(unittest.TestCase):
    def setUp(self):
        self.day = 2455519.5
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.day = -2455519.5

        with self.assertRaises(ValueError):
            ra_cio = cio_location(self.day)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra_cio = cio_location(self.day, self.accuracy)


class TestCIORa(unittest.TestCase):
    def setUp(self):
        self.day = 2455519.5
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.day = -2455519.5

        with self.assertRaises(ValueError):
            ra_cio = cio_ra(self.day)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ra_cio = cio_ra(self.day, self.accuracy)


class TestEECT(unittest.TestCase):
    def setUp(self):
        self.jd_high = 2455519.0
        self.jd_low = 0.5
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd_high = -2455519.0

        with self.assertRaises(ValueError):
            comp_terms = ee_ct(self.jd_high, self.jd_low)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            comp_terms = ee_ct(self.jd_high, self.jd_low, self.accuracy)


class TestEphemeris(unittest.TestCase):
    def setUp(self):
        self.jd = (2455519.0, 0.5)
        self.ss_body = make_object(0, 11, "Moon", None)
        self.origin = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd = (-2455519.0, 0.5)

        with self.assertRaises(ValueError):
            position, velocity = ephemeris(self.jd, self.ss_body, self.origin, self.accuracy)

    def test_origin_value(self):
        self.origin = 99

        with self.assertRaises(ValueError):
            position, velocity = ephemeris(self.jd, self.ss_body, self.origin, self.accuracy)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            position, velocity = ephemeris(self.jd, self.ss_body, self.origin, self.accuracy)


class TestETilt(unittest.TestCase):
    def setUp(self):
        self.jd = 2455519.5
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.jd = -2455519.5

        with self.assertRaises(ValueError):
            mobl, tobl, ee, dpsi, deps = e_tilt(self.jd)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            mobl, tobl, ee, dpsi, deps = e_tilt(self.jd, self.accuracy)


class TestFundArgs(unittest.TestCase):
    pass


class TestIRAEquinox(unittest.TestCase):
    def setUp(self):
        self.day = 2455519.5
        self.equinox = 0
        self.accuracy = 0

    def test_negative_julian_date(self):
        self.day = -2455519.5

        with self.assertRaises(ValueError):
            ira_eq = ira_equinox(self.day, self.equinox)

    def test_equinox_value(self):
        self.equinox = 99

        with self.assertRaises(ValueError):
            ira_eq = ira_equinox(self.day, self.equinox)

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            ira_eq = ira_equinox(self.day, self.equinox, self.accuracy)


class TestJulianDate(unittest.TestCase):
    def setUp(self):
        self.year = 2010
        self.month = 11
        self.day = 19
        self.hour = 0.0

    def test_year_lower_bound(self):
        self.year = -9999

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_month_range_high(self):
        self.month = 13

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_month_range_low(self):
        self.month = -1

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_day_range_high(self):
        self.day = 32

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_day_range_low(self):
        self.day = -1

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_hour_range_high(self):
        self.hour = 24.0

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)

    def test_hour_range_low(self):
        self.hour = -1.0

        with self.assertRaises(ValueError):
            jd = julian_date(self.year, self.month, self.day, self.hour)


class TestMeanObliq(unittest.TestCase):
    def setUp(self):
        self.day = 2455519.5

    def test_negative_julian_date(self):
        self.day = -2455519.5

        with self.assertRaises(ValueError):
            epsilon = mean_obliq(self.day)


class TestNormAng(unittest.TestCase):
    pass


class TestNutationAngles(unittest.TestCase):
    def setUp(self):
        self.t = 0.0
        self.accuracy = 0

    def test_accuracy_value(self):
        self.accuracy = 99

        with self.assertRaises(ValueError):
            dpsi, deps = nutation_angles(self.t, self.accuracy)


class TestRefract(unittest.TestCase):
    def setUp(self):
        self.location = make_on_surface(40.0, -105.0, 1800.0, 20.0, 1020.0)
        self.zenith = 89.0
        self.atmosphere = 1

    def test_atmosphere_value(self):
        self.atmosphere = 99

        with self.assertRaises(ValueError):
            refraction = refract(self.location, self.zenith, self.atmosphere)


class TestTDB2TT(unittest.TestCase):
    def setUp(self):
        self.day = 2455519.5

    def test_negative_julian_date(self):
        self.day = -2455519.5

        with self.assertRaises(ValueError):
            tt_jd, secdiff = tdb2tt(self.day)

if __name__ == '__main__':
    unittest.main()
