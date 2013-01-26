from unittest import TestCase

import novas.compat as novas

accuracy = 0

jd = 2454580.941875
ut1_utc = -0.387845
leap_secs = 33.0
delta_t = 32.184 + leap_secs - ut1_utc

latitude = 42.0
longitude = -70.0
height = 0.0
temperature = 10.0
pressure = 1010.0

class IssueTests(TestCase):

    def test_github_issue1_limb_angle(self):

        # https://github.com/brandon-rhodes/python-novas/pull/1
        # The _limb_angle() function was always return zero.

        observer = novas.make_observer_on_surface(
            latitude, longitude, height, temperature, pressure)
        position_observer, velocity_observer = novas.geo_posvel(
            jd, delta_t, observer, accuracy)

        position_obj = (-10.0, -0.3, 0.0)

        limb, nadir = novas.limb_angle(position_obj, position_observer)

        self.assertAlmostEqual(limb, -21.7796002325, 12)
        self.assertAlmostEqual(nadir, 0.758004441861, 12)

    def test_github_issue1_cal_day(self):

        # https://github.com/brandon-rhodes/python-novas/pull/1
        # The cal_date() routine was always raising ArgumentError

        year, month, day, hour = novas.cal_date(jd)
        self.assertEqual(year, 2008)
        self.assertEqual(month, 4)
        self.assertEqual(day, 24)
        self.assertAlmostEqual(hour, 10.605000000447035, 15)
