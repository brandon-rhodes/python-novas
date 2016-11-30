# -*- coding: utf-8 -*-

from novas.compat import topo_star, CatEntry, OnSurface
from novas.compat.eph_manager import ephem_open

jd_begin, jd_end, de_number = ephem_open()

print("JPL Ephemeris DE{de_number} open. jd_beg = {jd_begin:9.2f}  jd_end = {jd_end:9.2f}\n".format(de_number=de_number, jd_begin=jd_begin, jd_end=jd_end))

deltat = 60.0

tjds = [2450203.5, 2450203.5, 2450417.5, 2450300.5]

stars = [
    CatEntry("POLARIS".encode(), "HIP".encode(), 0, 2.530301028, 89.264109444, 44.22, -11.75, 7.56, -17.4),
    CatEntry("Delta ORI".encode(), "HIP".encode(), 1, 5.533444639, -0.299091944, 1.67, 0.56, 3.56, 16.0),
    CatEntry("Theta CAR".encode(), "HIP".encode(), 2, 10.715944806, -64.394450000, -18.87, 12.06, 7.43, 24.0)
]

geo_loc = OnSurface(45.0, -75.0, 0.0, 10.0, 1010.0)

for date in tjds:
    for star in stars:
        ra, dec = topo_star(date, deltat, star, geo_loc)

        print("JD = {jd:14.6f}  Star = {star}".format(jd=date, star=star.starname.decode('utf-8')))
        print("RA = {ra:>12.9f}  Dec = {dec: 12.8f}\n".format(ra=ra, dec=dec))

    print()
