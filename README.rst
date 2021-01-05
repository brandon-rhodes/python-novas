###################
Welcome to NOVAS Py
###################

.. important::

   This GitHub repository is merely a mirror of this official United
   States Naval Observatory software library, with the code tweaked to
   make possible its release on the Python Package Index as the `novas
   package`_.  None of the library’s actual authors are involved with
   this repository.  Its official documentation is here:

   * `User’s Guide to NOVAS`_

What is NOVAS?
==============

NOVAS is an integrated package of functions for computing various
commonly needed quantities in positional astronomy. The package can
supply, in one or two function calls, the instantaneous coordinates of
any star or solar system body in a variety of coordinate systems.  At a
lower level, NOVAS also provides astrometric utility transformations,
such as those for precession_, nutation_, aberration_, parallax_, and
gravitational `deflection of light`_.  The computations are accurate to
better than one milliarcsecond. The NOVAS library is an easy-to-use
facility that can be incorporated into data reduction programs,
telescope control systems, and simulations.  The U.S. parts of
*The Astronomical Almanac* are prepared using NOVAS.

What is NOVAS Py?
=================

With NOVAS Py, the USNO is expanding NOVAS to the Python programming
language. The NOVAS Py module is simply a wrapper around the NOVAS C code;
*all computations are still performed by the C code*. NOVAS Py makes use of
Python's `ctypes`_ module.

Installation
============

The NOVAS Py Module
-------------------

The NOVAS Py module is installed from the top-level source directory
with the command ``python setup.py install``. If there are multiple
versions of Python installed, the NOVAS Py package will be installed for
the version used to run the ``setup.py`` script. Note that you may need
superuser or administrator privileges to install the package.

Ephemerides
-----------

NOVAS requires access to a high-accuracy solar system ephemeris in order
to compute places of solar system bodies and the highest-accuracy star
places. Groups in the U.S., France, and Russia now construct
high-accuracy solar system ephemerides. NOVAS is able to use solar system
ephemerides that use JPL's export format, for example, the "developmental
ephemerides," designated as "DEnnn", which are produced by JPL in the
U.S. NOVAS C provides an implementation of JPL's ephemeris-access
software that enables reading and interpolating a binary, direct-access
ephemeris file.

Before using any NOVAS functionality that requires access to the
ephemerides, you must first open the file with the ``ephem_open``
function (``from novas.compat.eph_manager import ephem_open``). You must
either pass ``ephem_open`` the path to the binary ephemeris file you
wish to use or have the path set in an environment variable named
``EPHEMERIS_FILE``. If you choose to do the latter, you may then call
``ephem_open`` without passing any arguments to the function.

Consult the NOVAS C User Guide for directions on how to create the binary
ephemeris file.

Package Layout
--------------

Phase one of NOVAS Py is intended to feel like NOVAS C. For the most
part, function calls in Python match the function calls in C. All
results are returned from the function. Some function inputs have been
reordered so that the function can support optional inputs.

NOVAS functions can be found under the ``novas.compat`` namespace.
Functions from ``eph_manager.c``, ``solsys1.c``, and ``nutation.c`` can
be found under ``novas.compat.eph_manager``, ``novas.compat.solsys``,
and ``novas.compat.nutation``, respectively.

Note on constants
=================

NOVAS Py includes a constants file copied from the NOVAS C constants file; this
is only provided for consistency with the NOVAS C package. Since the wrapped
NOVAS C functions are where the constants are used, and those functions obtain
the constants from the NOVAS C constants file, any changes to constants in the
NOVAS Py constants file will have no effect.

Notes on error handling
=======================

Error return codes from NOVAS C functions are handled by raising an
appropriate exception. However, given the way NOVAS C does error codes,
it was found that some return codes from higher-level C functions can
represent more than one error. For these cases the developer choose what was
felt to be the most likely or most important error code to escalate to an
exception. All possible error codes are in the ``c_error`` dictionaries
attached to each C function object, within each wrapper function, but those
with duplicate numbers are commented out at this time.

Tests
=====

Some tests are available in the ``tests`` source directory; they are
designed to work with Python >= 2.7. To run the tests with Python 2.5 or
2.6, first install the `unittest2 module`_

Using NOVAS Py
==============

Once installation is complete, the NOVAS functions can be found under
the ``novas.compat`` namespace. Nutation models can be found under
``novas.nutation``, and constants under ``novas.constants``.

.. _novas package: https://pypi.org/project/novas/
.. _User’s Guide to NOVAS: https://github.com/brandon-rhodes/python-novas/raw/master/Cdist/NOVAS_C3.1_Guide.pdf
.. _precession: http://asa.usno.navy.mil/SecM/Glossary.html#precession
.. _nutation: http://asa.usno.navy.mil/SecM/Glossary.html#nutation
.. _aberration: http://asa.usno.navy.mil/SecM/Glossary.html#aberration
.. _parallax: http://asa.usno.navy.mil/SecM/Glossary.html#parallax
.. _deflection of light: http://asa.usno.navy.mil/SecM/Glossary.html#deflection-light
.. _ctypes: https://docs.python.org/3.4/library/ctypes.html
.. _webpage: http://ssd.jpl.nasa.gov/?planet_eph_export
.. _unittest2 module: https://pypi.org/project/unittest2/
