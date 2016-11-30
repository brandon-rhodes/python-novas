# -*- coding: utf-8 -*-

import os
import sys
from ctypes import CDLL


# Exception strings used when checking inputs
_neg_err = "{name} must be >= 0.0"
_hour_range_err = "{name} must be in the range 0.0 <= {name} < 24.0"
_elev_range_err = "{name} must be in the range -90.0 <= {name} <= 90.0"
_az180_range_err = "{name} must be in the range -180.0 <= {name} <= 180.0"
_az360_range_err = "{name} must be in the range 0.0 <= {name} <= 360.0"
_vector_len_err = "{name} must be a sequence of length 3"
_option_err = "{name} must be in {allowed}"
_jd_year_err = "{name} must be >= -4712"
_month_range_err = "{name} must be 1 <= {name} <= 12"
_day_range_err = "{name} must be 1 <= {name} <= {daysinmonth}"

_days_in_month = (None, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


class NonConvergentError(Exception):
    pass


class IndeterminateError(Exception):
    pass


class InitializationError(Exception):
    pass


def _check_c_errors(retval, func, args):
    """
    Function to check return values from wrapped C functions.

    This function is used in conjunction with the ctypes ``errcheck`` attribute
    attached to a C function object. Assigning this function to the
    ``errcheck`` attribute causes an automatic check of the return value from
    the C function. Every C function object that utilizes this ``errcheck``
    feature must also have an attribute called ``c_errors`` assigned to it.
    ``c_errors`` is a dictionary where the keys are possible return values and
    the values are length-2 tuples; each tuple consists of the exception to be
    raised and the message sent to that exception.

    https://docs.python.org/3.4/library/ctypes.html#ctypes._FuncPtr.errcheck

    """
    if retval:
        error = func.c_errors[retval]
        raise error[0](error[1])


def get_and_configure_library():
    for directory in sys.path:
        path_to_libnovas = os.path.join(directory, 'novas', 'libnovas.so')
        if os.path.isfile(path_to_libnovas):
            libnovas = CDLL(path_to_libnovas)

    return libnovas

novaslib = get_and_configure_library()
