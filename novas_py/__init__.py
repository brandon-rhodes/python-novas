# -*- coding: utf-8 -*-

import os
import sys
from ctypes import CDLL

# By Brandon Rhodes, for Python 3 compatibility:
def load_shared_library():
    global novaslib

    names = ['libnovas.so']
    if sys.version_info >= (3, 2):
        import sysconfig
        name = 'libnovas.{}.so'.format(sysconfig.get_config_var('SOABI'))
        names.append(name)
        if name.endswith('-33dm.so'):
            names.append(name.replace('-33dm.so', '-33m.so'))

    for directory in sys.path:
        for name in names:
            try:
                novaslib = CDLL(os.path.join(directory, 'novas', name))
                break
            except:
                pass

load_shared_library()
del load_shared_library
