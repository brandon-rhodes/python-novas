# -*- coding: utf-8 -*-

import os
import sys
from ctypes import CDLL

# By Brandon Rhodes, for Python 3 compatibility:
if sys.version_info >= (3, 2):
	import sysconfig
	soabi = sysconfig.get_config_var('SOABI')
	libname = 'libnovas.{}.so'.format(soabi)
else:
	libname = 'libnovas.so'

for directory in sys.path:
	try:
		novaslib = CDLL(os.path.join(directory, 'novas', libname))
		break
	except:
		pass