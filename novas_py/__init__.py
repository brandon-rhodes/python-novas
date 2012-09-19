# -*- coding: utf-8 -*-

import os
import sys
from ctypes import CDLL

for directory in sys.path:
	try:
		novaslib = CDLL(os.path.join(directory, 'novas', 'libnovas.so'))
		break
	except:
		pass