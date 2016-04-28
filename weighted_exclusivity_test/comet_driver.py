#!/usr/bin/env python

# Load required modules
import sys, os
from constants import cometMaxN

# Try and load CoMEt
this_dir = os.path.dirname(os.path.realpath(__file__))
default_comet_dir = os.path.normpath(this_dir  + '/../third-party/comet')
try:
    sys.path.append(default_comet_dir)
    import comet as C
    loaded = True
    C.precompute_factorials(cometMaxN)
    exact_test = C.exact_test
except ImportError:
    loaded    = False
    exact_test = None
