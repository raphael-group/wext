#!/usr/bin/env python

import numpy as np
import wext_exact_test
from constants import *

def exact_test( t, x, p, verbose=False ):
    k = len(x)
    if k == 2:
        return exact_test_k2( t, x, p, verbose )
    elif k == 3:
        return exact_test_k3( t, x, p, verbose )
    else:
        raise NotImplementedError("Exact test: only k=2,3 implemented.")

# Wrapper for k=3 exact test C function
def exact_test_k3(t, x, p, verbose):
    N = len(p[0])
    return wext_exact_test.triple_exact_test( N, t, x[0], x[1], x[2], p )

# Wrapper for k=2 exact test C function
def exact_test_k2(t, (x, y), (p_x, p_y), verbose):
	# Two-sided test
    N = len(p_x)
    z = (x + y - t)/2 # count number of co-occurrences
    tail_masses = wext_exact_test.conditional(N, range(z+1), x, y, p_x, p_y)
    obs_mass  = tail_masses[-1]
    pval = sum(tail_masses)
    return pval

