#!/usr/bin/env python

# Load required modules
import numpy as np
from constants import *
from exact import exact_test
import cpoibin
from saddlepoint import saddlepoint
from comet_exact_test import comet_exact_test
import warnings

# Check weighted test arguments
def weighted_check(t, x, p):
    # Check for equal numbers of genes.
    if len(x) != len(p): return False
    # Check for equal numbers of samples.
    elif min(len(a) for a in p) != max(len(a) for a in p): return False
    # Check that the probabilities are in (0, 1].
    elif not all(0<b<= 1 for a in p for b in a): return False
    # Check that the number of mutations in each gene is not greater than the number of samples.
    elif not all(a<=len(b) for a, b in zip(x, p)): return False
    # Check that the number of mutually exclusive mutations is not greater than the total number of mutations.
    elif t> sum(x): return False
    else: return True

# Subtype test
def subtype_test(s, x, p=None, N=None, method=SADDLEPOINT, verbose=0):
    # If p isn't provided, we make dummy weights
    assert( p or N )
    if not p:
        p = [ [ float(x_i)/N ] * N for x_i in x ]

    # Check arguments
    assert( method == SADDLEPOINT )
    assert( weighted_check(s, x, p) )

    # Ignore warnings
    with warnings.catch_warnings() as e:
        warnings.simplefilter("ignore")
        p_value = saddlepoint( s, x, p, condition=SUBTYPE, verbose=verbose )

    return p_value

# Perform the weighted-row exclusivity test (WR-test) using the given method.
# Note that EXACT refers to the WR-exclusivity recursive formula, and computes
# the p-value _exactly_.
def wre_test(t, x, p, method=EXACT, verbose=0):
    # Check we're using an appropriate method and that the arguments are as
    # expected
    assert( method in METHODS )
    assert( weighted_check(t, x, p) )

    #Check that we've implemented the given set size with the exact test
    if method == EXACT:
        assert( len(x) in WRE_EXACT_SET_SIZES_IMPLEMENTED )

    p = [ list(p_g) for p_g in p ]

    if method == EXACT:
        p_value = exact_test( t, x, p, condition=EXCLUSIVITY, verbose=verbose )
    if method == SADDLEPOINT:
        # Ignore warnings
        with warnings.catch_warnings() as e:
            warnings.simplefilter("ignore")
            p_value = saddlepoint( t, x, p, condition=EXCLUSIVITY, verbose=verbose )

    return p_value

# Perform the row-exclusivity test (RE-test) using the given method.
# Note that EXACT refers to the CoMEt tail enumeration scheme, and computes
# the p-value _exactly_.
def re_test(t, x, tbl, method=EXACT, verbose=0):
    N = sum(tbl)
    if method == SADDLEPOINT:
        p = [ [ float(x_i)/N ] * N for x_i in x ]
        # Ignore warnings
        with warnings.catch_warnings() as e:
            warnings.simplefilter("ignore")
            p_value = saddlepoint( t, x, p, condition=EXCLUSIVITY, verbose=verbose )

    elif method == EXACT:
        k = len(x)
        assert( tbl and len(tbl) == 2**k )
        p_value, mid_p_value = comet_exact_test( k, N, tbl, 1.1 )

    return p_value
