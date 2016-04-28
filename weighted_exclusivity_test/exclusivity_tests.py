#!/usr/bin/env python

# Load required modules
import numpy as np
from constants import *
from exact import exact_test, py_exact_test
import cpoibin
from saddlepoint import saddlepoint
import comet_driver as comet

# Perform the weighted exclusivity test using the given method.
def weighted_test(t, x, p, method=EXACT, tail=ONE_GREATER, check=True, verbose=0):
    if check:
        # Check we're using an appropriate tail and method
        assert( method in METHODS and tail in TAILS )
        # Check for equal numbers of genes.
        assert(len(x)==len(p))
        # Check for equal numbers of samples.
        assert(min(len(a) for a in p)==max(len(a) for a in p))
        # Check that the probabilities are in (0, 1].
        assert(all(0<b<= 1 for a in p for b in a))
        # Check that the number of mutations in each gene is not greater than the number of samples.
        assert(all(a<=len(b) for a, b in zip(x, p)))
        # Check that the number of mutually exclusive mutations is not greater than the total number of mutations.
        assert(t<=sum(x))

        p = [ list(p_g) for p_g in p ]
        
    if method == EXACT:
        p_value = exact_test( t, x, p, tail, verbose )
    if method == SADDLEPOINT:
        p_value = saddlepoint( t, x, p, tail, verbose )

    return p_value

# Perform the unweighted test
def unweighted_test(t, x, tbl=None, method=EXACT, tail=ONE_GREATER, verbose=0):
    if tbl: N = sum(tbl)
    if method == SADDLEPOINT:
        p = [ [1./x_i] * N for x_i in x ]
        p_value = saddlepoint( t, x, p, tail, verbose )
    elif method == EXACT:
        k = len(x)
        assert( tbl and len(tbl) == 2**k )
        if not comet.loaded:
            raise NotImplementedError("CoMEt is not available to compute the unweighted test exactly")
        else:
            k = len(x)
            if N >= cometMaxN:
                raise NotImplementedError("CoMEt is only initalized to compute P-values for N < {}".format(comet.maxN))
            num_tbls, p_value = comet.exact_test( k, N, tbl, 1.1 )

    return p_value
