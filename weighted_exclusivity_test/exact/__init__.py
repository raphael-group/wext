#!/usr/bin/env python

import numpy as np
from .. import cexact_test
from .. import cpoibin
from ..constants import *
conditional = cexact_test.conditional

def exact_test( t, x, p, tail=ONE_LESS, verbose=False ):
    k = len(x)
    if k == 2:
        return exact_test_k2( t, x, p, tail, verbose )
    elif k == 3:
        return exact_test_k3( t, x, p, tail, verbose )
    else:
        raise NotImplementedError("Exact test: only k=2,3 implemented.")

def py_exact_test( t, x, p, tail=ONE_LESS, verbose=False ):
    k = len(x)
    if k == 2:
        # Two-sided test
        N = len(p[0])
        x, y, z, p_x, p_y = x[0], x[1], (sum(x) - t)/2, p[0], p[1]
        if tail == TWO:
            tail_masses = conditional(N, range(N+1), x, y, p_x, p_y)
            obs_mass  = tail_masses[z]
            pval = sum( m for m in tail_masses if m <= obs_mass )
        # One-sided (less) test
        elif tail == ONE_LESS:
            tail_masses = conditional(N, range(z+1), x, y, p_x, p_y)
            obs_mass  = tail_masses[-1]
            pval = sum(tail_masses)
        # One-sided (greater) test
        elif tail == ONE_GREATER:
            tail_masses = conditional(N, range(z, N+1), x, y, p_x, p_y)
            obs_mass  = tail_masses[0]
            pval = np.array(tail_masses).sum()
        return pval
    elif k == 3:
        return py_exact_test_k3( t, x, p, tail, verbose )
    else:
        raise NotImplementedError("Exact test: only k=2,3 implemented.")

# Compute the exact test, reporting the given tail
def exact_test_k2(t, (x, y), (p_x, p_y), tail, verbose):
	# Two-sided test
    N = len(p_x)
    z = (x + y - t)/2 # count number of co-occurrences
    if tail == TWO:
		tail_masses = conditional(N, range(N+1), x, y, p_x, p_y)
		obs_mass  = tail_masses[z]
		pval = sum( m for m in tail_masses if m <= obs_mass )
	# One-sided (less exclusivity == greater co-occurrence) test
    elif tail == ONE_GREATER:
		tail_masses = conditional(N, range(z+1), x, y, p_x, p_y)
		obs_mass  = tail_masses[-1]
		pval = sum(tail_masses)
	# One-sided (greater exclusivity == less co-occurrence) test
    elif tail == ONE_LESS:
		tail_masses = conditional(N, range(z, N+1), x, y, p_x, p_y)
		obs_mass  = tail_masses[0]
		pval = np.array(tail_masses).sum()
    return pval

# Exact test k=3 just calls C directly, but only has the ONE_GREATER tail
# implemented
def exact_test_k3(t, x, p, tail, verbose):
    if tail != ONE_GREATER:
        raise NotImplementedError("Exact test: Only less tail implemented for k=3")
    N = len(p[0])
    return cexact_test.triple_exact_test( N, t, x[0], x[1], x[2], p )

# PMF for a weighted version of the hypergeometric (see http://goo.gl/QLDTUF)
def py_joint_mass( n, z, x, y, p_x, p_y, cache, verbose=False ):
	if verbose: print 'n={}, l={}, j={}, r={}'.format(n, z, x, y)
	# Use the cache if possible
	if (n, z, x, y) in cache: pass
	# Base cases (top of page 7)
	elif 0 > min([z, x, y]) or z > min([x, y]) or max([x, y]) > n:
		cache[(n, z, x, y)] = 0.
	elif n == z == x == y == 0:
		cache[(n, z, x, y)] = 1.
	# Recursive case
	else:
		# Equation at bottom of page 6
		i = n-1 # subtract one because of zero-based indexing
		val  = p_x[i]      * p_y[i]      * py_joint_mass(n-1, z-1, x-1, y-1, p_x, p_y, cache) + \
			   p_x[i]      * (1.-p_y[i]) * py_joint_mass(n-1, z, x-1, y, p_x, p_y, cache) + \
			   (1.-p_x[i]) * p_y[i]      * py_joint_mass(n-1, z, x, y-1, p_x, p_y, cache) + \
			   (1.-p_x[i]) * (1.-p_y[i]) * py_joint_mass(n-1, z, x, y, p_x, p_y, cache)

		cache[(n, z, x, y)] = val

	return cache[(n, z, x, y)]

# Compute the conditional distribution
def py_conditional( N, zs, x, y, p_x, p_y, cache=dict(), verbose=False):
	x_prob = cpoibin.pmf(x, p_x)
	y_prob = cpoibin.pmf(y, p_y)
	joints  = [ py_joint_mass(N, z, x, y, p_x, p_y, cache, verbose) for z in zs ]
	if verbose: print x_prob, y_prob, joints
	return  [ joint / (x_prob * y_prob) for joint in joints ]


# Define the triple exact test in Python
def py_exact_test_k3(T, X, p, tail, verbose):
    if tail != ONE_GREATER: raise ValueError("Exact test (Python): Only implemented for greater tail.")
    def conditional_k3(t, n, w, x, y, cache):
        joint = P(t, n, w, x, y, cache)
        marginals = cpoibin.pmf(w, p[0]) * cpoibin.pmf(x, p[1]) * cpoibin.pmf(y, p[2])
        return joint / marginals

    def P(t, n, w, x, y, cache):
        if (t, n, w, x, y) not in cache:
            if 0 > min([t, w, x, y]) or t > sum([w, x, y]) or max([w, x, y]) > n:
                cache[(t, n, w, x, y)] = 0.
            elif t == n == w == x == y == 0:
                cache[(t, n, w, x, y)] = 1.
            else:
                i = n - 1
                cache[(t, n, w, x, y)] = \
                  (1. - p[0][i]) * (1. - p[1][i]) * (1. - p[2][i]) * P(t,   n-1, w,   x,   y,   cache) + \
                  (1. - p[0][i]) * (1. - p[1][i]) * p[2][i]        * P(t-1, n-1, w,   x,   y-1, cache) + \
                  (1. - p[0][i]) * p[1][i]        * (1. - p[2][i]) * P(t-1, n-1, w,   x-1, y,   cache) + \
                  (1. - p[0][i]) * p[1][i]        * p[2][i]        * P(t,   n-1, w,   x-1, y-1, cache) + \
                  p[0][i]        * (1. - p[1][i]) * (1. - p[2][i]) * P(t-1, n-1, w-1, x,   y,   cache) + \
                  p[0][i]        * (1. - p[1][i]) * p[2][i]        * P(t,   n-1, w-1, x,   y-1, cache) + \
                  p[0][i]        * p[1][i]        * (1. - p[2][i]) * P(t,   n-1, w-1, x-1, y,   cache) + \
                  p[0][i]        * p[1][i]        * p[2][i]        * P(t,   n-1, w-1, x-1, y-1, cache)
        return cache[(t, n, w, x, y)]

    N = len(p[0])
    cache = dict()
    masses = [ conditional_k3(t, N, X[0], X[1], X[2], cache) for t in range(T, min([N, sum(X)])+1) ]

    return sum(masses)

