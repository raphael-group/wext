#!/usr/bin/env python

# PMF implementation using recursive formula from equation (9) in
# Hong [Computational Statistics & Data Analysis, 2013 ](http://goo.gl/m7smXj)
def pmf_recursion(k, j, ps, cache):
    if (k, j) in cache:
        return cache[(k, j)]
    # Boundary conditions
    elif k == 0 and j == 0:
        val = 1.
    elif k == -1 or k == j + 1:
        val = 0.
    # Recursion
    else:
        val = (1.-ps[j-1]) * pmf_recursion(k, j-1, ps, cache) + ps[j-1] * pmf_recursion(k-1, j-1, ps, cache)

    cache[(k, j)] = val
    return cache[(k, j)]

# Wrapper for PMF so we don't have to provide the length as
# input and can check the arguments
def pmf(k, ps):
    N = len(ps)
    assert( k <= N )
    assert( all(0 <= p <= 1) for p in ps )
    return pmf_recursion(k, N, ps, dict())
