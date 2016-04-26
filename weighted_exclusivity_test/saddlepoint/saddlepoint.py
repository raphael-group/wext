#!/usr/bin/env python

import numpy as np
from numpy.linalg import det
from scipy.optimize import fsolve
from scipy.stats import norm
import itertools

def condition(state):

    # Define a success condition for each state, e.g., mutual exclusivity.

    if sum(state)==1:
        return True
    else:
        return False

def enumeration(k):

    # Enumerate the states for the observed variables and identify the indices for the states that
    # satisfy the given condition, e.g., mutual exclusivity.

    states = list(itertools.product([0, 1], repeat=k))

    indices = []
    for i, state in enumerate(states):
        a = [j for j, s in enumerate(state) if s==1]
        if condition(state):
             a += [k]
        indices.append(a)

    # Identify the indices for each term of the gradient vector and the Hessian matrix.

    gradient_indices = []
    for i in range(k+1):
        b = [j for j, a in enumerate(indices) if i in a]
        gradient_indices.append(b)

    hessian_indices = []
    for i in range(k+1):
        b = []
        for j in range(k+1):
            c = [l for l, a in enumerate(indices) if i in a and j in a]
            b.append(c)
        hessian_indices.append(b)

    return states, indices, gradient_indices, hessian_indices

def saddlepoint(observed_t, observed_y, probabilities, enumeration_k=None):

    # Find the dimensions of the observations.

    k, n = np.shape(probabilities)

    # Enumerate the states for the observed variables and identify indices for the terms.

    if enumeration_k:
        states, indices, gradient_indices, hessian_indices = enumeration_k
    else:
        states, indices, gradient_indices, hessian_indices = enumeration(k)

    # Collect the observations and perform the continuity correction for t.

    y = np.zeros(k+1)
    y[0:k] = observed_y
    y[k] = observed_t-0.5

    # Precompute the products of the success and failure probabilities.

    p = np.zeros((2, k, n))
    p[0] = 1.0-probabilities
    p[1] = probabilities

    w = np.zeros((2**k, n))
    for i, state in enumerate(states):
        w[i, :] = np.product(p[state, range(k), :], axis=0)

    # Define the moment generating functions and cumulant generating functions.  These functions
    # use the above constants.

    def compute_terms(x):
        terms = np.zeros((2**k, n))
        for i, s in enumerate(indices):
            terms[i, :] = np.exp(np.sum(x[s]))*w[i, :]
        return terms

    # These are the joint moment and cumulant generating functions and their derivatives.

    def L(x):
        return np.sum(compute_terms(x), axis=0)

    def K(x):
        return np.sum(np.log(L(x)))

    def dK(x):
        terms = compute_terms(x)
        a = np.zeros(n)
        b = 1.0/np.sum(terms, axis=0)

        gradient_terms = np.zeros(k+1)

        for i in range(k+1):
            a[:] = np.sum(terms[gradient_indices[i], :], axis=0)
            gradient_terms[i] = np.dot(a, b)

        return gradient_terms

    def d2K(x):
        terms = compute_terms(x)

        a = np.zeros(n)
        b = 1.0/np.sum(terms, axis=0)
        c = np.zeros((k+1, n))

        for i in range(k+1):
            c[i, :] = np.sum(terms[gradient_indices[i], :], axis=0)*b

        hessian_terms = np.zeros((k+1, k+1))

        for i in range(k+1):
            for j in range(i, k+1):
                a[:] = np.sum(terms[hessian_indices[i][j], :], axis=0)
                hessian_terms[i, j] = np.sum(a*b - c[i]*c[j])
                if i!=j:
                    hessian_terms[j, i] = hessian_terms[i, j]

        return hessian_terms

    # These are the non-joint moment and cumulant generating functions and their derivatives.

    def Li(x, i):
        return np.exp(x)*p[1,i,:] + p[0,i,:]

    def Ki(x, i):
        return np.sum(np.log(Li(x, i)))

    def dKi(x, i):
        a = np.exp(x)*p[1,i,:]
        b = a + p[0,i,:]
        return np.array([np.sum(a/b)])

    def d2Ki(x, i):
        a = np.exp(x)*p[1,i,:]
        b = a + p[0,i,:]
        c = a/b

        return np.array([np.sum(c-c**2)])

    # Solve dKi(x_h) = y.

    x_h = np.zeros(k+1)

    for i in range(k):
        x_h[i] = fsolve(lambda a: dKi(a, i)-y[i], np.ones(1), fprime=lambda b: d2Ki(b, i))[0]

    # Compute Ki(x_h) and d2Ki(x_h).

    Ki_xh = np.zeros(k)
    d2Ki_xh = np.zeros(k)

    for i in range(k):
        Ki_xh[i] = Ki(x_h[i], i)
        d2Ki_xh[i] = d2Ki(x_h[i], i)

    # Solve dK(x_t) = y.

    x_t = fsolve(lambda a: dK(a)-y, np.ones(k+1), fprime=d2K)

    # Compute the saddlepoint approximation.

    w_t_inside = 2.0*(np.sum(Ki_xh)-K(x_t)-np.dot(y, x_h-x_t))
    w_t = np.sign(x_t[k])*np.sqrt(w_t_inside)

    u_t_inside = det(d2K(x_t))/np.product(d2Ki_xh)
    u_t = 2.0*np.sinh(0.5*x_t[k])*np.sqrt(u_t_inside)

    return norm.cdf(-w_t)-norm.pdf(w_t)*(1.0/w_t-1.0/u_t)
