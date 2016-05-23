#!/usr/bin/python

import sys, os, numpy as np
from collections import defaultdict
from time import time
from random import random, sample, choice, seed as random_seed

from constants import *
from enumerate_sets import observed_values
from exclusivity_tests import unweighted_test, weighted_test

@profile
def mcmc(ks, geneToCases, num_patients, method, test, geneToP, seed, annotations=set(), verbose=0, step_len=100, nchains=1, niters=1000, alpha=1):
    if verbose > 0: print '-' * 33, 'Running MCMC', '-' * 33

    # Set up a local version of the weight function
    if test == WEIGHTED:
        def _test(M, X, T, Z, tbl):
            return weighted_test(T, X, [ geneToP[g] for g in M ], method=method)
    elif test == UNWEIGHTED:
        def _test(M, X, T, Z, tbl):
            return unweighted_test(T, X, tbl, method=method)
    else:
        raise NotImplementedError('Test "{}" not implemented with MCMC'.format(testToName[test]))

    def _weight(M):
        if M not in setToPval:
            setToPval[M] = _test(M, *setToObs[M])
        return -np.log10(setToPval[M]**alpha)

    def _valid_set(M):
        # Compute or retrieve the observed statistics
        M = frozenset(M)
        if M not in setToObs:
            setToObs[M] = observed_values(M, num_patients, geneToCases)
        X, T, Z, tbl = setToObs[M]

        # We don't allow sets with T <= Z or with multiple annotations
        # (since these are often trivially exclusive)
        if T <= Z or len(M & annotations) > 1:
            return False
        else:
            return True

    def _collection_weight(collection):
        return sum( _weight(M) for M in collection )

    def _to_collection(solution):
        return frozenset( frozenset(M) for M in solution.values() )

    # Compute the acceptance ratio
    def _log_accept_ratio( W_current, W_next ): return W_next - W_current

    # Set up PRNG, sample space, and output
    random_seed(seed)
    t          = len(ks)
    genespace  = geneToCases.keys()
    setsToFreq = [ defaultdict(int) for _ in xrange(nchains) ]
    setToPval, setToObs =  dict(), dict()
    for c in xrange(nchains):
        if verbose > 0: print '- Experiment', c+1

        # Seed Markov chain
        soln, assigned = choose_random_set(ks, genespace)
        while not all( _valid_set(M) for M in _to_collection(soln) ):
            soln, assigned = choose_random_set(ks, genespace)
        weight = _collection_weight( _to_collection(soln) )

        # Run MCMC
        progress_step = float(np.ceil(niters / 72.))
        itera = 0
        while itera < niters:
            # Simple progress bar
            if verbose > 0 and (itera % progress_step == 0):
                sys.stdout.write("\r[%-72s] %d%%" % ('='*np.ceil(itera / progress_step) + '>', np.ceil(100.*itera / niters)))
                sys.stdout.flush()

            # Sample the next gene to swap in/around the set
            next_soln = dict( (index, set(M)) for index, M in soln.iteritems() )
            next_assigned = dict(assigned.items())
            next_gene = choice(genespace)

            # There are two possibilities for the next gene
            # 1) The gene we sampled is already in the current solution. In this
            #    case, we swap the gene with another gene in a different set.
            if next_gene in next_assigned:
                # if we only have one set, we can't swap between sets
                if t == 1: continue
                i = next_assigned[next_gene]
                swap_gene = choice([ g for g in next_assigned.keys() if g not in next_soln[i] ])
                j = next_assigned[swap_gene]
                next_assigned[swap_gene] = i
                next_soln[i].add(swap_gene)
                next_soln[i].remove( next_gene )
                next_assigned[next_gene] = j
                next_soln[j].remove(swap_gene)
                next_soln[j].add(next_gene)

                if not (_valid_set(next_soln[j]) and _valid_set(next_soln[i])):
                    continue

            # 2) The gene is not in the current solution. In this case, we choose
            #    a random gene in the solution to remove, and add the next gene.
            else:
                swap_gene = choice(next_assigned.keys())
                j = next_assigned[swap_gene]
                del next_assigned[swap_gene]
                next_assigned[next_gene] = j
                next_soln[j].remove(swap_gene)
                next_soln[j].add(next_gene)

                if not _valid_set(next_soln[j]): continue

            # Compare the current soln to the next soln
            next_weight = _collection_weight(_to_collection(next_soln))
            if _log_accept_ratio( weight, next_weight ) >= np.log10(random()):
                soln, weight, assigned = next_soln, next_weight, next_assigned

            # Freeze sets after thinning the chain a certain number of iterations
            itera += 1
            if (itera+1) % step_len == 0:
                setsToFreq[c][_to_collection(soln)] += 1

        if verbose > 0:
            print '\r[' + ('='*71) + '>] 100%'

    # Merge the various chains
    setsToTotalFreq = defaultdict(int)
    for counter in setsToFreq:
        for sets, freq in counter.iteritems():
            setsToTotalFreq[sets] += freq

    return setsToTotalFreq, setToPval, setToObs

# Choose a random set to initialize with
def choose_random_set(ks,  genespace):
    num_sampled_genes, t = sum(ks), len(ks)
    initial_genes = sample(genespace, num_sampled_genes)
    soln, assigned = dict(), dict()
    for i in range(t):
        soln[i] = set(initial_genes[sum(ks[:i]):sum(ks[:i+1])])
        for g in soln[i]: assigned[g] = i
    return soln, assigned
