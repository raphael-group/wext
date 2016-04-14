#!/usr/bin/env python

# Load required modules
import multiprocessing as mp
from time import time

# Load local modules
from exclusivity_tests import weighted_test, unweighted_test
from constants import *
from benjamini_hochberg import benjamini_hochberg_correction

# Construct a contingency table for an arbitrarily-sized gene set
def observed_values( M, N, geneToCases ):
    # Construct the table
    k = len(M)
    mutated_patients = set( p for g in M for p in geneToCases[g] )
    tbl = [0] * 2**k
    T, Z = 0, 0
    for p in mutated_patients:
        # Determine the binary index representing which genes
        # have mutations in this sample
        i = 0
        genes_mutated = 0
        for j, g in enumerate(M):
            if p in geneToCases[g]:
                genes_mutated += 1
                i = i | (1 << j)
        tbl[i] += 1

        # Increment the Z/T count if this mutation is exclusive or
        # co-occurring
        Z += int(genes_mutated > 1)
        T += int(genes_mutated == 1)

    # Finish the table and compute X
    tbl[0] = N - len(mutated_patients)
    X = [ len(geneToCases[g]) for g in M ]
    
    return X, T, Z, tbl

# Test the given sets with the given method and test
def test_set_group_wrapper(args): return test_set_group(*args)
def test_set_group( sets, geneToCases, num_patients, method, test, P=None, verbose=0 ):
    # Construct the arguments to test each set
    setToPval, setToTime, setToObs = dict(), dict(), dict()
    num_sets = len(sets)
    for i, M in enumerate(sets):
        # Simple progress bar
        if verbose > 1:
            sys.stdout.write('\r* Testing {}/{} triples...'.format(i+1, num_sets))
            sys.stdout.flush()

        # Do some simple mutation processing
        X, T, Z, tbl = observed_values(sorted(M), num_patients, geneToCases )
        
        # Ignore the opposite tail, where we have more co-occurrences than exclusivity
        if Z >= T or tbl[1] == 0 or tbl[2] == 0 or tbl[4] == 0 : continue

        # Compute the saddlepoint approximations
        start = time()
        if test == WEIGHTED:
            setToPval[M] = weighted_test( T, X, [ P[g] for g in M ], method )
        elif test == UNWEIGHTED:
            setToPval[M] = unweighted_test( T, X, tbl, method )
        else:
            raise NotImplementedError("Test {} not implemented".format(testToName[test]))
        
        setToTime[M] = time() - start
        setToObs[M]  = (T, Z, tbl)

    if verbose > 1: print

    return setToPval, setToTime, setToObs

def test_sets( sets, geneToCases, num_patients, method, test, P=None, num_cores=1, verbose=0):
    # Set up the multiprocessing
    num_cores = num_cores if num_cores != -1 else mp.cpu_count()
    if num_cores != 1:
        pool = mp.Pool(num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    # Split up the sets and run multiprocessing
    args = [ (sets[i::num_cores], geneToCases, num_patients, method, test, P, verbose)
             for i in range(num_cores) ]
    results = map_fn(test_set_group_wrapper, args)

    if num_cores != 1:
        pool.close()
        pool.join()

    # Combine the dictionaries
    setToPval, setToTime, setToObs = dict(), dict(), dict()
    for pval, time, obs in results:
        setToPval.update(pval.items())
        setToTime.update(time.items())
        setToObs.update(obs.items())

    # Compute the FDRs
    tested_sets = setToPval.keys()
    pvals = [ setToPval[M] for M in tested_sets ]
    setToFDR = dict(zip(tested_sets, benjamini_hochberg_correction(pvals)))
        
    return setToPval, setToTime, setToFDR, setToObs
