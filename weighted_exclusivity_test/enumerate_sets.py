#!/usr/bin/env python

# Load required modules
import multiprocessing as mp, json
from time import time
from collections import defaultdict, Counter
from math import ceil, isnan

# Load local modules
from exclusivity_tests import weighted_test, unweighted_test
from constants import *
from benjamini_hochberg import benjamini_hochberg_correction

################################################################################
# Permutational test
################################################################################
# Compute the mutual exclusivity T for the given gene set
def T(M, geneToCases):
    sampleToCount = Counter( s for g in M for s in geneToCases.get(g, []) )
    return sum( 1 for sample, count in sampleToCount.iteritems() if count == 1 )

# Compute the permutational 
def permutational_dist_wrapper( args ): return permutational_dist( *args )
def permutational_dist( sets, permuted_files ):
    setToDist, setToTime = defaultdict(list), defaultdict(int)
    for pf_group in permuted_files:
        # Load the file, keeping track of how long it takes
        reading_start = time()
        permutedGeneToCases = defaultdict(set)
        for pf in pf_group:
            with open(pf, 'r') as IN:
                for g, cases in json.load(IN)['geneToCases'].iteritems():
                    permutedGeneToCases[g] |= set(cases)
                    
        reading_time = time() - reading_start

        # Iterate through the sets, keeping track of how long
        # it takes to compute the test statistic
        for M in sets:
            start = time()
            setToDist[M].append( T(M, permutedGeneToCases) )
            setToTime[M] += reading_time + (time() - start)
            
    return setToDist, setToTime

def permutational_test(sets, geneToCases, num_patients, permuted_files, num_cores=1, verbose=0):
    # Set up the multi-core process
    num_cores = num_cores if num_cores != -1 else mp.cpu_count()
    if num_cores != 1:
        pool = mp.Pool(num_cores)
        map_fn = pool.map
    else:
        map_fn = map
        
    # Filter the sets based on the observed values
    k = len(next(iter(sets)))
    setToObs = dict( (M, observed_values(M, num_patients, geneToCases)) for M in sets )
    sets = set( M for M, (X, T, Z, tbl) in setToObs.iteritems() if testable_set(k, T, Z, tbl) )
    
    # Compute the distribution of exclusivity for each pair across the permuted files
    np    = float(len(permuted_files))
    args  = [ (sets, permuted_files[i::num_cores]) for i in range(num_cores) ]
    empirical_distributions = map_fn(permutational_dist_wrapper, args)

    if num_cores != 1:
        pool.close()
        pool.join()

    # Merge the different distributions
    setToDist, setToTime = defaultdict(list), dict()
    for dist, times in empirical_distributions:
        setToTime.update(times.items())
        for k, v in dist.iteritems():
            setToDist[k].extend(v)

    # Compute the observed values and then the P-values
    setToObs = dict( (M, setToObs[M]) for M in sets )
    setToPval = dict()
    for M, (X, T, Z, tbl) in setToObs.iteritems():
        # Compute the P-value.
        count = sum( 1. for d in setToDist[M] if d >= T )
        setToPval[M] = count / np

    # Compute FDRs
    tested_sets = setToPval.keys()
    pvals = [ setToPval[M] for M in tested_sets ]
    setToFDR = dict(zip(tested_sets, benjamini_hochberg_correction(pvals, independent=False)))
        
    return setToPval, setToTime, setToFDR, setToObs

################################################################################
# Weighted and unweighted tests
################################################################################
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
    k = len(next(iter(sets)))
    for i, M in enumerate(sets):
        # Simple progress bar
        if verbose > 1:
            sys.stdout.write('\r* Testing {}/{} triples...'.format(i+1, num_sets))
            sys.stdout.flush()

        # Do some simple mutation processing
        sorted_M = sorted(M)
        X, T, Z, tbl = setToObs[M] = observed_values(sorted_M, num_patients, geneToCases )
        
        # Ignore the opposite tail, where we have more co-occurrences than exclusivity
        if not testable_set(k, T, Z, tbl): continue

        # Compute the saddlepoint approximations
        start = time()
        if test == WEIGHTED:
            setToPval[M] = weighted_test( T, X, [ P[g] for g in sorted_M ], method )
        elif test == UNWEIGHTED:
            setToPval[M] = unweighted_test( T, X, tbl, method )
        else:
            raise NotImplementedError("Test {} not implemented".format(testToName[test]))
        
        setToTime[M] = time() - start

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

    # Make sure all P-values are numbers
    tested_sets = len(setToPval)
    invalid_sets = set( M for M, pval in setToPval.iteritems() if isnan(pval) or -PTOL > pval or pval > 1+PTOL )
    setToPval = dict( (M, pval) for M, pval in setToPval.iteritems() if not M in invalid_sets )
    setToTime = dict( (M, runtime) for M, runtime in setToTime.iteritems() if not M in invalid_sets )
    setToObs = dict( (M, obs) for M, obs in setToObs.iteritems() if not M in invalid_sets )
        
    if verbose > 0:
        print '* Output sets: {}'.format(len(setToPval))
        print '\t-Removed {} sets with NaN or invalid P-values'.format(len(invalid_sets))
        print '\t-Ignored {} sets with Z >= T or a gene with no exclusive mutations'.format(len(sets)-tested_sets)
        
    # Compute the FDRs
    tested_sets = setToPval.keys()
    pvals = [ setToPval[M] for M in tested_sets ]
    setToFDR = dict(zip(tested_sets, benjamini_hochberg_correction(pvals, independent=False)))
        
    return setToPval, setToTime, setToFDR, setToObs

################################################################################
# Helpers
################################################################################
# Testable set
def testable_set( k, T, Z, tbl ):
    return T > Z and all( tbl[2**i] > 0 for i in range(k) )
