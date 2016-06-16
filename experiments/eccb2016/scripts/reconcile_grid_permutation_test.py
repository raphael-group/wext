#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, multiprocessing as mp
from collections import defaultdict

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_directory', type=str, required=True)
parser.add_argument('-b', '--batch_size', type=int, required=True)
parser.add_argument('-c', '--num_cores', type=int, required=False, default=1)
parser.add_argument('-wd', '--wext_directory', type=str, required=True)
parser.add_argument('-o', '--output_prefix', type=str, required=True)
args = parser.parse_args( sys.argv[1:] )

# Load WExT
sys.path.append( args.wext_directory )
from wext import *

# Load and merge the JSON files
def load_json_files(( json_files )):
    setToCount       = defaultdict( int )
    setToRuntime     = defaultdict( float )
    setToObs         = dict()
    num_permutations = 0.
    for i, json_file in enumerate(json_files):
        # Parse the JSON file
        with open(json_file, 'r') as IN:
            obj = json.load(IN)
            for M, pval in obj['setToPval'].iteritems():
                frozen_M = frozenset(M.split('\t'))
                setToCount[frozen_M]   += int(round(pval * args.batch_size))
                setToRuntime[frozen_M] += obj['setToRuntime'][M]
                setToObs[frozen_M]      = obj['setToObs'][M]
                    
        num_permutations += args.batch_size

    return setToCount, setToRuntime, setToObs, num_permutations

# Extract the JSON files
for root, dirs, files in os.walk( args.input_directory ):
    json_files = [ '{}/{}'.format(root, f) for f in files if f.lower().endswith('.json') ]

# Set up the multiprocessing and run
print '* Loading {} JSON files...'.format(len(json_files))
num_cores = args.num_cores if args.num_cores != -1 else mp.cpu_count()
if num_cores != 1:
    pool = mp.Pool(num_cores)
    map_fn = pool.imap
else:
    map_fn = map

mp_args = [ (json_files[i::num_cores]) for i in range(num_cores) ]
results = map_fn(load_json_files, mp_args)

if num_cores != 1:
    pool.close()
    pool.join()

# Merge the results
print '\t- Merging results...'
setToCount       = defaultdict( int )
setToRuntime     = defaultdict( float )
setToObs         = dict()
num_permutations = 0.
for counts, runtimes, obs, N in results:
    num_permutations += N
    setToObs.update( obs.items() )
    for M, count in counts.iteritems():
        setToCount[M]   += count
        setToRuntime[M] += runtimes[M]
    
setToPval = dict( (M, count/num_permutations) for M, count in setToCount.iteritems() )

print '\t- Loaded {} sets with {} permutations'.format(len(setToPval), int(num_permutations))

# Compute FDR
print '* Computing FDRs...'
tested_sets = setToPval.keys()
pvals       = [ setToPval[M] for M in tested_sets ]
setToFDR    = dict(zip(tested_sets, benjamini_hochberg_correction(pvals, independent=False)))

# Output the merged file
print '* Outputting to file...'
k                  = len(tested_sets[0])
args.json_format   = True
args.test          = 'RCE'
output_enumeration_table( args, k, setToPval, setToRuntime, setToFDR, setToObs )

