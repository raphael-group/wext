#!/usr/bin/env python

# Load required modules
import sys, os, argparse
from itertools import combinations

# Parse arguments
job_id = os.environ['SGE_TASK_ID'] if 'SGE_TASK_ID' in os.environ else None
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mutation_file', type=str, required=True)
parser.add_argument('-f', '--min_freq', type=int, required=True)
parser.add_argument('-k', '--gene_set_size', type=int, required=True)
parser.add_argument('-b', '--batch_size', type=int, required=True)
parser.add_argument('-i', '--input_directory', type=str, required=True)
parser.add_argument('-np', '--num_permutations', type=int, required=True)
parser.add_argument('-o', '--output_prefix', type=str, required=True)
parser.add_argument('-w', '--wext_directory', type=str, required=True)
parser.add_argument('-j', '--job_id', type=int, required=job_id is None, default=job_id)
parser.add_argument('-v', '--verbose', type=int, required=False, default=0, choices=range(5))
args = parser.parse_args( sys.argv[1:] )

# Load weighted exclusivity test
sys.path.append(args.wext_directory)
from find_exclusive_sets import get_permuted_files
from wext import rce_permutation_test, load_mutation_data, output_enumeration_table

# Load the mutation data
if args.verbose > 0: print '* Loading mutation data..'
mutation_data = load_mutation_data( args.mutation_file, args.min_freq )
genes, all_genes, patients, geneToCases, _, params, _ = mutation_data
num_patients = len(patients)
sets = list( frozenset(t) for t in combinations(genes, args.gene_set_size) )

if args.verbose > 0: print '\t- Testing {} sets of size k={}'.format(len(sets), args.gene_set_size)

# Run the permutational test
if args.verbose > 0: print '* Running permutation test...'
start_index = (args.job_id-1) * args.batch_size
permuted_files = get_permuted_files([args.input_directory], args.num_permutations)[start_index:start_index + args.batch_size]
if args.verbose > 0: print '\t- Testing {} files'.format(len(permuted_files))
    
setToPval, setToRuntime, setToFDR, setToObs = rce_permutation_test( sets, geneToCases, num_patients, permuted_files, 1, 0 )

# Output to file
args.output_prefix += '-%s' % job_id
args.test = 'RCE'
args.json_format = True
output_enumeration_table( args, args.gene_set_size, setToPval, setToRuntime, setToFDR, setToObs )
