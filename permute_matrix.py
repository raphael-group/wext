#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, numpy as np, multiprocessing as mp
from collections import defaultdict

# Load the weighted exclusivity test
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(this_dir)
from weighted_exclusivity_test import *

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_file', type=str, required=True)
    parser.add_argument('-o', '--output_directory', type=str, required=True)
    parser.add_argument('-np', '--num_permutations', type=int, required=True)
    parser.add_argument('-si', '--start_index', type=int, required=False, default=1)
    parser.add_argument('-q', '--swap_multiplier', type=int, required=False, default=100)
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1, choices=range(5))
    return parser

def bipartite_edge_swap_wrapper(args):
    return bipartite_edge_swap(*args)

def run( args ):
    # Load mutation data
    if args.verbose > 0:
        print '* Loading mutation data...'

    mutation_data = load_mutation_data( args.mutation_file )
    genes, all_genes, patients, geneToCases, patientToMutations, params, hypermutators = mutation_data

    geneToIndex = dict( (g, i+1) for i, g in enumerate(all_genes) )
    indexToGene = dict( (i+1, g) for i, g in enumerate(all_genes) )
    patientToIndex = dict( (p, j+1) for j, p in enumerate(patients) )
    indexToPatient = dict( (j+1, p) for j, p in enumerate(patients) )

    edges = set()
    for gene, cases in geneToCases.iteritems():
        for patient in cases:
            edges.add( (geneToIndex[gene], patientToIndex[patient]) )

    edge_list = np.array(sorted(edges), dtype=np.int)

    # Run the bipartite edge swaps
    if args.verbose > 0:
        print '* Running biparite edge swap...'

    m = len(all_genes)
    n = len(patients)
    num_edges = len(edges)
    max_swaps = int(args.swap_multiplier*num_edges)
    max_tries = 10**9
    seeds = [ i+args.start_index for i in range(args.num_permutations) ]

    # Run the bipartite edge swaps in parallel if more than one core indicated
    num_cores = args.num_cores if args.num_cores != -1 else mp.cpu_count()
    if num_cores != 1:
        pool = mp.Pool(num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    wrapper_args = [ (edge_list, max_swaps, max_tries, seed, 0, m, n, num_edges) for seed in seeds ]
    permuted_edge_lists = map_fn(bipartite_edge_swap_wrapper, wrapper_args)

    if num_cores != 1:
        pool.close()
        pool.join()

    # Save the permuted mutation data
    if args.verbose > 0:
        print '* Saving permuted mutation data...'
    patients = set(patients)
    for permuted_edge_list, seed in zip(permuted_edge_lists, seeds):
        geneToCases  = defaultdict(set)

        for edge in permuted_edge_list:
            gene, patient = indexToGene[edge[0]], indexToPatient[edge[1]]
            geneToCases[gene].add(patient)

        # Output in adjacency list format
        with open('{}/permuted-mutations-{}.json'.format(args.output_directory, seed), 'w') as OUT:
            output = dict(params=params, permutation_number=i, geneToCases=dict( (g, list(cases)) for g, cases in geneToCases.iteritems()))
            json.dump( output, OUT )

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
