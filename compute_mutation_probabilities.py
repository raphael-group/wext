#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, numpy as np, multiprocessing as mp
from collections import defaultdict

# Load the weighted exclusivity test
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(this_dir)
from wext import *

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_file', type=str, required=True)
    parser.add_argument('-wf', '--weights_file', type=str, required=False, default=None)
    parser.add_argument('-pd', '--permutation_directory', type=str, required=False)
    parser.add_argument('-np', '--num_permutations', type=int, required=True)
    parser.add_argument('-si', '--start_index', type=int, required=False, default=1)
    parser.add_argument('-q', '--swap_multiplier', type=int, required=False, default=100)
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1, choices=range(5))
    return parser

def permute_matrices_wrapper(args): return permute_matrices(*args)
def permute_matrices(edge_list, max_swaps, max_tries, seeds, verbose,
                     m, n, num_edges, indexToGene, indexToPatient):
    # Initialize our output
    observed     = np.zeros((m, n))
    permutations = []
    for seed in seeds:
        # Permute the edge list
        permuted_edge_list = bipartite_edge_swap(edge_list, max_swaps, max_tries, seed, verbose,
                                                 m, n, num_edges)

        # Recover the mapping of mutations from the permuted edge list
        geneToCases  = defaultdict(set)
        indices = []
        for edge in permuted_edge_list:
            gene, patient = indexToGene[edge[0]], indexToPatient[edge[1]]
            geneToCases[gene].add(patient)
            indices.append( (edge[0]-1, edge[1]-1) )

        # Record the permutation
        observed[zip(*indices)] += 1.
        geneToCases = dict( (g, list(cases)) for g, cases in geneToCases.iteritems() )
        permutations.append( dict(geneToCases=geneToCases, permutation_number=seed) )

    return observed/float(len(seeds)), permutations

def run( args ):
    # Do some additional argument checking
    if not args.weights_file and not args.permutation_directory:
        sys.stderr.write('You must set the weights file or permutation directory, '\
                         'otherwise nothing will be output.')
        sys.exit(1)

    # Load mutation data
    if args.verbose > 0:
        print '* Loading mutation data...'

    mutation_data = load_mutation_data( args.mutation_file )
    genes, all_genes, patients, geneToCases, patientToMutations, params, hypermutators = mutation_data

    geneToObserved = dict( (g, len(cases)) for g, cases in geneToCases.iteritems() )
    patientToObserved = dict( (p, len(muts)) for p, muts in patientToMutations.iteritems() )
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
        print '* Permuting matrices...'

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

    wrapper_args = [ (edge_list, max_swaps, max_tries, seeds[i::num_cores], 0, m,
                      n, num_edges, indexToGene, indexToPatient) for i in range(num_cores) ]
    results = map_fn(permute_matrices_wrapper, wrapper_args)

    if num_cores != 1:
        pool.close()
        pool.join()

    # Create the weights file
    if args.weights_file:
        if args.verbose > 0:
            print '* Saving weights file...'

        # Merge the observeds
        observeds = [ observed for observed, _ in results ]
        P = np.add.reduce(observeds) / float(len(observeds))

        # Verify the weights
        for g, obs in geneToObserved.iteritems():
            assert( np.abs(P[geneToIndex[g]-1].sum() - obs) < 0.1)

        for p, obs in patientToObserved.iteritems():
            assert( np.abs(P[:, patientToIndex[p]-1].sum() - obs) < 0.1)

        # Add pseudocounts to entries with no mutations observed
        P[P == 0] = 1./(2. * args.num_permutations)

        # Output to file.
        # The rows/columns preserve the order given by the mutation file.
        np.save(args.weights_file, P)

    # Save the permuted mutation data
    if args.permutation_directory:
        output_prefix = args.permutation_directory + '/permuted-mutations-{}.json'
        if args.verbose > 0:
            print '* Saving permuted mutation data...'

        for _, permutation_list in results:
            for permutation in permutation_list:
                # Output in adjacency list format
                with open(output_prefix.format(permutation['permutation_number']), 'w') as OUT:
                    permutation['params'] = params
                    json.dump( permutation, OUT )

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
