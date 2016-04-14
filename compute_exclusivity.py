#!/usr/bin/env python

# Load required modules
import sys, os, argparse, numpy as np, json
from itertools import combinations

# Load the weighted enrichment test, ensuring that
# it is in the path (unless this script was moved)
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from weighted_exclusivity_test import *

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-f', '--min_frequency', type=int, default=1, required=False)
    parser.add_argument('-k', '--gene_set_size', type=int, required=True,
        choices=SET_SIZES_IMPLEMENTED)
    parser.add_argument('-tn', '--test', choices=TEST_NAMES, required=True)
    parser.add_argument('-mn', '--method', choices=METHOD_NAMES, required=True)
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1, choices=range(5))

    # Subparsers
    parser.add_argument('-wf', '--weights_file', type=str, required=True)
    parser.add_argument('-np', '--num_permutations', type=int, required=False)
    return parser

def run( args ):
    # Load the mutation data
    if args.verbose > 0:
        print '* Loading mutation data...'
        
    mutation_data = load_mutation_data( args.mutation_file, args.min_frequency )
    genes, all_genes, patients, geneToCases, params, hypermutators = mutation_data
    num_all_genes, num_genes, num_patients = len(all_genes), len(genes), len(patients)
    geneToIndex = dict(zip(all_genes, range(num_all_genes)))
    
    if args.verbose > 0:
        print '\tGenes:', num_all_genes
        print '\tPatients:', num_patients
        print '\tGenes mutated in >={} patients: {}'.format(args.min_frequency, num_genes)

    # Load the weights file
    P       = np.load(args.weights_file)
    geneToP = dict( (g, P[geneToIndex[g]]) for g in genes )
    
    # Create a list of sets to test
    sets = list( frozenset(t) for t in combinations(genes, args.gene_set_size) )
    num_sets = len(sets)
    
    # Test each set
    print '* Testing {} sets...'.format(num_sets)
    method = nameToMethod[args.method]
    test   = nameToTest[args.test]
    setToPval, setToRuntime, setToFDR, setToObs = test_sets(sets, geneToCases, num_patients,
                                                            method, test, geneToP, args.verbose,
                                                            args.num_cores)

    # Output to file
    output_table( args, setToPval, setToRuntime, setToFDR, setToObs )

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
