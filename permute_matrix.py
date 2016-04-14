#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json
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
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1, choices=range(5))
    parser.add_argument('-r', '--birewirer_script', type=str, required=False,
        default='{}/weighted_exclusivity_test/src/R/birewire.R'.format(this_dir))
    args = parser.parse_args(sys.argv[1:])
    return parser

def run( args ):
    # Load the mutation matrix
    if args.verbose > 0:
        print '* Loading mutation data...'
        
    mutation_data = load_mutation_data( args.mutation_file )
    genes, all_genes, patients, geneToCases, patientToMutations, params, hypermutators = mutation_data
    num_all_genes, num_genes, num_patients = len(all_genes), len(genes), len(patients)
    geneToIndex = dict(zip(all_genes, range(num_all_genes)))
    
    if args.verbose > 0:
        print '\tGenes:', num_all_genes
        print '\tPatients:', num_patients
        print '\tGenes mutated in at least one patient: {}'.format(num_genes)

    # Output an edge list representing each mutation in the dataset
    edges = []
    for gene, cases in geneToCases.iteritems():
        for patient in cases:
            edges.append( '{}\t{}'.format(gene, patient) )

    original_edgelist_file = '{}/original-edgelist.txt'.format(args.output_directory)
    with open(original_edgelist_file, 'w') as OUT:
        OUT.write('\n'.join(edges))

    # Run the R script to create permuted edge lists
    print '* Running the permutation...'
    os.system('Rscript {} {} {} {} {}'.format(args.birewirer_script, original_edgelist_file, args.output_directory, args.num_permutations, args.start_index))

    # Convert the permuted edge lists
    print '* Converting edge lists...'
    patients = set(patients)
    for i in range(args.start_index, args.start_index + args.num_permutations):
        # Load the edgelists
        permuted_edgelist_file = '{}/permuted-edgelist-{}.txt'.format(args.output_directory, i)
        with open(permuted_edgelist_file, 'r') as IN:
            edges = [ l.rstrip().split() for l in IN ]
            
        # Figure out which column has the samples and create a mapping of samples to genes
        patientIndex = 0 if edges[0][0] in patients else 1
        geneIndex    = 0 if edges[0][0] in genes else 1
        geneToCases  = defaultdict(set)

        for e in edges:
            geneToCases[e[geneIndex]].add( e[patientIndex])
    
        # Output in adjacency list format
        with open('{}/permuted-mutations-{}.json'.format(args.output_directory, i), 'w') as OUT:
            output = dict(params=params, permutation_number=i, geneToCases=dict( (g, list(cases)) for g, cases in geneToCases.iteritems()))
            json.dump( output, OUT )

        # Clean up
        if os.path.isfile(permuted_edgelist_file): os.remove(permuted_edgelist_file)

    # Clean up
    if os.path.isfile(original_edgelist_file): os.remove(original_edgelist_file)


if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
