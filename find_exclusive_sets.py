#!/usr/bin/env python

# Load required modules
import sys, os, argparse, numpy as np, json
from itertools import combinations
from collections import defaultdict
from time import time

# Load the weighted exclusivity test, ensuring that
# it is in the path (unless this script was moved)
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from weighted_exclusivity_test import *

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_files', type=str, required=True, nargs='*')
    parser.add_argument('-paf', '--patient_annotation_file', type=str, required=False, default=None)
    parser.add_argument('-o', '--output_prefix', type=str, required=True)
    parser.add_argument('-f', '--min_frequency', type=int, default=1, required=False)
    parser.add_argument('-c', '--num_cores', type=int, required=False, default=1)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1, choices=range(5))
    parser.add_argument('--json_format', action='store_true', default=False, required=False)

    # Search strategy
    parser.add_argument('-s', '--search_strategy', type=str, choices=['Enumerate', 'MCMC'], default='MCMC', required=False)
    parser.add_argument('-ks', '--gene_set_sizes', nargs="*", type=int, required=True)
    parser.add_argument('-N', '--num_iterations', type=int, default=pow(10, 3))
    parser.add_argument('-nc', '--num_chains', type=int, default=1)
    parser.add_argument('-sl', '--step_length', type=int, default=100)
    parser.add_argument('-a', '--alpha', type=float, default=2., required=False, help='CoMEt parameter alpha.')
    parser.add_argument('--mcmc_seed', type=int, default=int(time()), required=False)

    # Subparser: statistical test
    subparser1 = parser.add_subparsers(dest='test', help='Type of test')

    permutational_parser = subparser1.add_parser("Permutational")
    permutational_parser.add_argument('-np', '--num_permutations', type=int, required=True)
    permutational_parser.add_argument('-pf', '--permuted_matrix_files', type=str, nargs='*', required=True)

    weighted_parser = subparser1.add_parser("Weighted")
    weighted_parser.add_argument('-m', '--method', choices=METHOD_NAMES, type=str, required=True)
    weighted_parser.add_argument('-wf', '--weights_files', type=str, required=True, nargs='*')

    unweighted_parser = subparser1.add_parser("Unweighted")
    unweighted_parser.add_argument('-m', '--method', choices=METHOD_NAMES, type=str, required=True)

    return parser

def get_permuted_files(permuted_matrix_directories, num_permutations):
    # Group and restrict the list of files we're testing
    permuted_directory_files = []
    for permuted_matrix_dir in permuted_matrix_directories:
        files = os.listdir(permuted_matrix_dir)
        permuted_matrices = [ '{}/{}'.format(permuted_matrix_dir, f) for f in files ]
        permuted_directory_files.append( permuted_matrices[:num_permutations] )
    assert( len(files) == num_permutations for files in permuted_directory_files )

    return zip(*permuted_directory_files)

# Load a list of weights files, merging them at the patient and gene level.
# Note that if a (gene, patient) pair is present in more than one file, it will
# be overwritten.
def load_weight_files(weights_files, genes, patients, typeToGeneIndex, typeToPatientIndex, masterGeneToIndex, masterPatientToIndex):
    # Master matrix of all weights
    P = np.zeros((len(genes), len(patients)))
    for i, weights_file in enumerate(weights_files):
        # Load the weights matrix for this cancer type and update the entries appropriately.
        # Note that since genes/patients can be measured in multiple types, we need to map
        # each patient to the "master" index.
        type_P                 = np.load(weights_file)

        ty_genes               = set(typeToGeneIndex[i].keys()) & genes
        ty_gene_indices        = [ typeToGeneIndex[i][g] for g in ty_genes ]
        master_gene_indices    = [ masterGeneToIndex[g] for g in ty_genes ]
        
        ty_patients            = set(typeToPatientIndex[i].keys()) & patients
        ty_patient_indices     = [ typeToPatientIndex[i][p] for p in ty_patients ]
        master_patient_indices = [ masterPatientToIndex[p] for p in ty_patients ]

        master_mesh            = np.ix_(master_gene_indices, master_patient_indices)
        ty_mesh                = np.ix_(ty_gene_indices, ty_patient_indices)
        P[ master_mesh ]       = type_P[ ty_mesh  ]

    # Set any zero entries to the minimum (pseudocount). The only reason for zeros is if
    #  a gene wasn't mutated at all in a particular dataset.
    P[P == 0] = np.min(P[P > 0])
    
    return dict( (g, P[masterGeneToIndex[g]]) for g in genes )

def load_mutation_files(mutation_files):
    typeToGeneIndex, typeToPatientIndex = [], []
    genes, patients, geneToCases = set(), set(), defaultdict(set)
    for i, mutation_file in enumerate(mutation_files):
        # Load ALL the data, we restrict by mutation frequency later
        mutation_data = load_mutation_data( mutation_file, 0 )
        _, type_genes, type_patients, typeGeneToCases, _, params, _ = mutation_data

        # We take the union of all patients and genes
        patients |= set(type_patients)
        genes    |= set(type_genes)

        # Record the mutations in each gene
        for g, cases in typeGeneToCases.iteritems(): geneToCases[g] |= cases

        # Record the genes, patients, and their indices for later
        typeToGeneIndex.append(dict(zip(type_genes, range(len(type_genes)))))
        typeToPatientIndex.append(dict(zip(type_patients, range(len(type_patients)))))

    return genes, patients, geneToCases, typeToGeneIndex, typeToPatientIndex

def run( args ):
    # Provide additional checks on arguments
    if args.test == 'Permutational':
        assert( len(args.mutation_files) == len(args.permuted_matrix_directories ) )
        assert( args.strategy != "MCMC" ) # MCMC is not implemented for permutational
    elif args.test == 'Weighted':
        assert( len(args.mutation_files) == len(args.weights_files) )

    # Load the mutation data
    if args.verbose > 0:
        print ('-' * 30), 'Input Mutation Data', ('-' * 29)
    genes, patients, geneToCases, typeToGeneIndex, typeToPatientIndex = load_mutation_files( args.mutation_files )
    num_all_genes, num_patients = len(genes), len(patients)

    # Restrict to genes mutated in a minimum number of samples
    geneToCases = dict( (g, cases) for g, cases in geneToCases.iteritems() if g in genes and len(cases) >= args.min_frequency )
    genes     = set(geneToCases.keys())
    num_genes = len(genes)
    
    # Load patient annotations (if provided) and add per patient events
    if args.patient_annotation_file:
        annotationToPatients = load_patient_annotation_file(args.patient_annotation_file)
        annotations = set( annotationToPatients.keys() )
        genes |= annotations

        # Since we are looking for co-occurrence between exclusive sets with
        # an annotation A, we add events for each patient NOT annotated by
        # the given annotation
        for annotation, cases in annotationToPatients.iteritems():
            not_cases = patients - cases
            if len(not_cases) > 0:
                geneToCases[annotation] = not_cases
    else:
        annotations = set()

    if args.verbose > 0:
        print '- Genes:', num_all_genes
        print '- Patients:', num_patients
        print '- Genes mutated in >={} patients: {}'.format(args.min_frequency, num_genes)
        if args.patient_annotation_file:
            print '- Patient annotations:', len(annotations)

    # Load the weights (if necessary)
    test = nameToTest[args.test]
    if test == WEIGHTED:
        # Create master versions of the indices
        masterGeneToIndex    = dict(zip(sorted(genes), range(num_genes)))
        masterPatientToIndex = dict( zip(sorted(patients), range(num_patients)) )
        geneToP = load_weight_files(args.weights_files, genes, patients, typeToGeneIndex, typeToPatientIndex, masterGeneToIndex, masterPatientToIndex)
    else:
        geneToP = None

    # Find the permuted matrices (if necessary)
    if test == PERMUTATIONAL:
        permuted_files = get_permuted_files(args.permuted_matrix_directories, args.num_permutations)
        if args.verbose > 0:
            print '* Using {} permuted matrix files'.format(len(permuted_files))

    # Search for exclusive sets
    method = nameToMethod[args.method]

    #Enumeration
    if args.search_strategy == 'Enumerate':
        if args.verbose > 0: print ('-' * 31), 'Enumerating Sets', ('-' * 31)
        for k in set( args.gene_set_sizes ): # we don't need to enumerate the same size more than once
            # Create a list of sets to test
            sets = list( frozenset(t) for t in combinations(genes, k) )
            num_sets = len(sets)

            if args.verbose  > 0: print 'k={}: {} sets...'.format(k, num_sets)
            if test == PERMUTATIONAL:
                # Run the permutational
                setToPval, setToRuntime, setToFDR, setToObs = permutational_test( sets, geneToCases, num_patients, permuted_files, args.num_cores, args.verbose )
            else:
                # Run the test
                method = nameToMethod[args.method]
                setToPval, setToRuntime, setToFDR, setToObs = test_sets(sets, geneToCases, num_patients, method, test, geneToP, args.num_cores, args.verbose)
            output_enumeration_table( args, k, setToPval, setToRuntime, setToFDR, setToObs )

    # MCMC
    elif args.search_strategy == 'MCMC':
        mcmc_params = dict(annotations=annotations, niters=args.num_iterations, nchains=args.num_chains, step_len=args.step_length, verbose=args.verbose, seed=args.mcmc_seed)
        setsToFreq, setToPval, setToObs = mcmc(args.gene_set_sizes, geneToCases, num_patients, method, test, geneToP, **mcmc_params)
        output_mcmc(args, setsToFreq, setToPval, setToObs)
    else:
        raise NotImplementedError("Strategy '{}' not implemented.".format(args.strategy))

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
