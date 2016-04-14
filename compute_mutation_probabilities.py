#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, numpy as np, multiprocessing as mp

def get_parser():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-pf', '--permuted_files', required=True, type=str, nargs='*')
    parser.add_argument('-mf', '--mutation_file', required=True, type=str)
    parser.add_argument('-np', '--num_permutations', required=True, type=int)
    parser.add_argument('-nc', '--num_cores', required=False, default=1, type=int)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    return parser

def compute_observed(mutation_files, geneToIndex, patientToIndex):
    # Compute the exclusivity for each mutation file
    observed = np.zeros((len(geneToIndex), len(patientToIndex)))

    for mutation_file in mutation_files:
        # Load the mutations
        with open(mutation_file, 'r') as IN:
            permutedGeneToCases = json.load(IN)['geneToCases']

        # Increase the count for all the observed mutations in this
        # permuted file
        for g, cases in permutedGeneToCases.iteritems():
            indices = [ (geneToIndex[g], patientToIndex[p]) for p in cases ]
            observed[zip(*indices)] += 1.

    return observed / float(len(mutation_files))

def compute_observed_wrapper((mutation_files, geneToIndex, patientToIndex)):
    return compute_observed(mutation_files, geneToIndex, patientToIndex)
    
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
        
    # Get the list of files we want to consider
    permuted_files = args.permuted_files[:args.num_permutations]
    print '* Analyzing {} permuted files...'.format(len(permuted_files))
                                                    
    # Set up the jobs
    num_cores = args.num_cores if args.num_cores != -1 else mp.cpu_count()
    if num_cores != 1:
        pool = mp.Pool(num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    compute_args = [ (permuted_files[i::args.num_cores], geneToIndex, patientToIndex)
                     for i in range(args.num_cores) ]

    observeds = map_fn(compute_observed_wrapper, compute_args)

    if num_cores != 1:
        pool.close()
        pool.join()

    # Merge the observeds
    P = np.add.reduce(observeds) / float(len(observeds))

    # Verify the weights
    for g, obs in geneToObserved.iteritems():
        assert( np.abs(P[geneToIndex[g]].sum() - obs) < 0.1)
        
    for p, obs in patientToObserved.iteritems():
        assert( np.abs(P[:, patientToIndex[p]].sum() - obs) < 0.1)
    
    # Output to file.
    # The rows/columns preserve the order given by the mutation file.
    np.save(args.output_file, P)

if __name__ == '__main__': run(get_parser().parse_args(sys.argv[1:]))
