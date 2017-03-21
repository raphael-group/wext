#!/usr/bin/env python

# Load required modules
import sys, os, argparse, numpy as np, random
from collections import defaultdict

# Get argument parser
def get_parser():
    # General arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_prefix', type=str, required=True)
    parser.add_argument('-rs', '--random_seed', type=int, required=True)

    # Mode-specific arguments
    subparser = parser.add_subparsers(dest='mode', help='Data generation mode')
    pancan_parser = subparser.add_parser("pancan")
    pancan_parser.add_argument('-ns', '--num_samples', dest='N', type=int, required=True, help='We need at least 25.')
    pancan_parser.add_argument('-ng', '--num_genes', dest='M', type=int, required=True, help='We need at least four.')
    pancan_parser.add_argument('-b', '--bmr', type=float, required=True)
    return parser

# Generate Pan-Cancer data:
# - 2 datasets/cancer types
# - (G1, G2) with "across" exclusivity
# - (G3, G4) with "between" exclusivity
def generate_pancan_data(M, N, bmr, output_prefix):
    # Generate samples
    assert(N >= 25)
    N_1 = N/2
    N_2 = N-N_1
    samples1 = [ 'd1-%s' % (i+1) for i in range(N_1)]
    samples2 = [ 'd2-%s' % (i+1) for i in range(N_2)]
    samples  = samples1 + samples2

    # Generate genes (we need at least 4)
    assert( M >= 4 )
    genes = [ 'G%s' % (i+1) for i in range(M) ]

    # Generate mutation data, first implanting our known patterns
    sampleToMutations = defaultdict(set)

    across_samples = set(random.sample(samples, N/4))
    across_cases1  = set(random.sample(across_samples, N/8))
    for s in across_cases1:
        sampleToMutations[s].add( genes[0] )
    for s in across_samples - across_cases1:
        sampleToMutations[s].add( genes[1] )

    for s in set(random.sample(samples1, N/8)):
        sampleToMutations[s].add( genes[2] )
    for s in set(random.sample(samples2, N/8)):
        sampleToMutations[s].add( genes[3] )

    # Then adding noise (add default one mutation per sample)
    for s in samples: sampleToMutations[s].add( random.choice(genes) )
    for g in genes:
        for s in samples1:
            if np.random.rand() < bmr:
                sampleToMutations[s].add( g )
        for s in samples2:
            if np.random.rand() < 2*bmr:
                sampleToMutations[s].add( g )

    # Output to file each of the datasets
    with open('%s1-aberrations.tsv' % output_prefix, 'w') as OUT:
        OUT.write('\n'.join('%s\t%s' % (s, '\t'.join(sorted(sampleToMutations[s]))) for s in samples1))

    with open('%s2-aberrations.tsv' % output_prefix, 'w') as OUT:
        OUT.write('\n'.join('%s\t%s' % (s, '\t'.join(sorted(sampleToMutations[s]))) for s in samples2))


def run(args):
    # Seed the random generators
    np.random.seed(args.random_seed)
    random.seed(args.random_seed)

    # Run the given data generation mode
    if args.mode == 'pancan':
        generate_pancan_data(args.M, args.N, args.bmr, args.output_prefix)
    else:
        raise NotImplementedError('Data generation mode "%s" is not implemented.' % args.mode)
    return

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )