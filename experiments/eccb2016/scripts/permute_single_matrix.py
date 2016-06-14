#!/usr/bin/env python

# Import modules.
import numpy as np, os, sys, argparse, json
from time import time
from collections import defaultdict

# Parse arguments.
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_file', type=str, required=True)
    parser.add_argument('-wd', '--wext_dir', type=str, required=True)
    parser.add_argument('-o', '--output_prefix', type=str, required=True)
    parser.add_argument('-q', '--swaps_multiplier', type=int, required=False, default=100)
    parser.add_argument('-s', '--seed', type=int, required=False, default=int(time()))
    parser.add_argument('-j', '--job_id', type=int, required=False,
                        default=os.environ.get('SGE_TASK_ID', 0))
    return parser

def run( args ):
    # Load WExT
    sys.path.append(args.wext_dir)
    from weighted_exclusivity_test import load_mutation_data, bipartite_edge_swap
    from compute_mutation_probabilities import permute_matrices
        
    # Load edge list.
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

    # Find sizes of sets and number of edges in graph.
    m, n = np.max(edge_list, axis=0)
    num_edges = len(edge_list)

    # Perform bipartite edge swap.
    max_swaps = int(round(args.swaps_multiplier*num_edges))
    permuted_edge_list = bipartite_edge_swap(edge_list, max_swaps=max_swaps,
                                             max_tries=2**30, seed=args.seed,
                                             verbose=0, m=m, n=n, num_edges=num_edges)

    # Convert the permuted edge list back to a dictionary
    permutedGeneToCases, permutedPatientToMutations = defaultdict( set ), defaultdict( set )
    for edge in permuted_edge_list:
        gene, patient = indexToGene[edge[0]], indexToPatient[edge[1]]
        permutedGeneToCases[gene].add(patient)
        permutedPatientToMutations[patient].add(gene)
        
    # Verify the number of mutations per gene/patient is preserved
    for g, cases in geneToCases.iteritems():
        assert( len(cases) == len(permutedGeneToCases[g]) )

    for p, muts in patientToMutations.iteritems():
        assert( len(muts) == len(permutedPatientToMutations[p]) )

    # Save edge list.
    output_file = '{}-{}.json'.format(args.output_prefix, args.job_id)
    permutation = dict(params=params, permutation_number=args.job_id,
                       geneToCases=dict( (g, list(cases)) for g, cases in permutedGeneToCases.iteritems()))
    with open(output_file, 'w') as OUT: json.dump( permutation, OUT )
    
if __name__ == '__main__':
    run( get_parser().parse_args(sys.argv[1:]) )
