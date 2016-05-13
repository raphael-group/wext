#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-mf', '--mutation_file', type=str, required=True)
parser.add_argument('-lf', '--lengths_file', type=str, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

# Load the lengths file
with open(args.lengths_file, 'r') as IN:
    arrs = [ l.rstrip().split('\t') for l in IN if not l.startswith('#') ]
    geneToLength = dict( (arr[0], float(arr[1])) for arr in arrs )

# Load the mutation file
with open(args.mutation_file, 'r') as IN:
    obj = json.load(IN)
    geneToCases = obj['geneToCases']

# Remove genes without a length
original_genes  = set(obj['genes'])
remaining_genes = original_genes & set(geneToLength.keys())
obj['geneToCases'] = dict( (g, cases) for g, cases in obj['geneToCases'].iteritems() if g in geneToLength )
obj['genes'] = sorted(obj['geneToCases'].keys())
obj['num_genes'] = len(obj['genes'])
obj['params']['lengths_file'] = os.path.abspath(args.lengths_file)
obj['genes_with_no_length_removed'] = sorted(original_genes - set(obj['genes']))
obj['patientToMutations'] = dict( (p, sorted(set(muts) & remaining_genes)) for p, muts in obj['patientToMutations'].iteritems() )
print 'Removed {} genes with no length'.format(len(obj['genes_with_no_length_removed']))

# Output the new file
with open(args.output_file, 'w') as OUT:
    json.dump( obj, OUT )
