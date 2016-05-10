#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json
from collections import defaultdict
from math import isnan
from rank import rank

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-uf', '--unweighted_exact_file', type=str, required=True)
parser.add_argument('-wf', '--weighted_saddlepoint_file', type=str, required=True)
parser.add_argument('-mf', '--mutation_file', type=str, required=True)
parser.add_argument('-lf', '--lengths_file', type=str, required=True)
parser.add_argument('-o', '--output_prefix', type=str, required=True)
parser.add_argument('-nt', '--num_triples', type=int, required=True)
parser.add_argument('-lt', '--length_threshold', type=int, required=True)
parser.add_argument('-f', '--fdr_cutoff', type=float, required=True)
args = parser.parse_args(sys.argv[1:])

# Load the lengths
with open(args.lengths_file, 'r') as IN:
    arrs = [ l.rstrip().split('\t') for l in IN if not l.startswith('#') ]
    geneToLength = dict( (arr[0], float(arr[1])) for arr in arrs )
    lengths = geneToLength.values()
    length_ranks = rank(lengths, reverse=True)
    geneToLengthRank = defaultdict( lambda : args.length_threshold + 1 )
    geneToLengthRank.update(zip(geneToLength.keys(), length_ranks))
    threshold_gene = sorted(geneToLength.keys(), key=lambda g: geneToLengthRank[g])[args.length_threshold]
    print 'Length of {} longest gene: {}'.format(args.length_threshold, geneToLength[threshold_gene])
    
# Load the mutations
with open(args.mutation_file, 'r') as IN:
    obj = json.load(IN)
    genes, patients = obj['genes'], obj['patients']
    hypermutators = set(obj['hypermutators'])
    geneToCases = dict( (g, set(cases)) for g, cases in obj['geneToCases'].iteritems() )

# Load the triples
with open(args.unweighted_exact_file, 'r') as IN:
    obj            = json.load(IN)
    unweightedPval = dict( (frozenset(t.split('\t')), pval) for t, pval in obj['setToPval'].iteritems()  )
    assert( all( not(isnan(pval)) for pval in unweightedPval.values() ))
    unweightedFDR  = dict( (frozenset(t.split('\t')), fdr) for t, fdr in obj['setToFDR'].iteritems() )

with open(args.weighted_saddlepoint_file, 'r') as IN:
    obj          = json.load(IN)
    weightedPval = dict( (frozenset(t.split('\t')), pval) for t, pval in obj['setToPval'].iteritems() )
    assert( all( not(isnan(pval)) for pval in weightedPval.values() ))
    weightedFDR  = dict( (frozenset(t.split('\t')), fdr) for t, fdr in obj['setToFDR'].iteritems() )

print 'Triples with weighted FDR < {}: {}/{}'.format(args.fdr_cutoff, sum(1 for t, fdr in weightedFDR.iteritems() if fdr < args.fdr_cutoff), len(weightedFDR))
print 'Triples with unweighted FDR < {}: {}/{}'.format(args.fdr_cutoff, sum(1 for t, fdr in unweightedFDR.iteritems() if fdr < args.fdr_cutoff), len(unweightedFDR))
    
# Rank triples by P-value
triples = sorted(set(weightedPval.keys()) & set(unweightedPval.keys()))
top_weighted_triples = sorted(triples, key=lambda t: weightedPval[t])
weightedRank = dict(zip(triples, rank([ weightedPval[t] for t in triples ])))
top_unweighted_triples = sorted(triples, key=lambda t: unweightedPval[t])
unweightedRank = dict(zip(triples, rank([ unweightedPval[t] for t in triples ])))

# Create tables
def length_indicate(g):
    if geneToLengthRank[g] > args.length_threshold: return g
    else: return '\\textbf{%s}' % g

header = ['CoMEt rank', 'Weighted rank', 'Triple', 'Phi(M)', 'Psi(M)', 'Hypermutator mutations']
tbl = [ header ]
latex_tbl = [ header ]
for t in top_unweighted_triples[:args.num_triples] + top_weighted_triples[:args.num_triples]:
    cases = set( p for g in t for p in geneToCases[g] )
    tbl.append([ unweightedRank[t]+1, weightedRank[t]+1,
                 ', '.join(sorted(t)), '{:.2E}'.format(unweightedPval[t]),
                 '{:.2E}'.format(weightedPval[t]), len(cases & hypermutators) ])
    latex_tbl.append([ unweightedRank[t]+1, weightedRank[t]+1,
                       ', '.join([ length_indicate(g) for g in sorted(t)]),
                       '{:.2E}'.format(unweightedPval[t]),
                       '{:.2E}'.format(weightedPval[t]),
                       len(cases & hypermutators) ])

# Output tables to file
with open(args.output_prefix + '.tex', 'w') as OUT:
    OUT.write('\\\\\n'.join([ ' & '.join(map(str, row)) for row in latex_tbl ]) )
with open(args.output_prefix + '.tsv', 'w') as OUT:
    OUT.write('\n'.join([ '\t'.join(map(str, row)) for row in tbl ]) )
