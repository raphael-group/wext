#!/usr/bin/env python

import sys, os, json, numpy as np
from constants import *

# Load mutation data from one of our processed JSON file
def load_mutation_data( mutation_file, min_freq=1 ):
    with open(mutation_file, 'r') as IN:
        obj         = json.load(IN)
        all_genes   = obj['genes']
        patients    = obj['patients']
        geneToCases = dict( (g, set(cases)) for g, cases in obj['geneToCases'].iteritems() )
        patientToMutations = dict( (p, set(muts)) for p, muts in obj['patientToMutations'].iteritems() )
        hypermutators = set(obj['hypermutators'])
        params      = obj['params']

    # Restrict the genes based on the minimum frequency
    genes = set( g for g, cases in geneToCases.iteritems() if len(cases) >= min_freq )

    return genes, all_genes, patients, geneToCases, patientToMutations, params, hypermutators

# Converts keys from an iterable to tab-separated, so the dictionary can be
# output as JSON
def convert_dict_for_json( setToVal, sep='\t' ):
    return dict( (sep.join(sorted(M)), val) for M, val in setToVal.iteritems() )

# Converts tab-separated keys back to frozensets
def convert_dict_from_json( setToVal, sep='\t', iterable=frozenset ):
    return dict( (iterable(M.split(sep)), val) for M, val in setToVal.iteritems() )

# Create the header strings for a contingency table
def create_tbl_header( k ):
    bin_format = 't{0:0%sb}' % k
    return '\t'.join([ bin_format.format(i) for i in range(2**k) ])

# Output a run to file as a table or JSON file
def output_enumeration_table(args, k, setToPval, setToRuntime, setToFDR, setToObs ):
    is_permutational = nameToTest[args.test] == PERMUTATIONAL
    extension = 'json' if args.json_format else 'tsv'
    with open('{}-k{}.{}'.format(args.output_prefix, k, extension), 'w') as OUT:
        # Tab-separated
        if not args.json_format:
            # Construct the rows
            rows = []
            for M, pval in setToPval.iteritems():
                X, T, Z, tbl = setToObs[M]
                row = [ ', '.join(sorted(M)), pval, setToFDR[M], setToRuntime[M], T, Z ] + tbl
                rows.append( row )
            rows.sort(key=lambda row: row[1]) # sort ascending by P-value

            # Create the header
            method_paren = '' if is_permutational else ' ({})'.format(args.method)
            k = len(rows[0][0].split(', '))
            tbl_header = create_tbl_header( k )
            header = 'Gene set\t{0}{1} P-value\t{0}{1} FDR\t{0}{1} '\
                     'Runtime\tT\tZ\t{2}'.format(args.test, method_paren, tbl_header)

            # Output to file
            OUT.write('#{}\n'.format(header))
            OUT.write( '\n'.join([ '\t'.join(map(str, row)) for row in rows ]) )

        # JSON
        else:
            # Record the parameters
            params = vars(args)

            # Output to file
            output = dict(params=params, setToPval=convert_dict_for_json(setToPval),
                          setToObs=convert_dict_for_json(setToObs),
                          setToFDR=convert_dict_for_json(setToFDR),
                          setToRuntime=convert_dict_for_json(setToRuntime))
            json.dump( output, OUT )

# Output MCMC
def output_mcmc(args, setsToFreq, setToPval, setToObs):
    if args.json_format:
        params = vars(args)
        output = dict(params=params, setToPval=convert_dict_for_json(setToPval),
                      setToObs=convert_dict_for_json(setToObs),
                      setsToFreq=dict( (' '.join([ ','.join(sorted(M)) for M in sets ]), freq) for sets, freq in setsToFreq.iteritems() ))
        with open(args.output_prefix + '.json', 'w') as OUT:
            json.dump( output, OUT )
    else:
        # Output a gene set file
        with open(args.output_prefix + '-sampled-collections.tsv', 'w') as OUT:
            rows = []
            for sets, freq in setsToFreq.iteritems():
                row = [ ' '.join([ ','.join(M) for M in sets ]), freq ]
                row.append( sum( -np.log10(setToPval[M] ** args.alpha) for M in sets ))
                rows.append(row)
            rows.sort(key=lambda r: r[1], reverse=True)

            OUT.write('#Gene sets\tSampling Frequency\tCombined Weight\n')
            OUT.write( '\n'.join([ '\t'.join(map(str, row)) for row in rows ]) )

        # Output each of the sample gene sets
        with open(args.output_prefix + '-sampled-sets.tsv', 'w') as OUT:
            rows = []
            for M, pval in setToPval.iteritems():
                X, T, Z, tbl = setToObs[M]
                rows.append([ ','.join(sorted(M)), pval, T, Z] + tbl )
            rows.sort(key=lambda r: r[1])
            k = max(args.gene_set_sizes)

            tbl_header = create_tbl_header( k )
            OUT.write('#Gene set\t{} ({}) P-value\tT\tZ\t{}\n'.format(args.test, args.method, tbl_header))
            OUT.write( '\n'.join([ '\t'.join(map(str, row)) for row in rows ]) )
