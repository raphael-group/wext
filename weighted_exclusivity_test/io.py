#!/usr/bin/env python

import os, json
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
def output_table(args, setToPval, setToRuntime, setToFDR, setToObs ):
    with open(args.output_file, 'w') as OUT:
        # Tab-separated
        if args.output_file.lower().endswith('.tsv'):
            # Construct the rows
            rows = []
            for M, pval in setToPval.iteritems():
                T, Z, tbl = setToObs[M]
                row = [ ', '.join(sorted(M)), pval, setToFDR[M], setToRuntime[M], T, Z ] + tbl
                rows.append( row )
            rows.sort(key=lambda row: row[1]) # sort ascending by P-value
            
            # Create the header
            k = len(rows[0][0].split(', '))
            tbl_header = create_tbl_header( k )
            header = 'Gene set\t{0} ({1}) P-value\t{0} ({1}) FDR\t{0} ({1}) '\
                     'Runtime\tT\tZ\t{2}'.format(args.test, args.method, tbl_header)
                     
            # Output to file
            OUT.write('#{}\n'.format(header))
            OUT.write( '\n'.join([ '\t'.join(map(str, row)) for row in rows ]) )
            
        # JSON
        else:
            params = dict(method=args.method, test=args.test, min_frequency=args.min_frequency,
                          mutation_file=os.path.abspath(args.mutation_file),
                          k=args.gene_set_size, weights_file=os.path.abspath(args.weights_file))
            output = dict(params=params, setToPval=convert_dict_for_json(setToPval),
                          setToObs=convert_dict_for_json(setToObs),
                          setToFDR=convert_dict_for_json(setToFDR),
                          setToRuntime=convert_dict_for_json(setToRuntime))
            json.dump( output, OUT )

    
