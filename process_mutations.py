#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, numpy as np
from collections import defaultdict

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--maf_file', type=str, required=True)
    parser.add_argument('-hf', '--hypermutators_file', type=str, required=False, default=None)
    parser.add_argument('-ivc', '--ignored_variant_classes', type=str, required=False, nargs='*',
                        default=["Silent", "Intron", "3'UTR", "5'UTR", "IGR", "lincRNA"])
    parser.add_argument('-ivt', '--ignored_variant_types', type=str, required=False, nargs='*',
                        default=['Germline'])
    parser.add_argument('-ivs', '--ignored_validation_statuses', type=str, required=False, nargs='*',
                        default=['Wildtype', 'Invalid'])
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, default=1, required=False, choices=range(5))
    return parser

def run( args ):
    # Load the mutations from the MAF
    if args.verbose > 0: print '* Loading the MAF...'
    with open(args.maf_file, 'r') as IN:
        geneToCases, patientToMutations = defaultdict( set ), defaultdict( set )
        seenHeader = False
        ignored_var_classes  = set(map(str.lower, args.ignored_variant_classes))
        ignored_var_types    = set(map(str.lower, args.ignored_variant_types))
        ignored_val_statuses = set(map(str.lower, args.ignored_validation_statuses))
        genes, patients, var_classes, var_types, val_statuses = set(), set(), set(), set(), set()
        for l in IN:
            arr = l.rstrip('\n').split('\t')
            # Parse the header if we haven't seen it yet
            if not seenHeader and arr[0].lower() == 'hugo_symbol':
                arr              = map(str.lower, arr)
                seenHeader       = True
                gene_index       = 0
                var_class_index  = arr.index('variant_classification')
                var_type_index   = arr.index('variant_type')
                patient_index    = arr.index('tumor_sample_barcode')
                val_status_index = arr.index('validation_status')
            # Process the mutation
            else:
                # Record the patients and genes, even if we ignore their mutations
                patient, gene = '-'.join(arr[patient_index].split('-')[:3]), arr[gene_index]
                
                patients.add(patient)
                genes.add(gene)
                
                # Ignore certain mutations
                var_classes.add( arr[var_class_index].lower())
                var_types.add( arr[var_type_index].lower())
                val_statuses.add( arr[val_status_index].lower() )
                if arr[var_class_index].lower() in ignored_var_classes \
                   or arr[var_type_index].lower() in ignored_var_types \
                   or arr[val_status_index].lower() in ignored_val_statuses:
                    continue
                
                # Record the mutation
                geneToCases[gene].add( patient )
                patientToMutations[patient].add( gene )
                
    patients     = sorted(patients)
    genes        = sorted(genes)
    num_patients = len(patients)
    num_genes    = len(genes)

    # Load the hypermutators    
    if args.hypermutators_file:
        with open( args.hypermutators_file, 'r') as IN:
            hypermutators = set( l.rstrip('\n') for l in IN if not l.startswith('#') )
    else:
        hypermutators = set()

    # Summarize the data
    if args.verbose > 0:
        print '\tGenes: {}'.format(num_genes)
        print '\tPatients: {} ({} hypermutators)'.format(num_patients, len(hypermutators))
        print '\tVariant classes:', ', '.join(sorted(var_classes))
        print '\tVariant types:', ', '.join(sorted(var_types))
        print '\tValidation statuses:', ', '.join(sorted(val_statuses))
    
    # Output to file
    with open(args.output_file, 'w') as OUT:
        params = dict(maf_file=os.path.abspath(args.maf_file),
                      hypermutators_file=os.path.abspath(args.hypermutators_file) if args.hypermutators_file else None)
        output = dict(params=params, patients=patients, genes=genes, hypermutators=list(hypermutators),
                      geneToCases=dict( (g, list(cases)) for g, cases in geneToCases.items()),
                      patientToMutations=dict( (p, list(muts)) for p, muts in patientToMutations.items()),
                      num_genes=num_genes, num_patients=num_patients)
        json.dump( output, OUT )
    
if __name__ == '__main__': run( get_parser().parse_args( sys.argv[1:]) )
