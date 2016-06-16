#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, numpy as np
from collections import defaultdict

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutation_file_groups', type=str, required=True,
                        action='append', nargs='*')
    parser.add_argument('-ct', '--cancer_types', type=str, required=True, nargs='*')
    parser.add_argument('-pw', '--patient_whitelist', type=str, required=False)
    parser.add_argument('-hf', '--hypermutators_file', type=str, required=False, default=None)
    parser.add_argument('-ivc', '--ignored_variant_classes', type=str, required=False, nargs='*',
                        default=["Silent", "Intron", "3'UTR", "5'UTR", "IGR", "lincRNA", "RNA"])
    parser.add_argument('-ivt', '--ignored_variant_types', type=str, required=False, nargs='*',
                        default=['Germline'])
    parser.add_argument('-ivs', '--ignored_validation_statuses', type=str, required=False, nargs='*',
                        default=['Wildtype', 'Invalid'])
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, default=1, required=False, choices=range(5))
    return parser

def process_maf( maf_file, patientWhitelist, geneToCases, patientToMutations, vc, vt, vs, ivc, ivt, ivs, verbose ):
    if verbose > 1: print '\tLoading MAF:', maf_file
    genes, patients = set(), set()
    with open(maf_file, 'r') as IN:
        seenHeader = False
        for i, l in enumerate(IN):
            arr = l.rstrip('\n').split('\t')
            # Parse the header if we haven't seen it yet
            if not seenHeader and arr[0].lower() == 'hugo_symbol':
                arr              = map(str.lower, arr)
                seenHeader       = True
                gene_index       = 0
                patient_index    = arr.index('tumor_sample_barcode')
                var_class_index  = arr.index('variant_classification')  if 'variant_classification' in arr else None
                var_type_index   = arr.index('variant_type')  if 'variant_type' in arr else None
                val_status_index = arr.index('validation_status') if 'validation_status' in arr else None

            # Process the mutation
            elif seenHeader:
                # Record the patients and genes, even if we ignore their mutations
                patient, gene = '-'.join(arr[patient_index].split('-')[:3]), arr[gene_index]

                if not patientWhitelist[patient]: continue

                patients.add(patient)
                genes.add(gene)

                # Ignore certain mutations. If we aren't given certain information, we
                # will include the mutation by default
                if var_class_index:
                    var_class = arr[var_class_index].lower()
                    vc_check = not (var_class in ivc)
                else:
                    var_class = ''
                    vc_check = True

                if var_type_index:
                    var_type = arr[var_type_index].lower()
                    vt_check = not (var_type in ivt)
                else:
                    var_type = ''
                    vt_check = True

                if val_status_index:
                    val_status = arr[val_status_index].lower()
                    vs_check = not (val_status in ivs)
                else:
                    val_status = ''
                    vs_check = True

                # Record the mutation
                if vc_check and vt_check and vs_check:
                    vt.add( var_type )
                    vs.add( val_status )
                    vc.add( var_class )
                    geneToCases[gene].add( patient )
                    patientToMutations[patient].add( gene )

    return genes, patients

def process_events_file( events_file, patientWhitelist, geneToCases, patientToMutations, verbose ):
    if verbose > 1: print '\tProcessing events file:', events_file

    # Parse the events file
    events, patients = set(), set()
    with open(events_file, 'r') as IN:
        arrs = [ l.rstrip('\n').split('\t') for l in IN if not l.startswith('#') ]
        for arr in arrs:
            # Skip patients that aren't whitelisted
            patient, mutations = arr[0], set(arr[1:])
            if not patientWhitelist[patient]: continue

            # Record the events and mutations
            patients.add(patient)
            events |= mutations
            patientToMutations[patient] |= mutations
            for event in mutations:
                geneToCases[event].add( patient )

    return events, patients

def run( args ):
    # Do some additional argument checking
    assert( len(args.mutation_file_groups) == len(args.cancer_types) )
    ivc = set(map(str.lower, args.ignored_variant_classes))
    ivt = set(map(str.lower, args.ignored_variant_types))
    ivs = set(map(str.lower, args.ignored_validation_statuses))

    # Load the patient whitelist (if supplied)
    if args.patient_whitelist:
        if args.verbose > 0: print '* Loading patient whitelist...'
        patientWhitelist = defaultdict( lambda : False )
        with open(args.patient_whitelist, 'r') as IN:
            patientWhitelist.update( (l.rstrip('\n').split()[0], True) for l in IN if not l.startswith('#') )
    else:
        if args.verbose > 0: print '* No patient whitelist provided, including all patients...'
        patientWhitelist = defaultdict( lambda : True )

    # Load the mutations from each MAF
    if args.verbose > 0: print '* Loading and combining {} datasets...'.format(len(args.cancer_types))
    geneToCases, patientToMutations = defaultdict( set ), defaultdict( set )
    genes, patients = set(), set()
    vc, vt, vs = set(), set(), set() # variant classes/types and validation statuses
    patientToType = dict()
    cancerToFiles = defaultdict(list)

    for cancer_type, mutation_file_group in zip(args.cancer_types, args.mutation_file_groups):
        for mutation_file in mutation_file_group:
            cancerToFiles[cancer_type].append(os.path.abspath(mutation_file))
            if mutation_file.lower().endswith('.maf'):
                per_type_genes, per_type_patients = process_maf( mutation_file, patientWhitelist, geneToCases, patientToMutations, vc, vt, vs, ivc, ivt, ivs, args.verbose)
            else:
                per_type_genes, per_type_patients = process_events_file( mutation_file, patientWhitelist, geneToCases, patientToMutations, args.verbose )

            patientToType.update( (p, cancer_type) for p in per_type_patients )
            patients |= per_type_patients
            genes    |= per_type_genes

    patients     = sorted(patients)
    genes        = sorted(genes)
    num_patients = len(patients)
    num_genes    = len(genes)

    # Load the hypermutators/
    if args.hypermutators_file:
        with open( args.hypermutators_file, 'r') as IN:
            hypermutators = set( l.rstrip('\n') for l in IN if not l.startswith('#') )
    else:
        hypermutators = set()

    # Summarize the data
    if args.verbose > 0:
        print '* Summary of mutation data...'
        print '\tGenes: {}'.format(num_genes)
        print '\tPatients: {} ({} hypermutators)'.format(num_patients, len(hypermutators))
        print '\tUsed variant classes:', ', '.join(sorted(vc))
        print '\tUsed variant types:', ', '.join(sorted(vt))
        print '\tUsed validation statuses:', ', '.join(sorted(vs))

    # Output to file
    with open(args.output_file, 'w') as OUT:
        params = dict(cancerToFiles=cancerToFiles,
                      cancer_types=args.cancer_types,
                      ignored_variant_classes=args.ignored_variant_classes,
                      ignored_variant_types=args.ignored_variant_types,
                      ignored_validation_statuses=args.ignored_validation_statuses,
                      patient_whitelist_file=os.path.abspath(args.patient_whitelist) if args.patient_whitelist else None,
                      hypermutators_file=os.path.abspath(args.hypermutators_file) if args.hypermutators_file else None)
        output = dict(params=params, patients=patients, genes=genes, hypermutators=list(hypermutators),
                      geneToCases=dict( (g, list(cases)) for g, cases in geneToCases.items()),
                      patientToType=patientToType,
                      patientToMutations=dict( (p, list(muts)) for p, muts in patientToMutations.items()),
                      num_genes=num_genes, num_patients=num_patients)
        json.dump( output, OUT )

if __name__ == '__main__': run( get_parser().parse_args( sys.argv[1:]) )
