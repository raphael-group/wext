#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json
from collections import defaultdict

# Load the weighted exclusivity test code
current_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.normpath(current_dir + '/../'))
from weighted_exclusivity_test import *

# Argument parser
def get_parser():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-rf', '--results_files', type=str, required=True, nargs='*')
    parser.add_argument('-mf', '--mutation_file', type=str, required=True)
    parser.add_argument('-wf', '--weights_file', type=str, required=False, default=None)
    parser.add_argument('-N', '--num_sets', type=int, required=False, default=None)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    return parser

# Main
def run( args ):
    # Load the input files
    setToFDR, setToPval    = defaultdict(dict), defaultdict(dict)
    setToRuntime, setToObs = defaultdict(dict), defaultdict(dict)
    sets, methods = set(), set()
    for results_file in args.results_files:
        with open(results_file, 'r') as IN:
            obj = json.load(IN)
            params = obj['params']
            min_frequency = params['min_frequency']
            is_permutational = nameToTest[params['test']] == PERMUTATIONAL
            method_paren = '' if is_permutational else ' ({})'.format(params['method'])
            run_name = '{}{}'.format(params['test'], method_paren)
            methods.add( run_name )
            setToPval[run_name].update( obj['setToPval'].items() )
            setToRuntime[run_name].update( obj['setToRuntime'].items() )
            setToFDR[run_name].update( obj['setToFDR'].items() )
            setToObs[run_name].update(obj['setToObs'].items() )
            sets |= set(obj['setToPval'].keys())

    # Load the mutation data
    mutation_data = load_mutation_data( args.mutation_file, min_frequency )
    genes, _, patients, geneToCases, patientToMutations, params, hypermutators = mutation_data
    num_genes, num_patients = len(genes), len(patients)
    geneToIndex = dict(zip(genes, range(num_genes)))
    patientToType = dict( (p, "Hypermutator" if p in hypermutators else "Non-hypermutator")
                          for p in patients )

    # Load the weights
    P = np.load(args.weights_file)
    P = dict( (g, dict(zip(patients, P[geneToIndex[g]]))) for g in genes )

    # Restrict the sets (if necessary)
    if args.num_sets:
        new_sets = set()
        for run_name in methods:
            new_sets |= set(sorted( setToPval[run_name].keys(), key=lambda M: setToPval[M] )[:args.num_sets])

        sets = new_sets
        setToPval       = dict( (run_name, dict( (M, pval) for M, pval in setToPval[run_name].iteritems() if M in new_sets)) for run_name in methods )
        setToRuntime    = dict( (run_name, dict( (M, pval) for M, pval in setToRuntime[run_name].iteritems() if M in new_sets)) for run_name in methods )
        setToObs        = dict( (run_name, dict( (M, pval) for M, pval in setToObs[run_name].iteritems() if M in new_sets)) for run_name in methods )
        setToFDR        = dict( (run_name, dict( (M, pval) for M, pval in setToFDR[run_name].iteritems() if M in new_sets)) for run_name in methods )

    # Restrict the weights
    genes_in_sets = set( g for M in sets for g in M.split('\t') )
    P = dict( (g, P[g]) for g in genes_in_sets )
    geneToCases = dict( (g, cases) for g, cases in geneToCases.iteritems() if g in genes_in_sets )

    print '* Considering {} sets...'.format(len(new_sets))

    # Output the JSON file
    with open(args.output_file, 'w') as OUT:
        # Params
        params = dict(min_freq=min_frequency)
        params['results_files'] = [ os.path.abspath(f) for f in args.results_files ]
        params['mutation_file'] = os.path.abspath(args.mutation_file)
        params['weights_file'] = os.path.abspath(args.weights_file)

        # Output
        output = dict(params=params, geneToCases=dict( (g, list(cases)) for g, cases in geneToCases.iteritems() ),
                      setToPval=setToPval, methods=sorted(methods),
                      patientToType=patientToType, setToFDR=setToFDR,
                      setToRuntime=setToRuntime, setToObs=setToObs, sets=list(sets),
                      genes=list(genes), patients=patients, P=P)
        json.dump( output, OUT )

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
