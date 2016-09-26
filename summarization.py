#!/usr/bin/env python

# Load required modules
import sys, os, argparse, networkx as nx
from math import log10
from itertools import combinations
from collections import defaultdict

# Argument parser
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True)
    parser.add_argument('-w', '--weight_column', type=int, required=False, default=1)
    parser.add_argument('-t', '--threshold', type=float, required=False, default=None)
    parser.add_argument('-N', '--num_top_edges', type=int, required=True)
    parser.add_argument('--negative-log10', action='store_true', required=False, default=False, help='Flag to take negative log10 of weight column.')
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, choices=range(5), required=False, default=1)
    return parser

def cut_graph():
    return 

def run( args ):
    # Load the input file
    if args.verbose > 0: print '* Loading input file...'
    edges = defaultdict( float )
    with open(args.input_file, 'r') as IN:
        for i, l in enumerate(IN):
            if l.startswith('#'): continue
            
            # Parse theline
            arr = l.rstrip('\n').split('\t')
            M   = arr[0].split(', ')

            # Compute the weight, stopping our iteration if the
            # weight is above our threshold
            if args.negative_log10:
                weight = -log10(float(arr[1]))
                if args.threshold and weight < args.threshold:
                    break
            else:
                weight = float(arr[1])
                if args.threshold and weight > args.threshold:
                    break

            # Add the weight for each 
            for u, v in combinations(M, 2):
                edges[frozenset([u, v])] += weight

    # Print summarization
    nodes = set( n for e in edges for n in e )
    if args.verbose > 0:
        print '\t- Nodes: {}'.format(len(nodes))
        print '\t- Edges: {}'.format(len(edges))
                
    # Output the summarization
    rows = [ (u, v, w) for (u, v), w in edges.iteritems() ]
    rows.sort(key=lambda (u, v, w): w, reverse=True)
    with open(args.output_file, 'w') as OUT:
        OUT.write('#Gene 1\tGene 2\tWeight\n')
        OUT.write( '\n'.join( '\t'.join(map(str, row)) for row in rows ))

    # Print connected components with top N edges
    G = nx.Graph()
    G.add_edges_from( [ (u, v) for u, v, w in rows[:args.num_top_edges] ] )
    print '\n'.join( ', '.join(sorted(M)) for M in nx.connected_components(G) )

if __name__ == '__main__': run( get_parser().parse_args( sys.argv[1:] ))
