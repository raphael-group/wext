#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, matplotlib.pyplot as plt, numpy as np, json
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-mf', '--mutation_files', type=str, required=True, nargs='*')
parser.add_argument('-wf', '--weights_files', type=str, required=True, nargs='*')
parser.add_argument('-c', '--cancers', type=str, required=True, nargs='*')
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

assert( len(args.cancers) == len(args.weights_files) == len(args.mutation_files) )

# Load the mutation file
print '* Loading mutation files...'
cancerToWeights, cancerToPatients, cancerToGenes, cancerToHypermutators, patientToMutations, geneToCases = dict(), dict(), dict(), dict(), dict(), dict()
for cancer, mutation_file, weights_file in zip(args.cancers, args.mutation_files, args.weights_files):
    with open(mutation_file, 'r') as IN:
        obj = json.load(IN)
        num_genes, num_patients = len(obj['genes']), len(obj['patients'])
        cancerToPatients[cancer] = obj['patients']
        cancerToGenes[cancer] = obj['genes']
        cancerToHypermutators[cancer] = set(obj['hypermutators'])
        geneToCases[cancer] = obj['geneToCases']
        patientToMutations[cancer] = dict( (p, set()) for p in obj['patients'] )
        for g, cases in geneToCases[cancer].iteritems():
            for p in cases:
                patientToMutations[cancer][p].add( g )
    cancerToWeights[cancer] = np.load(weights_file)
    print '\t{}\n\t\t- Genes: {}\n\t\t- Patients: {}'.format(cancer, num_genes, num_patients)

# Set up the figure
fig, axes = plt.subplots( 1, len(args.cancers))
fig.set_size_inches( len(args.cancers) * 5, 5)
min_weight = min([ np.min(W) for W in cancerToWeights.values() ])
print 'Min weight:', min_weight
for ax, cancer in zip(axes, args.cancers):
    # Sort the weights so that hypermutators are all on one side
    patients = cancerToPatients[cancer]
    genes = cancerToGenes[cancer]
    hypermutators = cancerToHypermutators[cancer]
    num_non_hypermutators = len(set(patients) - hypermutators)
    patient_indices = sorted(range(len(patients)), key=lambda p: (patients[p] in hypermutators, len(patientToMutations[cancer][patients[p]])))
    gene_indices = sorted([ i for i, g in enumerate(genes) if g in geneToCases[cancer]], key=lambda g: len(geneToCases[cancer].get(genes[g], [])), reverse=True)
    weights = [ row[patient_indices] for row in cancerToWeights[cancer][gene_indices] ]

    # Plot the weights matrix
    img = ax.matshow(weights, norm=matplotlib.colors.LogNorm(), vmin=min_weight, vmax=1)
    ax.plot([num_non_hypermutators, num_non_hypermutators], [0, len(weights)], '--', c='k')
    ax.set_xlabel('Samples')
    ax.set_ylabel('Genes')
    ax.set_xlim([0, len(weights[0])])
    ax.set_ylim([0, len(weights)])
    ax.set_title(cancer)
    ax.set_aspect('auto')
    ax.tick_params(labelsize=6)

# Add a colorbar
plt.tight_layout()
fig.subplots_adjust(right=0.88)
cbar_ax = fig.add_axes([0.90, 0.08, 0.05, 0.80])
fig.colorbar(img, cax=cbar_ax, format=matplotlib.ticker.LogFormatter(10, labelOnlyBase=False))

plt.savefig(args.output_file)
