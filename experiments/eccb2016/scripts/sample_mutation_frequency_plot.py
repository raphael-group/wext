#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, json, matplotlib.pyplot as plt, numpy as np, seaborn as sns, pandas as pd
from ragrpyfu.io import aligned_plaintext_table
plt.style.use('ggplot')

# Parse required arguments
parser = argparse.ArgumentParser()
parser.add_argument('-mf', '--mutation_files', type=str, required=True, nargs='*')
parser.add_argument('-c', '--cancers', type=str, required=True, nargs='*')
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

assert( len(args.cancers) == len(args.mutation_files) )

# Load the data into a data frame
items = []
for mutation_file, cancer in zip(args.mutation_files, args.cancers):
    with open(mutation_file, 'r') as IN:
        obj = json.load(IN)
        patients = set(obj['patients'])
        hypermutators = set(obj['hypermutators'])

        # Make a map of patients to their mutated genes
        patientToMutations = dict( (p, set()) for p in patients )
        for g, cases in obj['geneToCases'].iteritems():
            for p in cases:
                patientToMutations[p].add( g )

        # Assemble the data into dictionaries for Pandas
        for p, mutations in patientToMutations.iteritems():
            ty = "Hypermutator" if p in hypermutators else "Non-hypermutator"
            items.append({ "Patient": p, "Mutated genes": len(mutations), "Type": ty, "Cancer": cancer })

df = pd.DataFrame(items)

# Create the plot
sns.boxplot(x="Cancer", y="Mutated genes", hue="Type", data=df)
plt.yscale('log')
plt.savefig(args.output_file)

# Output table
tbl = [['Cancer', 'All', 'Hypermutators', 'Non-hypermutators']]
for c in args.cancers:
    all_rates = list(df.loc[df['Cancer'] == c]['Mutated genes'])
    hyper_rates = list(df.loc[(df['Cancer'] == c) & (df['Type'] == "Hypermutator")]['Mutated genes'])
    non_hyper_rates = list(df.loc[(df['Cancer'] == c) & (df['Type'] == "Non-hypermutator")]['Mutated genes'])
    tbl.append([ c, np.median(all_rates), np.median(hyper_rates) if len(hyper_rates) > 0 else '--', np.median(non_hyper_rates)])

print aligned_plaintext_table('\n'.join([ '\t'.join(map(str, row)) for row in tbl ]))
