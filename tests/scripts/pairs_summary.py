#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, json, matplotlib.pyplot as plt, numpy as np, pandas as pd, seaborn as sns
from collections import defaultdict
from scipy.stats import spearmanr
plt.style.use('ggplot')

# Load the weighted exclusivity test
sys.path.append(os.getcwd())
from weighted_exclusivity_test import *
from helper import add_y_equals_x

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pf', '--pairs_files', type=str, required=True, nargs='*')
parser.add_argument('-c', '--cancers', type=str, required=True, nargs='*')
parser.add_argument('-np', '--num_permutations', type=int, required=True)
parser.add_argument('-ff', '--figure_file', type=str, required=True)
parser.add_argument('-tf', '--table_files', type=str, required=True, nargs=2)
args = parser.parse_args( sys.argv[1:] )

# Load the pairs files
setToPval, setToRuntime, sets = defaultdict(dict), defaultdict(dict), defaultdict(set)
for cancer, pairs_file in zip(args.cancers, args.pairs_files):
    with open(pairs_file, 'r') as IN:
        obj = json.load(IN)
        params = obj['params']
        is_permutational = nameToTest[params['test']] == PERMUTATIONAL
        method_paren = '' if is_permutational else ' ({})'.format(params['method'])
        run_name = '{}{}'.format(params['test'], method_paren)
        setToPval[cancer][run_name] = obj['setToPval']
        setToRuntime[cancer][run_name] = obj['setToRuntime']
        sets[cancer] |= set(obj['setToPval'].keys())

# Construct the data
cancers = sorted(set(args.cancers))
runs = sorted(setToRuntime[cancers[0]].keys())
comet_pvals = []
weighted_exact_pvals = []
weighted_saddlepoint_pvals = []
permutational_pvals = []
comet_nt_pvals = [] #nt = not tail
weighted_exact_nt_pvals = []
weighted_saddlepoint_nt_pvals = []
permutational_nt_pvals = []
items = []
for cancer in cancers:
    for M in sets[cancer]:
        # Record the P-values (setting zeros in the permutational appropriately)
        if not all(M in setToPval[cancer][r] for r in runs ): continue
        comet_pvals.append(setToPval[cancer]["Unweighted (Exact)"][M])
        weighted_exact_pvals.append(setToPval[cancer]["Weighted (Exact)"][M])
        weighted_saddlepoint_pvals.append(setToPval[cancer]["Weighted (Saddlepoint)"][M])

        pval = setToPval[cancer]["Permutational"][M]
        if pval == 0:
            pval = 1./args.num_permutations
        else:
            comet_nt_pvals.append(comet_pvals[-1])
            weighted_exact_nt_pvals.append(weighted_exact_pvals[-1])
            weighted_saddlepoint_nt_pvals.append(weighted_saddlepoint_pvals[-1])
            permutational_nt_pvals.append(pval)
        permutational_pvals.append(pval)

        # Record the runtimes
        for method in ["Weighted (Exact)", "Weighted (Saddlepoint)", "Unweighted (Exact)"]:
            items.append({ "Method": method, "Runtime (secs)": setToRuntime[cancer][method][M],
                        "Cancer": cancer})
df = pd.DataFrame(items)

print 'Testing {} pairs...'.format(len(weighted_exact_pvals))

# Set up the figure
fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4)
fig.set_size_inches(20, 5)

# Plot the permutational against CoMEt (eq. to Fisher's exact test)
ax1.plot(permutational_pvals, comet_pvals, 'o', c='r', alpha=0.5, mew=0.0)
ax1.set_xlabel("Permutational P-value (N={})".format(args.num_permutations))
ax1.set_ylabel("Unweighted P-value (Fisher's exact test)")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title('(a)')
add_y_equals_x(ax1)

# Plot the permutational against the weighted exact test
ax2.plot(permutational_pvals, weighted_exact_pvals, 'o', c='b', alpha=0.5, mew=0.0)
ax2.set_xlabel("Permutational P-value (N={})".format(args.num_permutations))
ax2.set_ylabel("Weighted P-value (exact test)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title('(b)')
add_y_equals_x(ax2)

# Plot the weighted exact test against the
ax3.plot(weighted_exact_pvals, weighted_saddlepoint_pvals, 'o', c='m', alpha=0.5, mew=0.0)
ax3.set_xlabel("Weighted P-value (exact test)")
ax3.set_ylabel("Weighted P-value (saddlepoint)")
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_title('(c)')
add_y_equals_x(ax3)

# Plot the runtimes of the saddlepoint vs. weighted exact test
df_weighted = df.loc[df['Method'].isin(set(['Weighted (Exact)', 'Weighted (Saddlepoint)']))]
sns.boxplot(x="Cancer", y="Runtime (secs)", hue="Method", data=df_weighted, ax=ax4)
ax4.set_yscale('log')
ax4.set_title('(d)')
ax4.legend(loc='upper left')

# Output a table sumamrizing the correlations (Table 2)
with open(args.table_files[0], 'w') as OUT:
    OUT.write('#Pairs\tCoMEt\tWeighted Exact\tWeighted Saddlepoint\n')
    OUT.write('All\t{}\t{}\t{}\n'.format(spearmanr(permutational_pvals, comet_pvals)[0], spearmanr(permutational_pvals, weighted_exact_pvals)[0], spearmanr(permutational_pvals, weighted_saddlepoint_pvals)[0]))
    OUT.write('H(M) >= 10^-4\t{}\t{}\t{}\n'.format(spearmanr(permutational_nt_pvals, comet_nt_pvals)[0], spearmanr(permutational_nt_pvals, weighted_exact_nt_pvals)[0],
                                                   spearmanr(permutational_nt_pvals, weighted_saddlepoint_nt_pvals)[0]))

# Output a table summarizing the runtimes (Table 3)
tbl = ['#Method\tMinimum\tMedian\tMaximum\tTotal']
for method in ["Weighted (Exact)", "Weighted (Saddlepoint)"]:
    runtimes = list(df.loc[df['Method'] == method]['Runtime (secs)'])
    minimum, median, maximum = np.min(runtimes), np.median(runtimes), np.max(runtimes)
    total = np.sum(runtimes)
    tbl.append('{}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}'.format(method, minimum, median, maximum, total))

with open(args.table_files[1], 'w') as OUT:
    OUT.write('\n'.join(tbl) )

# Output to file
plt.tight_layout()
plt.savefig(args.figure_file)
