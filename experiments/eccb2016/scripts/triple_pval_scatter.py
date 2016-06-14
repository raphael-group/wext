#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, matplotlib.pyplot as plt, numpy as np, json
from scipy.stats import spearmanr
from helper import add_y_equals_x
plt.style.use('ggplot')

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-uwf', '--unweighted_files', type=str, required=True, nargs='*')
parser.add_argument('-wf', '--weighted_files', type=str, required=True, nargs='*')
parser.add_argument('-pf', '--permutational_files', type=str, required=True, nargs='*')
parser.add_argument('-np', '--num_permutations', type=int, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

if args.num_permutations < 1e6:
    sys.stderr.write('Warning: Need at least 1 million permutations to accurately estimate correlations\n')

# Load the sets and their p-values
setToUnweighted = dict()
for unweighted_file in args.unweighted_files:
    with open(unweighted_file, 'r') as IN:
        setToUnweighted.update( json.load(IN)['setToPval'] )

setToWeighted = dict()
for weighted_file in args.weighted_files:
    with open(weighted_file, 'r') as IN:
        setToWeighted.update( json.load(IN)['setToPval'] )

setToPermuted = dict()
for permuted_file in args.permutational_files:
    with open(permuted_file, 'r') as IN:
        setToPermuted.update( json.load(IN)['setToPval'] )

for M, pval in setToPermuted.iteritems():
    if pval == 0:
        setToPermuted[M] = 1./args.num_permutations

sets = set(setToWeighted.keys()) & set(setToUnweighted.keys()) & set(setToPermuted.keys())

print '* Loaded weighted/unweighted P-values in {} triples...'.format(len(setToWeighted))
print '\t- Weighted range: [{}, {}]'.format(np.min(setToWeighted.values()), np.max(setToWeighted.values()))
print '\t- Unweighted range: [{}, {}]'.format(np.min(setToUnweighted.values()), np.max(setToUnweighted.values()))
print '* Loaded permuted P-values for {} sets ({} intersection)...'.format(len(setToPermuted), len(sets))

# Create two scatter plots
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(10, 5)

# First, permutation vs. unweighted
xs = [ setToPermuted[t] for t in sets ]
ys = [ setToUnweighted[t] for t in sets ]
unweighted_rho = spearmanr(xs, ys)[0]
unweighted_tail_rho = spearmanr([ x for x in xs if 1e-5 < x < 1e-3 ], [y for x, y in zip(xs,ys) if 1e-5 < x < 1e-3 ])[0]
ax1.plot( xs, ys, 'o', c='b', alpha=0.5, mew=0.0, ms=3)
ax1.set_xlabel('RC-exclusivity $p$-value (permutations, $N={}$)'.format(args.num_permutations))
ax1.set_ylabel('R-exclusivity $p$-value (tail enumeration/CoMEt)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([1./(10*args.num_permutations), 1.])
ax1.plot(ax1.get_xlim(), ax1.get_xlim(), ls="--", c=".3")

# Second, permutation vs weighted
ys = [ setToWeighted[t] for t in sets ]
weighted_rho = spearmanr(xs, ys)[0]
weighted_tail_rho = spearmanr([ x for x in xs if 1e-5 < x < 1e-3 ], [y for x, y in zip(xs,ys) if 1e-5 < x < 1e-3 ])[0]
ax2.plot( xs, ys, 'o', c='r', alpha=0.5, mew=0.0, ms=3)
ax2.set_xlabel('RC-exclusivity $p$-value (permutations, $N={}$)'.format(args.num_permutations))
ax2.set_ylabel('WR-exclusivity $p$-value (saddlepoint approximation)')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([1./(10*args.num_permutations), 1.])
ax2.plot(ax2.get_xlim(), ax2.get_xlim(), ls="--", c=".3")

# Output maximum deviation and correlations
print 'Max deviation permutational vs. weighted (1E-3 to 1E-5):',
deviations = [ (x, y, np.abs(y/x)) for x, y in zip(xs, ys) if 1e-3 > x > 1e-5 ]
if deviations:
    print max(deviations, key=lambda (x, y, z): z)
else:
    print 'None in p-value interval'

print 'Unweighted correlation (all): \\rho={}'.format(unweighted_rho)
print 'Unweighted correlation (P<0.001): \\rho={}'.format(unweighted_tail_rho)
print 'Weighted correlation (all): \\rho={}'.format(weighted_rho)
print 'Weighted correlation (P<0.001): \\rho={}'.format(weighted_tail_rho)

# Output to file
plt.tight_layout()
plt.savefig(args.output_file)
