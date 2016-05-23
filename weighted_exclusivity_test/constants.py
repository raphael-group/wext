#!/usr/bin/env python

# P-values are called invalid if P > 1+PTOL or P < -PTOL
PTOL = 10**-3

# Set sizes implemented
WEIGHTED_SADDLEPOINT_SET_SIZES_IMPLEMENTED = set([2, 3])

# Different methods for computing the P-value
EXACT        = 1
SADDLEPOINT  = 2

nameToMethod = dict(Exact=EXACT, Saddlepoint=SADDLEPOINT)
methodToName = { EXACT: "Exact", SADDLEPOINT: "Saddlepoint" }
METHODS      = set(methodToName.keys())
METHOD_NAMES = set(methodToName.values())

WEIGHTED      = 1
UNWEIGHTED    = 2
PERMUTATIONAL = 3
testToName    = { UNWEIGHTED: "Unweighted", WEIGHTED: "Weighted", PERMUTATIONAL: "Permutational"}
nameToTest    = dict(Unweighted=UNWEIGHTED, Weighted=WEIGHTED, Permutational=PERMUTATIONAL)
TESTS         = set(testToName.keys())
TEST_NAMES    = set(testToName.values())
