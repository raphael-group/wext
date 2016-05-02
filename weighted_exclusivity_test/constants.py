#!/usr/bin/env python

# Set sizes implemented
SET_SIZES_IMPLEMENTED = set([2, 3])

# Sides of the weighted test
ONE_LESS    = 1
ONE_GREATER = 2
TWO         = 3
TAILS       = set([ ONE_LESS, ONE_GREATER, TWO ])
tailConstants = dict(one_less=ONE_LESS, one_greater=ONE_GREATER, two=TWO)

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

