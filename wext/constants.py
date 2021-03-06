#!/usr/bin/env python

# P-values are called invalid if P > 1+PTOL or P < -PTOL
PTOL = 10**-3

# Set sizes implemented
WRE_EXACT_SET_SIZES_IMPLEMENTED = set([2, 3])

# Different methods for computing the P-value
EXACT        = 1
SADDLEPOINT  = 2

nameToMethod = dict(Exact=EXACT, Saddlepoint=SADDLEPOINT)
methodToName = { EXACT: "Exact", SADDLEPOINT: "Saddlepoint" }
METHODS      = set(methodToName.keys())
METHOD_NAMES = set(methodToName.values())

RE         = 1
RCE        = 2
WRE        = 3
testToName = { RE: "RE", WRE: "WRE", RCE: "RCE"}
nameToTest = dict(RE=RE, WRE=WRE, RCE=RCE)
TESTS      = set(testToName.keys())
TEST_NAMES = set(testToName.values())

EXCLUSIVITY       = 'exclusivity'
ANY_CO_OCCURRENCE = 'any-co-occurrence'
ALL_CO_OCCURRENCE = 'all-co-occurrence'
nameToStatistic   = { EXCLUSIVITY: 'exclusivity', ANY_CO_OCCURRENCE: 'any-co-occurrence', ALL_CO_OCCURRENCE: 'all-co-occurrence'}
statisticToName   = dict(EXCLUSIVITY=EXCLUSIVITY, ANY_CO_OCCURRENCE=ANY_CO_OCCURRENCE, ALL_CO_OCCURRENCE=ALL_CO_OCCURRENCE)
STATISTICS        = set(statisticToName.keys())
STATISTIC_NAMES   = set(statisticToName.values())
