#!/usr/bin/env bash

num_permutations=1000
num_cores=4

# Preprocess mutations.
python ../../process_mutations.py \
    -m  adjacency_list.tsv \
    -ct NA \
    -o  data.json

# Compute mutation probabilities.
python ../../compute_mutation_probabilities.py \
    -mf data.json \
    -np $num_permutations \
    -nc $num_cores \
    -wf weights.npy \
    -s  12345 \
    -v  1

# Find sets using mutual exclusivity test statistic.
python ../../find_sets.py \
    -mf data.json \
    -wf weights.npy \
    -s  exclusivity \
    -k  2 \
    -c  $num_cores \
    -f  2 \
    -o  exclusivity_results \
    -v  0

# Find sets using a co-occurrence test statistic (any co-occurrence).
python ../../find_sets.py \
    -mf data.json \
    -wf weights.npy \
    -s  any-co-occurrence \
    -k  2 \
    -c  $num_cores \
    -f  2 \
    -o  any-co-occurrence_results \
    -v  0

# Find sets using another co-occurrence test statistic (all co-occurrence).
python ../../find_sets.py \
    -mf data.json \
    -wf weights.npy \
    -s  all-co-occurrence \
    -k  2 \
    -c  $num_cores \
    -f  2 \
    -o  all-co-occurrence_results \
    -v  0
