#!/bin/sh

################################################################################
# SETTINGS                                                                     #
################################################################################
NUM_CORES=25
COADREAD_MIN_FREQ=20
THCA_MIN_FREQ=5
UCEC_MIN_FREQ=30
LENGTH_THRESHOLD=600
FDR_CUTOFF=0.001
TOTAL_PERMUTATIONS=10000 # set to 1000000 for ECCB 2016 submission
PAIRS_PERMUTATIONS=10000
WEIGHTS_PERMUTATIONS=1000

################################################################################
# SET UP DIRECTORIES                                                           #
################################################################################
CODE_DIR=../../
EXPERIMENT_DIR=experiments/eccb2016
DATA_DIR=$EXPERIMENT_DIR/data
SCRIPTS_DIR=$EXPERIMENT_DIR/scripts
FIGURES_DIR=$EXPERIMENT_DIR/figures
TABLES_DIR=$EXPERIMENT_DIR/tables
MUTATIONS_DIR=$DATA_DIR/mutations
PERMUTED_DIR=$DATA_DIR/permuted
WEIGHTS_DIR=$DATA_DIR/weights

OUTPUT_DIR=$EXPERIMENT_DIR/output
PAIRS_DIR=$OUTPUT_DIR/pairs
TRIPLES_DIR=$OUTPUT_DIR/triples
GENE_LENGTH_FILE=$DATA_DIR/gene-lengths.tsv

cd $CODE_DIR
mkdir -p $MUTATIONS_DIR $WEIGHTS_DIR $PERMUTED_DIR $OUTPUT_DIR $PAIRS_DIR
mkdir -p $TRIPLES_DIR $TABLES_DIR $FIGURES_DIR

################################################################################
# DOWNLOAD AND UNPACK DATA                                                     #
################################################################################
cd $EXPERIMENT_DIR
wget http://compbio-research.cs.brown.edu/projects/weighted-exclusivity-test/data/eccb2016.tar
tar -xvf eccb2016.tar && rm eccb2016.tar

################################################################################
# PRE-PROCESS MUTATIONS                                                        #
################################################################################
cd $CODE_DIR

# Colorectal
COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/coadread-mutations-wout-length-filtering.json
python process_mutations.py -mf $DATA_DIR/mafs/COADREAD.maf \
       -hf $DATA_DIR/sample_lists/coadread-hypermutators.txt \
       -o $COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING

COADREAD_MUTATIONS=$MUTATIONS_DIR/coadread-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $COADREAD_MUTATIONS

# # Thyroid
THCA_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/thca-mutations-wout-length-filtering.json
python process_mutations.py -mf $DATA_DIR/mafs/THCA.maf \
       -o $THCA_MUTATIONS_WOUT_LENGTH_FILTERING

THCA_MUTATIONS=$MUTATIONS_DIR/thca-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $THCA_MUTATIONS

# # Endometrial
UCEC_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/ucec-mutations-wout-length-filtering.json
python process_mutations.py -mf $DATA_DIR/mafs/UCEC.maf \
       -hf $DATA_DIR/sample_lists/ucec-hypermutators.txt \
       -o $UCEC_MUTATIONS_WOUT_LENGTH_FILTERING

UCEC_MUTATIONS=$MUTATIONS_DIR/ucec-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $UCEC_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $UCEC_MUTATIONS

################################################################################
# GENERATE PERMUTED MATRICES AND COMPUTE WEIGHTS                               #
################################################################################
WEIGHTS_SUFFIX=`printf %.E $WEIGHTS_PERMUTATIONS`
PAIRS_PERMUTATION_SUFFIX=`printf %.E $PAIRS_PERMUTATIONS`
TOTAL_PERMUTATION_SUFFIX=`printf %.E $PAIRS_PERMUTATIONS`

# Set up directories
COADREAD_PERMUTATIONS=$PERMUTED_DIR/coadread
THCA_PERMUTATIONS=$PERMUTED_DIR/thca
UCEC_PERMUTATIONS=$PERMUTED_DIR/ucec
mkdir -p $COADREAD_PERMUTATIONS $THCA_PERMUTATIONS $UCEC_PERMUTATIONS

# Colorectal
COADREAD_WEIGHTS=$WEIGHTS_DIR/coadread-np${WEIGHTS_SUFFIX}.npy
python permute_matrix.py -mf $COADREAD_MUTATIONS -np $TOTAL_PERMUTATIONS -nc $NUM_CORES \
       -o $COADREAD_PERMUTATIONS

python compute_mutation_probabilities.py -pf $COADREAD_PERMUTATIONS/*.json \
       -mf $COADREAD_MUTATIONS -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES \
       -o $COADREAD_WEIGHTS

# Thyroid
THCA_WEIGHTS=$WEIGHTS_DIR/thca-np${WEIGHTS_SUFFIX}.npy
python permute_matrix.py -mf $THCA_MUTATIONS -np $TOTAL_PERMUTATIONS -nc $NUM_CORES \
       -o $THCA_PERMUTATIONS

python compute_mutation_probabilities.py -pf $THCA_PERMUTATIONS/*.json \
       -mf $THCA_MUTATIONS -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES \
       -o $THCA_WEIGHTS

# Endometrial
UCEC_WEIGHTS=$WEIGHTS_DIR/ucec-np${WEIGHTS_SUFFIX}.npy
python permute_matrix.py -mf $UCEC_MUTATIONS -np $TOTAL_PERMUTATIONS -nc $NUM_CORES \
       -o $UCEC_PERMUTATIONS

python compute_mutation_probabilities.py -pf $UCEC_PERMUTATIONS/*.json \
       -mf $UCEC_MUTATIONS -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES \
       -o $UCEC_WEIGHTS

################################################################################
# ENUMERATE PAIRS                                                              #
################################################################################
# Colorectal
COADREAD_PERMUTATIONAL_PAIRS=$PAIRS_DIR/coadread-pairs-permutational-np${PAIRS_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -k 2 -nc $NUM_CORES \
       -o $COADREAD_PERMUTATIONAL_PAIRS  -f $COADREAD_MIN_FREQ Permutational \
       -np $PAIRS_PERMUTATIONS -pf $COADREAD_PERMUTATIONS/*.json

COADREAD_WEIGHTED_EXACT_PAIRS=$PAIRS_DIR/coadread-pairs-weighted-exact-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 2\
       -f $COADREAD_MIN_FREQ -o $COADREAD_WEIGHTED_EXACT_PAIRS \
       Weighted -m Exact -wf $COADREAD_WEIGHTS

COADREAD_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/coadread-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 2\
       -f $COADREAD_MIN_FREQ -o $COADREAD_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $COADREAD_WEIGHTS

COADREAD_UNWEIGHTED_EXACT_PAIRS=$PAIRS_DIR/coadread-pairs-unweighted-exact.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 2\
       -f $COADREAD_MIN_FREQ -o $COADREAD_UNWEIGHTED_EXACT_PAIRS Unweighted -m Exact


# THYROID
THCA_PERMUTATIONAL_PAIRS=$PAIRS_DIR/thca-pairs-permutational-np${PAIRS_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -k 2 -nc $NUM_CORES \
       -o $THCA_PERMUTATIONAL_PAIRS  -f $THCA_MIN_FREQ Permutational \
       -np $PAIRS_PERMUTATIONS -pf $THCA_PERMUTATIONS/*.json

THCA_WEIGHTED_EXACT_PAIRS=$PAIRS_DIR/thca-pairs-weighted-exact-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 2\
       -f $THCA_MIN_FREQ -o $THCA_WEIGHTED_EXACT_PAIRS \
       Weighted -m Exact -wf $THCA_WEIGHTS

THCA_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/thca-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 2\
       -f $THCA_MIN_FREQ -o $THCA_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $THCA_WEIGHTS

THCA_UNWEIGHTED_EXACT_PAIRS=$PAIRS_DIR/thca-pairs-unweighted-exact.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 2\
       -f $THCA_MIN_FREQ -o $THCA_UNWEIGHTED_EXACT_PAIRS Unweighted -m Exact


# Endometrial
UCEC_PERMUTATIONAL_PAIRS=$PAIRS_DIR/ucec-pairs-permutational-np${PAIRS_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -k 2 -nc $NUM_CORES \
       -o $UCEC_PERMUTATIONAL_PAIRS  -f $UCEC_MIN_FREQ Permutational \
       -np $PAIRS_PERMUTATIONS -pf $UCEC_PERMUTATIONS/*.json

UCEC_WEIGHTED_EXACT_PAIRS=$PAIRS_DIR/ucec-pairs-weighted-exact-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 2\
       -f $UCEC_MIN_FREQ -o $UCEC_WEIGHTED_EXACT_PAIRS \
       Weighted -m Exact -wf $UCEC_WEIGHTS

UCEC_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/ucec-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 2\
       -f $UCEC_MIN_FREQ -o $UCEC_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $UCEC_WEIGHTS

UCEC_UNWEIGHTED_EXACT_PAIRS=$PAIRS_DIR/ucec-pairs-unweighted-exact.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 2\
       -f $UCEC_MIN_FREQ -o $UCEC_UNWEIGHTED_EXACT_PAIRS Unweighted -m Exact

################################################################################
# ENUMERATE TRIPLES                                                            #
################################################################################

# Colorectal
COADREAD_PERMUTATIONAL_TRIPLES=$TRIPLES_DIR/coadread-triples-permutational-np${TOTAL_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_PERMUTATIONAL_TRIPLES \
       Permutational -pf $COADREAD_PERMUTATIONS/*.json -np $TOTAL_PERMUTATIONS

COADREAD_UNWEIGHTED_EXACT_TRIPLES=$TRIPLES_DIR/coadread-triples-unweighted-exact.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_UNWEIGHTED_EXACT_TRIPLES \
       Unweighted -m Exact

COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $COADREAD_WEIGHTS

# Thyroid
THCA_PERMUTATIONAL_TRIPLES=$TRIPLES_DIR/thca-triples-permutational-np${TOTAL_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_PERMUTATIONAL_TRIPLES \
       Permutational -pf $THCA_PERMUTATIONS/*.json -np $TOTAL_PERMUTATIONS

THCA_UNWEIGHTED_EXACT_TRIPLES=$TRIPLES_DIR/thca-triples-unweighted-exact.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_UNWEIGHTED_EXACT_TRIPLES \
       Unweighted -m Exact

THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

THCA_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $THCA_WEIGHTS

# Endometrial
UCEC_PERMUTATIONAL_TRIPLES=$TRIPLES_DIR/ucec-triples-permutational-np${TOTAL_PERMUTATION_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_PERMUTATIONAL_TRIPLES \
       Permutational -pf $UCEC_PERMUTATIONS/*.json -np $TOTAL_PERMUTATIONS

UCEC_UNWEIGHTED_EXACT_TRIPLES=$TRIPLES_DIR/ucec-triples-unweighted-exact.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_UNWEIGHTED_EXACT_TRIPLES \
       Unweighted -m Exact

UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

UCEC_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $UCEC_WEIGHTS

################################################################################
# CREATE FIGURES AND TABLES                                                    #
################################################################################

# FIGURE 1: Distribution of number of mutations per sample
python $SCRIPTS_DIR/sample_mutation_frequency_plot.py -mf $THCA_MUTATIONS \
       $COADREAD_MUTATIONS $UCEC_MUTATIONS -c THCA COADREAD UCEC \
       -o $FIGURES_DIR/Figure1.pdf

# FIGURE 2: Mutation probability matrices
python $SCRIPTS_DIR/weights_matrix.py -mf $THCA_MUTATIONS $COADREAD_MUTATIONS \
       $UCEC_MUTATIONS -wf $THCA_WEIGHTS $COADREAD_WEIGHTS $UCEC_WEIGHTS \
       -c THCA COADREAD UCEC -o $FIGURES_DIR/Figure2.pdf

# FIGURE 3 and TABLE 1: Comparison of exact and saddlepoint for unweighted test
python $SCRIPTS_DIR/unweighted_comparison.py -c THCA COADREAD UCEC \
       -sf $THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       $COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       $UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES\
       -ef $THCA_UNWEIGHTED_EXACT_TRIPLES $COADREAD_UNWEIGHTED_EXACT_TRIPLES \
       $UCEC_UNWEIGHTED_EXACT_TRIPLES\
       -ff $FIGURES_DIR/Figure3.png -tf $TABLES_DIR/Table1.tsv

# FIGURE 4 and TABLE 3: Pairs
python $SCRIPTS_DIR/pairs_summary.py -c THCA THCA THCA THCA COADREAD \
       COADREAD COADREAD COADREAD UCEC UCEC UCEC UCEC -np $PAIRS_PERMUTATIONS \
       -pf $THCA_PERMUTATIONAL_PAIRS $THCA_WEIGHTED_EXACT_PAIRS \
       $THCA_WEIGHTED_SADDLEPOINT_PAIRS $THCA_UNWEIGHTED_EXACT_PAIRS \
       $COADREAD_PERMUTATIONAL_PAIRS $COADREAD_WEIGHTED_EXACT_PAIRS \
       $COADREAD_WEIGHTED_SADDLEPOINT_PAIRS $COADREAD_UNWEIGHTED_EXACT_PAIRS \
       $UCEC_PERMUTATIONAL_PAIRS $UCEC_WEIGHTED_EXACT_PAIRS \
       $UCEC_WEIGHTED_SADDLEPOINT_PAIRS $UCEC_UNWEIGHTED_EXACT_PAIRS \
       -ff $FIGURES_DIR/Figure4.pdf -tf $TABLES_DIR/Table2.tsv $TABLES_DIR/Table3.tsv

# FIGURE 5: Unweighted vs. Permutational (N=10^6)and Weighted (N=10^3) vs.
# Permutational

python $SCRIPTS_DIR/triple_pval_scatter.py -np $TOTAL_PERMUTATIONS \
       -o $FIGURES_DIR/Figure5.png -uwf $COADREAD_UNWEIGHTED_EXACT_TRIPLES \
       $THCA_UNWEIGHTED_EXACT_TRIPLES $UCEC_UNWEIGHTED_EXACT_TRIPLES \
       -wf $COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES \
       $THCA_WEIGHTED_SADDLEPOINT_TRIPLES $UCEC_WEIGHTED_SADDLEPOINT_TRIPLES \
       -pf $COADREAD_PERMUTATIONAL_TRIPLES $THCA_PERMUTATIONAL_TRIPLES \
       $UCEC_PERMUTATIONAL_TRIPLES

# TABLE 4: THCA results
echo "THCA"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS -wf $THCA_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $THCA_UNWEIGHTED_EXACT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table4

# TABLE 5: COADREAD results
echo "COADREAD"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $COADREAD_MUTATIONS -wf $COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $COADREAD_UNWEIGHTED_EXACT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table5

# TABLE 6: UCEC results
echo "UCEC"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $UCEC_MUTATIONS -wf $UCEC_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $UCEC_UNWEIGHTED_EXACT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table6
