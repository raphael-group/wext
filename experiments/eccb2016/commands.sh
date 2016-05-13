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

# Thyroid
THCA_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/thca-mutations-wout-length-filtering.json
python process_mutations.py -mf $DATA_DIR/mafs/THCA.maf \
       -o $THCA_MUTATIONS_WOUT_LENGTH_FILTERING

THCA_MUTATIONS=$MUTATIONS_DIR/thca-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $THCA_MUTATIONS

# Endometrial
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
python compute_mutation_probabilities.py \
    -mf $COADREAD_MUTATIONS \
    -wf $COADREAD_WEIGHTS \
    -np $TOTAL_PERMUTATIONS \
    -nc $NUM_CORES

# Thyroid
THCA_WEIGHTS=$WEIGHTS_DIR/thca-np${WEIGHTS_SUFFIX}.npy
python compute_mutation_probabilities.py \
    -mf $THCA_MUTATIONS \
    -wf $THCA_WEIGHTS \
    -np $TOTAL_PERMUTATIONS \
    -nc $NUM_CORES

# Endometrial
UCEC_WEIGHTS=$WEIGHTS_DIR/ucec-np${WEIGHTS_SUFFIX}.npy
python compute_mutation_probabilities.py \
    -mf $UCEC_MUTATIONS \
    -wf $UCEC_WEIGHTS \
    -np $TOTAL_PERMUTATIONS \
    -nc $NUM_CORES

################################################################################
# ENUMERATE PAIRS                                                              #
################################################################################

# Colorectal
COADREAD_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/coadread-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 2\
       -f $COADREAD_MIN_FREQ -o $COADREAD_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $COADREAD_WEIGHTS


# THYROID
THCA_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/thca-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 2\
       -f $THCA_MIN_FREQ -o $THCA_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $THCA_WEIGHTS


# Endometrial
UCEC_WEIGHTED_SADDLEPOINT_PAIRS=$PAIRS_DIR/ucec-pairs-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 2\
       -f $UCEC_MIN_FREQ -o $UCEC_WEIGHTED_SADDLEPOINT_PAIRS \
       Weighted -m Saddlepoint -wf $UCEC_WEIGHTS

################################################################################
# ENUMERATE TRIPLES                                                            #
################################################################################

# Colorectal
COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $COADREAD_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $COADREAD_WEIGHTS

# Thyroid
THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

THCA_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $THCA_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $THCA_WEIGHTS

# Endometrial
UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-triples-unweighted-saddlepoint.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES \
       Unweighted -m Saddlepoint

UCEC_WEIGHTED_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-triples-weighted-saddlepoint-nw${WEIGHTS_SUFFIX}.json
python compute_exclusivity.py -mf $UCEC_MUTATIONS -nc $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_WEIGHTED_SADDLEPOINT_TRIPLES \
       Weighted -m Saddlepoint -wf $UCEC_WEIGHTS

################################################################################
# SUMMARIZE RESULTS                                                            #
################################################################################

# THCA results
echo "THCA"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS -wf $THCA_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $THCA_UNWEIGHTED_SADDLEPOINT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table4

# COADREAD results
echo "COADREAD"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $COADREAD_MUTATIONS -wf $COADREAD_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $COADREAD_UNWEIGHTED_SADDLEPOINT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table5

# UCEC results
echo "UCEC"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $UCEC_MUTATIONS -wf $UCEC_WEIGHTED_SADDLEPOINT_TRIPLES \
       -uf $UCEC_UNWEIGHTED_SADDLEPOINT_TRIPLES -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table6
