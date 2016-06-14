#!/bin/sh

################################################################################
# SETTINGS                                                                     #
################################################################################

# Parameters for analysis
COADREAD_MIN_FREQ=20
THCA_MIN_FREQ=5
UCEC_MIN_FREQ=30
LENGTH_THRESHOLD=600
FDR_CUTOFF=0.001
EXTRA_PERMUTATIONS=1000000
PAIRS_PERMUTATIONS=10000
WEIGHTS_PERMUTATIONS=1000

# Parallel processing
NUM_CORES=25
GRID_ENGINE=true
NUM_GRID_NODES=400
COADREAD_GRID_BATCH_SIZE=250

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
WEXT_DIR=`pwd`

# Colorectal
COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/coadread-mutations-wout-length-filtering.json
python process_mutations.py -m $DATA_DIR/mafs/COADREAD.maf -ct COADREAD \
       -hf $DATA_DIR/sample_lists/coadread-hypermutators.txt \
       -o $COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING

COADREAD_MUTATIONS=$MUTATIONS_DIR/coadread-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $COADREAD_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $COADREAD_MUTATIONS

# # Thyroid
THCA_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/thca-mutations-wout-length-filtering.json
python process_mutations.py -m $DATA_DIR/mafs/THCA.maf -ct THCA \
       -o $THCA_MUTATIONS_WOUT_LENGTH_FILTERING

THCA_MUTATIONS=$MUTATIONS_DIR/thca-mutations.json
python $SCRIPTS_DIR/remove_genes_with_no_length.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS_WOUT_LENGTH_FILTERING \
       -o $THCA_MUTATIONS

# # Endometrial
UCEC_MUTATIONS_WOUT_LENGTH_FILTERING=$MUTATIONS_DIR/ucec-mutations-wout-length-filtering.json
python process_mutations.py -m $DATA_DIR/mafs/UCEC.maf -ct UCEC \
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
EXTRA_PERMUTATION_SUFFIX=`printf %.E $EXTRA_PERMUTATIONS`

# Set up directories
COADREAD_PAIRS_PERMUTATIONS=$PERMUTED_DIR/coadread-$PAIRS_PERMUTATION_SUFFIX
COADREAD_EXTRA_PERMUTATIONS=$PERMUTED_DIR/coadread-$EXTRA_PERMUTATION_SUFFIX
THCA_PAIRS_PERMUTATIONS=$PERMUTED_DIR/thca-$PAIRS_PERMUTATION_SUFFIX
UCEC_PAIRS_PERMUTATIONS=$PERMUTED_DIR/ucec-$PAIRS_PERMUTATION_SUFFIX
mkdir -p $COADREAD_PAIRS_PERMUTATIONS $THCA_PAIRS_PERMUTATIONS $UCEC_PAIRS_PERMUTATIONS $COADREAD_EXTRA_PERMUTATIONS

# Colorectal
COADREAD_WEIGHTS=$WEIGHTS_DIR/coadread-np${WEIGHTS_SUFFIX}.npy
python compute_mutation_probabilities.py -mf $COADREAD_MUTATIONS -v 1 \
       -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES -wf $COADREAD_WEIGHTS

if [ "$GRID_ENGINE" = true ]; then
    qsub -t 1-$PAIRS_PERMUTATIONS -tc $NUM_GRID_NODES -sync y -V -cwd -N per \
         -o $COADREAD_PAIRS_PERMUTATIONS/command.out \
         -e $COADREAD_PAIRS_PERMUTATIONS/command.err \
         $SCRIPTS_DIR/permute_single_matrix.py -mf $COADREAD_MUTATIONS \
         -o $COADREAD_PAIRS_PERMUTATIONS/coadread-permuted -wd $WEXT_DIR
    
    qsub -t 1-$EXTRA_PERMUTATIONS -tc $NUM_GRID_NODES -sync y -V -cwd -N per \
         -o $COADREAD_EXTRA_PERMUTATIONS/command.out \
         -e $COADREAD_EXTRA_PERMUTATIONS/command.err \
         $SCRIPTS_DIR/permute_single_matrix.py -mf $COADREAD_MUTATIONS \
         -o $COADREAD_EXTRA_PERMUTATIONS/coadread-permuted -wd $WEXT_DIR
else
    python compute_mutation_probabilities.py -mf $COADREAD_MUTATIONS -v 1 \
           -np $PAIRS_PERMUTATIONS -nc $NUM_CORES -pd $COADREAD_PAIRS_PERMUTATIONS
    python compute_mutation_probabilities.py -mf $COADREAD_MUTATIONS -v 1 \
           -np $EXTRA_PERMUTATIONS -nc $NUM_CORES -pd $COADREAD_EXTRA_PERMUTATIONS
fi

# Thyroid
THCA_WEIGHTS=$WEIGHTS_DIR/thca-np${WEIGHTS_SUFFIX}.npy
python compute_mutation_probabilities.py -mf $THCA_MUTATIONS -v 1 \
       -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES -wf $THCA_WEIGHTS

if [ "$GRID_ENGINE" = true ]; then
    qsub -t 1-$PAIRS_PERMUTATIONS -tc $NUM_GRID_NODES -sync y -V -cwd -N per \
         -o $THCA_PAIRS_PERMUTATIONS/command.out \
         -e $THCA_PAIRS_PERMUTATIONS/command.err \
         $SCRIPTS_DIR/permute_single_matrix.py -mf $THCA_MUTATIONS \
         -o $THCA_PAIRS_PERMUTATIONS/thca-permuted -wd $WEXT_DIR
else
    python compute_mutation_probabilities.py -mf $THCA_MUTATIONS -v 1 \
           -np $PAIRS_PERMUTATIONS -nc $NUM_CORES -pd $THCA_PAIRS_PERMUTATIONS
fi


# Endometrial
UCEC_WEIGHTS=$WEIGHTS_DIR/ucec-np${WEIGHTS_SUFFIX}.npy
python compute_mutation_probabilities.py -mf $UCEC_MUTATIONS -v 1 \
       -np $WEIGHTS_PERMUTATIONS -nc $NUM_CORES -wf $UCEC_WEIGHTS

if [ "$GRID_ENGINE" = true ]; then
    qsub -t 1-$PAIRS_PERMUTATIONS -tc $NUM_GRID_NODES -sync y -V -cwd -N per \
         -o $UCEC_PAIRS_PERMUTATIONS/command.out \
         -e $UCEC_PAIRS_PERMUTATIONS/command.err \
         $SCRIPTS_DIR/permute_single_matrix.py -mf $UCEC_MUTATIONS \
         -o $UCEC_PAIRS_PERMUTATIONS/ucec-permuted -wd $WEXT_DIR
else
    python compute_mutation_probabilities.py -mf $UCEC_MUTATIONS -v 1 \
           -np $PAIRS_PERMUTATIONS -nc $NUM_CORES -pd $UCEC_PAIRS_PERMUTATIONS/

fi

################################################################################
# ENUMERATE PAIRS                                                              #
################################################################################
# Colorectal
COADREAD_RCE_PERMUTATIONS_PAIRS=$PAIRS_DIR/coadread-RCE-permutations-np${PAIRS_PERMUTATION_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 2 \
       -c $NUM_CORES -o $COADREAD_RCE_PERMUTATIONS_PAIRS \
       -f $COADREAD_MIN_FREQ --json_format \
       RCE -np $PAIRS_PERMUTATIONS -pd $COADREAD_PAIRS_PERMUTATIONS/

COADREAD_WRE_EXACT_PAIRS=$PAIRS_DIR/coadread-WRE-exact-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 2 \
       -c $NUM_CORES -f $COADREAD_MIN_FREQ -o $COADREAD_WRE_EXACT_PAIRS \
       --json_format WRE -m Exact -wf $COADREAD_WEIGHTS

COADREAD_WRE_SADDLEPOINT_PAIRS=$PAIRS_DIR/coadread-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 2 \
        -c $NUM_CORES -f $COADREAD_MIN_FREQ -o $COADREAD_WRE_SADDLEPOINT_PAIRS\
       --json_format WRE -m Saddlepoint -wf $COADREAD_WEIGHTS

COADREAD_RE_EXACT_PAIRS=$PAIRS_DIR/coadread-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 2 \
       -c $NUM_CORES -f $COADREAD_MIN_FREQ -o $COADREAD_RE_EXACT_PAIRS \
       --json_format RE -m Exact

# THYROID
THCA_RCE_PERMUTATIONS_PAIRS=$PAIRS_DIR/thca-RCE-permutations-np${PAIRS_PERMUTATION_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -k 2 \
       -c $NUM_CORES -o $THCA_RCE_PERMUTATIONS_PAIRS  -f $THCA_MIN_FREQ \
       --json_format RCE -np $PAIRS_PERMUTATIONS -pd $THCA_PAIRS_PERMUTATIONS/

THCA_WRE_EXACT_PAIRS=$PAIRS_DIR/thca-WRE-exact-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES \
       -k 2 -f $THCA_MIN_FREQ -o $THCA_WRE_EXACT_PAIRS --json_format \
       WRE -m Exact -wf $THCA_WEIGHTS

THCA_WRE_SADDLEPOINT_PAIRS=$PAIRS_DIR/thca-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES -k 2 \
       -f $THCA_MIN_FREQ -o $THCA_WRE_SADDLEPOINT_PAIRS --json_format \
       WRE -m Saddlepoint -wf $THCA_WEIGHTS

THCA_RE_EXACT_PAIRS=$PAIRS_DIR/thca-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES \
       -k 2 -f $THCA_MIN_FREQ -o $THCA_RE_EXACT_PAIRS --json_format RE -m Exact

# Endometrial
UCEC_RCE_PERMUTATIONS_PAIRS=$PAIRS_DIR/ucec-RCE-permutations-np${PAIRS_PERMUTATION_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES \
       -k 2 -o $UCEC_RCE_PERMUTATIONS_PAIRS  -f $UCEC_MIN_FREQ --json_format \
       RCE -np $PAIRS_PERMUTATIONS -pd $UCEC_PAIRS_PERMUTATIONS/

UCEC_WRE_EXACT_PAIRS=$PAIRS_DIR/ucec-WRE-exact-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES \
       -k 2 -f $UCEC_MIN_FREQ -o $UCEC_WRE_EXACT_PAIRS --json_format \
       WRE -m Exact -wf $UCEC_WEIGHTS

UCEC_WRE_SADDLEPOINT_PAIRS=$PAIRS_DIR/ucec-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES \
       -k 2 -f $UCEC_MIN_FREQ -o $UCEC_WRE_SADDLEPOINT_PAIRS --json_format \
       WRE -m Saddlepoint -wf $UCEC_WEIGHTS

UCEC_RE_EXACT_PAIRS=$PAIRS_DIR/ucec-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES \
       -k 2 -f $UCEC_MIN_FREQ -o $UCEC_RE_EXACT_PAIRS --json_format \
       RE -m Exact

################################################################################
# ENUMERATE TRIPLES                                                            #
################################################################################

# Colorectal
COADREAD_RCE_PERMUTATIONS_TRIPLES=$TRIPLES_DIR/coadread-RCE-permutations-np${EXTRA_PERMUTATION_SUFFIX}
if [ "$GRID_ENGINE" = true ]; then
    COADREAD_RCE_PERMUTATIONS_GRID=$OUTPUT_DIR/grid/coadread-np${EXTRA_PERMUTATION_SUFFIX}
    mkdir -p $COADREAD_RCE_PERMUTATIONS_GRID
    
    qsub -t 1-$(($EXTRA_PERMUTATIONS / $COADREAD_GRID_BATCH_SIZE)) \
         -tc $NUM_GRID_NODES -sync y -V -cwd -N per \
         -o $COADREAD_RCE_PERMUTATIONS_GRID/command.out \
         -e $COADREAD_RCE_PERMUTATIONS_GRID/command.err \
         $SCRIPTS_DIR/permutation_test_helper.py -m $COADREAD_MUTATIONS -k 3 \
         -i $COADREAD_EXTRA_PERMUTATIONS/ -b $COADREAD_GRID_BATCH_SIZE -np $EXTRA_PERMUTATIONS \
         -f $COADREAD_MIN_FREQ -o $COADREAD_RCE_PERMUTATIONS_GRID/coadread-permutation-test \
         -w $WEXT_DIR

    python $SCRIPTS_DIR/reconcile_grid_permutation_test.py \
           -i $COADREAD_RCE_PERMUTATIONS_GRID -b $COADREAD_GRID_BATCH_SIZE \
           -o $COADREAD_RCE_PERMUTATIONS_TRIPLES -wd $WEXT_DIR -c $NUM_CORES

else    
    python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 3 -v 4 \
           -c $NUM_CORES -f $COADREAD_MIN_FREQ -o $COADREAD_RCE_PERMUTATIONS_TRIPLES \
           --json_format RCE -pd $COADREAD_EXTRA_PERMUTATIONS/ -np $EXTRA_PERMUTATIONS
fi

COADREAD_RE_EXACT_TRIPLES=$TRIPLES_DIR/coadread-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -k 3 \
       -c $NUM_CORES -f $COADREAD_MIN_FREQ -o $COADREAD_RE_EXACT_TRIPLES \
       --json_format RE -m Exact

COADREAD_RE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-RE-saddlepoint
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -c $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_RE_SADDLEPOINT_TRIPLES \
       --json_format RE -m Saddlepoint

COADREAD_WRE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/coadread-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $COADREAD_MUTATIONS -c $NUM_CORES -k 3 \
       -f $COADREAD_MIN_FREQ -o $COADREAD_WRE_SADDLEPOINT_TRIPLES \
       --json_format WRE -m Saddlepoint -wf $COADREAD_WEIGHTS

# Thyroid
THCA_RE_EXACT_TRIPLES=$TRIPLES_DIR/thca-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_RE_EXACT_TRIPLES --json_format \
       RE -m Exact

THCA_RE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-RE-saddlepoint
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_RE_SADDLEPOINT_TRIPLES --json_format \
       RE -m Saddlepoint

THCA_WRE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/thca-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $THCA_MUTATIONS -c $NUM_CORES -k 3 \
       -f $THCA_MIN_FREQ -o $THCA_WRE_SADDLEPOINT_TRIPLES --json_format \
       WRE -m Saddlepoint -wf $THCA_WEIGHTS

# Endometrial
UCEC_RE_EXACT_TRIPLES=$TRIPLES_DIR/ucec-RE-exact
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_RE_EXACT_TRIPLES --json_format \
       RE -m Exact

UCEC_RE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-RE-saddlepoint
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_RE_SADDLEPOINT_TRIPLES --json_format \
       RE -m Saddlepoint

UCEC_WRE_SADDLEPOINT_TRIPLES=$TRIPLES_DIR/ucec-WRE-saddlepoint-nw${WEIGHTS_SUFFIX}
python find_exclusive_sets.py -s Enumerate -mf $UCEC_MUTATIONS -c $NUM_CORES -k 3 \
       -f $UCEC_MIN_FREQ -o $UCEC_WRE_SADDLEPOINT_TRIPLES --json_format \
       WRE -m Saddlepoint -wf $UCEC_WEIGHTS

################################################################################
# CREATE FIGURES AND TABLES                                                    #
################################################################################

# FIGURE S1: Distribution of number of mutations per sample
python $SCRIPTS_DIR/sample_mutation_frequency_plot.py -mf $THCA_MUTATIONS \
       $COADREAD_MUTATIONS $UCEC_MUTATIONS -c THCA COADREAD UCEC \
       -o $FIGURES_DIR/FigureS1.pdf

# FIGURE 1: Mutation probability matrices
python $SCRIPTS_DIR/weights_matrix.py -mf $THCA_MUTATIONS $COADREAD_MUTATIONS \
       $UCEC_MUTATIONS -wf $THCA_WEIGHTS $COADREAD_WEIGHTS $UCEC_WEIGHTS \
       -c THCA COADREAD UCEC -o $FIGURES_DIR/Figure1.pdf

# FIGURE S2: Comparison of exact and saddlepoint for R-Exclusivity test
python $SCRIPTS_DIR/unweighted_comparison.py -c THCA COADREAD UCEC \
       -sf ${THCA_RE_SADDLEPOINT_TRIPLES}-k3.json \
       ${COADREAD_RE_SADDLEPOINT_TRIPLES}-k3.json \
       ${UCEC_RE_SADDLEPOINT_TRIPLES}-k3.json \
       -ef $THCA_RE_EXACT_TRIPLES-k3.json \
       ${COADREAD_RE_EXACT_TRIPLES}-k3.json \
       ${UCEC_RE_EXACT_TRIPLES}-k3.json \
       -ff $FIGURES_DIR/Figure3.pdf

# FIGURE 2 and TABLE 2: Pairs
python $SCRIPTS_DIR/pairs_summary.py -c THCA THCA THCA THCA COADREAD \
       COADREAD COADREAD COADREAD UCEC UCEC UCEC UCEC -np $PAIRS_PERMUTATIONS \
       -pf ${THCA_RCE_PERMUTATIONS_PAIRS}-k2.json \
           ${THCA_WRE_EXACT_PAIRS}-k2.json\
           ${THCA_WRE_SADDLEPOINT_PAIRS}-k2.json \
           ${THCA_RE_EXACT_PAIRS}-k2.json \
           ${COADREAD_RCE_PERMUTATIONS_PAIRS}-k2.json \
           ${COADREAD_WRE_EXACT_PAIRS}-k2.json \
           ${COADREAD_WRE_SADDLEPOINT_PAIRS}-k2.json \
           ${COADREAD_RE_EXACT_PAIRS}-k2.json \
           ${UCEC_RCE_PERMUTATIONS_PAIRS}-k2.json \
           ${UCEC_WRE_EXACT_PAIRS}-k2.json \
           ${UCEC_WRE_SADDLEPOINT_PAIRS}-k2.json \
           ${UCEC_RE_EXACT_PAIRS}-k2.json \
       -ff $FIGURES_DIR/Figure2.pdf -tf $TABLES_DIR/Table2.tsv

# FIGURE 3: Unweighted vs. Permutational (N=10^6) and Weighted (N=10^3) vs.
# Permutational
python $SCRIPTS_DIR/triple_pval_scatter.py -np $EXTRA_PERMUTATIONS \
       -o   $FIGURES_DIR/Figure3.pdf \
       -uwf ${COADREAD_RE_EXACT_TRIPLES}-k3.json \
       -wf  ${COADREAD_WRE_SADDLEPOINT_TRIPLES}-k3.json \
       -pf  ${COADREAD_RCE_PERMUTATIONS_TRIPLES}-k3.json 

# TABLE S1: THCA results
echo "THCA"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $THCA_MUTATIONS -wf ${THCA_WRE_SADDLEPOINT_TRIPLES}-k3.json \
       -uf ${THCA_RE_EXACT_TRIPLES}-k3.json -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/TableS1

# TABLE 3: COADREAD results
echo "COADREAD"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $COADREAD_MUTATIONS -wf ${COADREAD_WRE_SADDLEPOINT_TRIPLES}-k3.json \
       -uf ${COADREAD_RE_EXACT_TRIPLES}-k3.json -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table3

# TABLE 4: UCEC results
echo "UCEC"
python $SCRIPTS_DIR/results_table.py -lf $GENE_LENGTH_FILE \
       -mf $UCEC_MUTATIONS -wf ${UCEC_WRE_SADDLEPOINT_TRIPLES}-k3.json \
       -uf ${UCEC_RE_EXACT_TRIPLES}-k3.json -nt 5 -lt $LENGTH_THRESHOLD \
       -f $FDR_CUTOFF -o $TABLES_DIR/Table4
