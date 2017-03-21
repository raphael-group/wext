###############################################################################
# SETTINGS, FILES, AND DIRECTORIES
###############################################################################

# Settings
ifndef K
	K=2
endif

ifndef NP
	NP=100
endif

ifndef RANDOM_SEED
	RANDOM_SEED=2
endif

ifndef PANCAN_NUM_SAMPLES
	PANCAN_NUM_SAMPLES=100
endif

ifndef PANCAN_NUM_GENES
	PANCAN_NUM_GENES=50
endif

ifndef PANCAN_BMR
	PANCAN_BMR=0.01
endif

# Directories
REPO_DIR=$(realpath ..)
EXAMPLES_DIR=$(realpath .)

PANCAN_DIR=$(EXAMPLES_DIR)/pancan
PANCAN_DATA_DIR=$(PANCAN_DIR)/data
PANCAN_OUTPUT_DIR=$(PANCAN_DIR)/output

OUTPUT_DIRS=$(PANCAN_OUTPUT_DIR)
DATA_DIRS=$(PANCAN_DATA_DIR)

# Data files
D1_ABERRATIONS=$(PANCAN_DATA_DIR)/dataset1-aberrations.tsv
D2_ABERRATIONS=$(PANCAN_DATA_DIR)/dataset2-aberrations.tsv

# Output files
PANCAN_MUTATIONS=$(PANCAN_DATA_DIR)/pancan-dataset-mutations.json
D1_MUTATIONS=$(PANCAN_DATA_DIR)/dataset1-mutations.json
D2_MUTATIONS=$(PANCAN_DATA_DIR)/dataset2-mutations.json

PANCAN_WEIGHTS=$(PANCAN_DATA_DIR)/pancan-dataset-weights.npy
D1_WEIGHTS=$(PANCAN_DATA_DIR)/dataset1-weights.npy
D2_WEIGHTS=$(PANCAN_DATA_DIR)/dataset2-weights.npy

PANCAN_ACROSS_SETS=$(PANCAN_OUTPUT_DIR)/pancan-across-sets
PANCAN_BETWEEN_SETS=$(PANCAN_OUTPUT_DIR)/pancan-between-sets

# Scripts
PROCESS_MUTATIONS=$(REPO_DIR)/process_mutations.py
COMPUTE_MUT_PROBS=$(REPO_DIR)/compute_mutation_probabilities.py
FIND_EXCLUSIVE_SETS=$(REPO_DIR)/find_exclusive_sets.py

###############################################################################
# GENERAL COMMANDS (MAKE ALL EXAMPLES)
###############################################################################
all: pancan

setup:
	mkdir -p $(OUTPUT_DIRS) $(DATA_DIRS)

clean:
	rm -r $(PANCAN_DIR)

###############################################################################
# PAN-CANCER COMMANDS
###############################################################################
pancan: $(PANCAN_ACROSS_SETS)-k$K.tsv $(PANCAN_BETWEEN_SETS)-k$K.tsv

# Generate the Pan-Cancer datasets
pancan-aberrations: $(PANCAN_DATA_DIR)/dataset%-aberrations.tsv
$(PANCAN_DATA_DIR)/dataset%-aberrations.tsv: | setup
	python generate_data.py \
		-rs $(RANDOM_SEED) \
		-o $(PANCAN_DATA_DIR)/dataset \
		pancan \
		-ns $(PANCAN_NUM_SAMPLES) \
		-ng $(PANCAN_NUM_GENES) \
		-b $(PANCAN_BMR)

# Process the mutation events files into separate and combined mutation datasets
pancan-mutations: $(PANCAN_MUTATIONS) $(D1_MUTATIONS) $(D2_MUTATIONS)
$(PANCAN_MUTATIONS): pancan-aberrations
	python $(PROCESS_MUTATIONS) \
		-m $(D1_ABERRATIONS) \
		-m $(D2_ABERRATIONS) \
		-ct Dataset1 Dataset2 \
		-o $(PANCAN_MUTATIONS)

$(D1_MUTATIONS): pancan-aberrations
	python $(PROCESS_MUTATIONS) \
		-m $(D1_ABERRATIONS) \
		-ct Dataset1 \
		-o $(D1_MUTATIONS)

$(D2_MUTATIONS): pancan-aberrations
	python $(PROCESS_MUTATIONS) \
		-m $(D2_ABERRATIONS) \
		-ct Dataset2 \
		-o $(D2_MUTATIONS)

# Compute mutation probabilities per cancer type independently _and_ together
pancan-weights: $(PANCAN_WEIGHTS) $(D1_WEIGHTS) $(D2_WEIGHTS)
$(PANCAN_WEIGHTS): $(PANCAN_MUTATIONS)
	python $(COMPUTE_MUT_PROBS) \
		-mf $(PANCAN_MUTATIONS) \
		-wf $(PANCAN_WEIGHTS) \
		-np $(NP)

$(D1_WEIGHTS): $(D1_MUTATIONS)
	python $(COMPUTE_MUT_PROBS) \
		-mf $(D1_MUTATIONS) \
		-wf $(D1_WEIGHTS) \
		-np $(NP)

$(D2_WEIGHTS): $(D2_MUTATIONS)
	python $(COMPUTE_MUT_PROBS) \
		-mf $(D2_MUTATIONS) \
		-wf $(D2_WEIGHTS) \
		-np $(NP)

# Find exclusive sets _between_ cancer types (possible subtype exclusivity; BETWEEN_ME)
$(PANCAN_BETWEEN_SETS)-k$K.tsv: $(PANCAN_WEIGHTS)
	python $(FIND_EXCLUSIVE_SETS) \
		-o $(PANCAN_BETWEEN_SETS) \
		-mf $(PANCAN_MUTATIONS) \
		-ks $(K) \
		-s Enumerate \
		WRE \
		-m Saddlepoint \
		-wf $(PANCAN_WEIGHTS)

# Find exclusive sets _across_ cancer types (no subtype exclusivity; ACROSS_ME)
$(PANCAN_ACROSS_SETS)-k$K.tsv: $(D1_WEIGHTS) $(D2_WEIGHTS)
	python $(FIND_EXCLUSIVE_SETS) \
		-o $(PANCAN_ACROSS_SETS) \
		-mf $(D1_MUTATIONS) $(D2_MUTATIONS) \
		-ks $(K) \
		-s Enumerate \
		WRE \
		-m Saddlepoint \
		-wf $(D1_WEIGHTS) $(D2_WEIGHTS)
