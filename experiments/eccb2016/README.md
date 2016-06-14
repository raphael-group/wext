# EXPERIMENTS FOR ECCB/BIOINFORMATICS 2016 PUBLICATION #

### DESCRIPTION ###

This directory contains the scripts for reproducing the experiments, tables, and figures for the 2016 ECCB/Bioinformatics publication. You can see the commands for reproducing the in `commands.sh`. This will perform the following:

1. Download and unpack the raw mutation data used in the experiments into `data/`.
2. Preprocess the mutation data and generate permuted matrices into `data/mutations`, `data/permuted`, and `data/weights`.
3. Reproduce the tables and figures from the paper to `tables/` and `figures/`, with intermediate output stored in `output`.

Note that in the paper we restricted our analysis to genes with lengths, which was not mentioned in the paper.

### SETUP ###

Our experiments require additional Python modules to be installed. You can install them with:

    pip install -r requirements.txt

### SETTINGS ###

We include several settings that can be changed in `commands.sh`. Currently, they are set to reproduce the results from the ECCB/Bioinformatics 2016 paper. They will also take advantage of a GridEngine (`GRID_ENGINE` and `NUM_GRID_NODES=400`) installation and to use multiple cores (`NUM_CORES`) on the machine on which `commands.sh` is run.

#### GRID ENGINE ####
We include several commands that will take advantage of a [GridEngine](http://gridscheduler.sourceforge.net/) cluster if available. To use GridEngine, set the bash variable `GRID_ENGINE=true`.

### DATA ###

In order to reproduce the experiments, you will need to download the COADREAD, THCA, and UCEC MAFs and the list of hypermutator samples in the tarball we've posted on our website. We've included the relevant commands in `commands.sh`, but you can also find the data at the URL below:

> http://compbio-research.cs.brown.edu/projects/weighted-exclusivity-test/data/eccb2016.tar

### REFERENCE ###

M.D.M. Leiserson, M.A. Reyna, B.J. Raphael. A Weighted Exact Test for Mutually Exclusive Mutations in Cancer. ECCB/Bioinformatics (2016). To appear.
