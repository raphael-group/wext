# Tests #

### Description ###

This directory contains the scripts for running experiments on our ECCB 2016 submission data. You can see the commands in `commands.sh`. This will perform the following:

1. Download and unpack raw mutation data used in the experiments into `data/`.
2. Process mutation data into `data/mutations` and `data/weights`.
3. Produce tables in `tables/` with intermediate output stored in `output`.

### Set up ###

Our experiments require additional Python modules to be installed. You can install them with:

    pip install -r requirements.txt

### Data ###

In order to run these experiments, you will need to download the COADREAD, THCA, and UCEC MAFs and the list of hypermutator samples in the tarball we've posted on our website. We've included the relevant commands in `commands.sh`, but you can also find the data at the URL below:

> http://compbio-research.cs.brown.edu/projects/weighted-exclusivity-test/data/eccb2016.tar
