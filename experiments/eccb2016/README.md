# Experiments for ECCB 2016 submission #

### Description ###

This directory contains the scripts for reproducing the experiments, tables, and figures for the 2016 ECCB submission. You can see the commands for reproducing the in `commands.sh`. This will perform the following:

1. Download and unpack the raw mutation data used in the experiments into `data/`.
2. Preprocess the mutation data and generate permuted matrices into `data/mutations`, `data/permuted`, and `data/weights`.
3. Reproduce the tables and figures from the paper to `tables/` and `figures/`, with intermediate output stored in `output`.

### Set up ###

Our experiments require additional Python modules to be installed. You can install them with:

    pip install -r requirements.txt

### Data ###

In order to reproduce the experiments, you will need to download the COADREAD, THCA, and UCEC MAFs and the list of hypermutator samples in the tarball we've posted on our website. We've included the relevant commands in `commands.sh`, but you can also find the data at the URL below:

> http://compbio-research.cs.brown.edu/projects/weighted-exclusivity-test/data/eccb2016.tar

### Differences ###

There are only a few, relatively small differences in the experiments here and those in the actual submission:

* The mutation probability matrices in `figures/Figure 2.pdf` will include the "pseudocounts" of 1/2N (where N is the number of permutations) for cells with zeros.  These were missing from the ECCB 2016 submission and were shown as zeros.
* Obviously, the permutations are randomly generated, so the weighted and permutational results may look different, though we wouldn't expect this difference to be large. For example, when I ran this code to reproduce the experiments, we found 49, 5305, and 6804 triples with FDR < 0.001 for THCA, COADREAD, and UCEC, respectively (compared to 48, 5286, 6790 in the paper). Related, we also found that the top triples (ranked 4-10) in THCA all had very similar weighted p-values so were reorded when regenerating the weight matrices.
* The `commands.sh` script currently only generates 10,000 permutations for comparing the permutational vs. the saddlepoint (see Figure 5 in the submission), while we used 1,000,000 in the paper. This is because 1,000,000 permutations is computationally intensive (it will take days to generate without parallelization) and storage intensive (terabytes). Thus, the correlations between the weighted test and permutational reported in Section 3.4 and Figure 5 were not reproduced. We believe the code to be correct, however, and so the results can be reproduced by changing the `TOTAL_PERMUTATIONS` variable in `commands.sh`.
