# Running WExT WRE on Pan-Cancer data (`pancan`) #

### Description ###
This simple example generates two datasets -- which can be interpreted as data from two cancer types -- and then runs the WExT weighted row exclusivity (WRE) test on the combined, pan-cancer dataset in two modes. These two modes were inspired by [Y.A. Kim, et al. (_ISMB/Bioinformatics_, 2015)](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv247).

In the first mode, we are looking for pan-cancer exclusivity where we are controlling for cancer type. In other words, we _don't_ want to find exclusivity that is driven by cancer type, i.e. where one gene is only mutated in cancer type A and a second gene is mutated only in cancer type B. Kim, et al. call this _across_ mutual exclusivity ("ACROSS\_ME"), although their ACROSS\_ME also requires that the exclusivity isn't present in a single cancer type.

To run WExT in this mode, we process mutations and generate mutation probabilities for each dataset/cancer type separately, and then combine these when computing the WRE test across both datasets. 

In the second mode, we are looking for pan-cancer exclusivity without controlling for cancer type. This is called _between mutual exclusivity_ ("BETWEEN\_ME") by Kim, et al. To perform this type of run in WExT, the datasets need to be combined 

### Data ###
To demonstrate these two modes, we generate simulated data. The data includes two implanted gene sets mutated in approximately 25% of samples.
* `(G1, G2)`: this set has "across" mutual exclusivity, i.e. exclusivity where both genes are mutated in both datasets/cancer types.
* `(G3, G4)`: this set has "between" mutual exclusivity, i.e. each gene is mutated in only one dataset/cancer type.

The datasets also include other genes that are mutated in each sample at a fixed background mutation rate. We ensure that every sample has at least one mutation.

### Usage ###

Run `make` to process the mutations, compute mutation probabilities, and run the WRE test in both modes described above. You can find the two output TSV files in `output/`. 

You can pass in a few optional parameters to `make`, including:
* `K`: gene set size (default: 2)
* `NP`: number of permutations (default: 100)
* `NUM_GENES`: number of genes in simulated dataset (default: 50)
* `NUM_SAMPLES`: number of samples in simulated dataset (default: 100)
* `RANDOM_SEED`: random seed for PRNG (default: 2)
* `BMR`: background mutation rate (default: 0.01)

For example, to run the example searching for sets of size three with 100 permutations:

    make K=3 NP=100

Both these examples -- using `K=2` or `K=3` -- take less than 5 minutes when running on a single core on a modern machine.

The easiest way to interpret the results is to compare the set P-value `(G3, G4)` with between exclusivity in the two output files in `output/`. The P-value for `(G3, G4)` should be smaller (more significant) in the "between" output file.