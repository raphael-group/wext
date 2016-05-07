# Weighted Exclusivity Test #

The weighted exclusivity test was developed by the [Raphael research group](http://compbio.cs.brown.edu/) at Brown University.

### Requirements ###

Latest tested version in parentheses.

1. Python (2.7.9)

    a. NumPy (1.11.0)

    b. SciPy (0.17.0)

2. R (3.1.1)

    a. BiReWire (1.8.0)
3. gcc (4.9.2)

We recommend using [`virtualenv`](https://virtualenv.pypa.io/en/latest/) to install the Python requirements. After installing `virtualenv`, you can install the Python requirements for the weighted exclusivity test as follows:

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt

You can install the BiRewire R module following the instructions on [their website](https://www.bioconductor.org/packages/release/bioc/html/BiRewire.html).

See the wiki for additional instructions on [Setup and installation](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Setup-and-installation).

### Setup ###

The C extensions must be compiled before running the weighted exclusivity test:

    cd weighted_exclusivity_test
    python setup.py build

### Usage ###

#### Data preprocessing ####
Before computing the weighted test, you need to process the input mutation data and generated permuted matrices.

1. Process mutation data in MAF format with `process_mutations.py`. See [Process mutations](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Process-mutations) on the wiki for details on usage and input.
2. Generate permuted versions of the mutation data with `permute_matrix.py`. See [Permute mutation data](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Permute-mutation-data) on the wiki for details on usage and input.
3. Compute mutation probabilities using the permuted mutation data with `compute_mutation_probabilities.py`. See [Compute mutation probabilities](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Compute-mutation-probabilities) on the wiki for details on usage and input.

#### Computing the weighted test ####

Given the mutation data, we compute the exclusivity of sets _M_ of genes with `compute_exclusivity.py`. Users can choose which test (unweighted, weighted, or permutational) and, for the unweighted and weighted tests, which method (exact or saddlepoint) is used to compute the _p_-values. See [Compute exclusivity](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Compute-exclusivity) on the wiki for details on usage and input.

### Visualization ###

We provide scripts to run an interactive web application to view the output of `compute_exclusivity.py`, including both the mutations and mutation probabilities for each set. See [`viz/README.md`](https://github.com/raphael-group/weighted-exclusivity-test/blob/master/viz/README.md) and [the wiki](https://github.com/raphael-group/weighted-exclusivity-test/wiki/VIsualization) for additional instructions and details.

### Testing ###
[Add testing instructions here.]

### Support ###

Please visit the [Dendrix Google Group](https://groups.google.com/forum/#!forum/dendrix) to post questions and view discussions from other users about our methods for identifying mutually exclusive mutations, or contact us through our research group's website.

### Reference ###

Mark D.M. Leiserson, Matthew A. Reyna, and Benjamin J. Raphael. (2016) A Weighted Exact Test for Mutually Exclusive Mutations in Cancer. *In submission*.
