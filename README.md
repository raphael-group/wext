# Weighted Exclusivity Test (WExT) #

The Weighted Exclusivity Test (WExT) was developed by the [Raphael research group](http://compbio.cs.brown.edu/) at Brown University.

### Requirements ###

Latest tested version in parentheses.

1. Python (2.7.9)

    a. NumPy (1.11.0)

    b. SciPy (0.17.0)

2. gcc (4.9.2)

We recommend using [`virtualenv`](https://virtualenv.pypa.io/en/latest/) to install the Python requirements. After installing `virtualenv`, you can install the Python requirements for the weighted exclusivity test as follows:

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt

See the wiki for additional instructions on [Setup and installation](https://github.com/raphael-group/weighted-exclusivity-test/wiki/Setup-and-installation).

### Setup ###

The C and Fortran extensions must be compiled before running the weighted exclusivity test:

    cd wext
    python setup.py build
    f2py -c src/fortran/bipartite_edge_swap_module.f95 -m bipartite_edge_swap_module

### Usage ###

#### Data preprocessing ####
Before computing the weighted test, you need to process the input mutation data and generated permuted matrices.

1. Process mutation data in MAF format with `process_mutations.py`. See [Process mutations](https://github.com/raphael-group/wext/wiki/Process-mutations) on the wiki for details on usage and input.
2. Generate permuted versions of the mutation data -- fixing the number of mutations per gene and per patient/sample -- and compute mutation probabilities with `compute_mutation_probabilities.py`. See [Compute mutation probabilities](https://github.com/raphael-group/wext/wiki/Compute-mutation-probabilities) on the wiki for details on usage and input.

#### Searching for exclusive sets ####

Given the mutation data, we compute the exclusivity of sets _M_ of genes with `compute_exclusivity.py`. Users can choose which test (unweighted, weighted, or permutational) and, for the unweighted and weighted tests, which method (exact or saddlepoint) is used to compute the _p_-values. See [Compute exclusivity](https://github.com/raphael-group/wext/wiki/Compute-exclusivity) on the wiki for details on usage and input.

### Visualization ###

We provide scripts to run an interactive web application to view the output of `compute_exclusivity.py`, including both the mutations and mutation probabilities for each set. See [`viz/README.md`](https://github.com/raphael-group/wext/blob/master/viz/README.md) and [the wiki](https://github.com/raphael-group/wext/wiki/Visualization) for additional instructions and details.

### Testing ###
[Add testing instructions here.]

### Support ###

Please visit the [Dendrix Google Group](https://groups.google.com/forum/#!forum/dendrix) to post questions and view discussions from other users about our methods for identifying mutually exclusive mutations, or contact us through our research group's website.

### Reference ###

Mark D.M. Leiserson, Matthew A. Reyna, and Benjamin J. Raphael. (2016) A Weighted Exact Test for Mutually Exclusive Mutations in Cancer. _ECCB/Bioinformatics_ (to appear). [arXiv preprint](http://arxiv.org/abs/1607.02447)
