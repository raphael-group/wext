# Weighted Exclusivity Test #

The weighted exclusivity test was developed by the [Raphael research group](http://compbio.cs.brown.edu/) at Brown University.

### Requirements ###

Latest tested version in parentheses.

1. Python (2.7.9)

    a. NumPy (1.11.0)

    b. SciPy (0.17.0)

    c. Matplotlib (1.5.1)

    d. Seaborn (0.7.0)

    e. Pandas (0.18.0)

    f. NetworkX (1.11) -- optional as it is required for CoMEt, see below

2. R (3.1.1)

    a. BiReWire (1.8.0)
3. gcc (4.9.2)

We recommend using [`virtualenv`](https://virtualenv.pypa.io/en/latest/) to install the Python requirements. After installing `virtualenv`, you can install the Python requirements for the weighted exclusivity test as follows:

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt

You can install the BiRewire R module following the instructions on [their website](https://www.bioconductor.org/packages/release/bioc/html/BiRewire.html).

[Point to Wiki for additional information]

#### Optional ####

The [CoMEt Python module](https://github.com/raphael-group/comet) is required to compute the exact version of the unweighted test. By default, the weighted exclusivity test expects CoMEt to be install in `third-party/comet`. You can install and setup the latest version of CoMEt using the following commands:

    mkdir third-party
    cd third-party
    git clone https://github.com/raphael-group/comet.git
    cd comet
    python setup.py build

### Setup ###

The C extensions must be compiled before running the weighted exclusivity test:

    cd weighted_exclusivity_test
    python setup.py build

### Usage ###
[Add usage here. Point to Wiki.]

### Testing ###
[Add testing instructions here.]

### Support ###

Please visit the [Dendrix Google Group](https://groups.google.com/forum/#!forum/dendrix) to post questions and view discussions from other users about our methods for identifying mutually exclusive mutations, or contact us through our research group's website.

### Reference ###

Mark D.M. Leiserson, Matthew A. Reyna, and Benjamin J. Raphael. (2016) A Weighted Exact Test for Mutually Exclusive Mutations in Cancer. *In submission*.

