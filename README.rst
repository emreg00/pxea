
# PxEA: ProXimal pathway Enrichment Analysis
--------------------------------------------

Code for the ProXimal pathway enrichment analysis introduced in the "Targeting comorbid diseases via network endopharmacology" manuscript.


## Requirements

- Python 2 or 3
- numpy
- scipy
- networkx


## Installing & running tests

Download (i.e. clone) the files to your computer, you can use the package as a bare package (without setup.py install) or install it
using the following command:

>>> python setup.py install

Several test cases for the methods are provided in `test/test_pxea.py`. 
To run these, on the parent directory (where this README file resides) type

>>> python -m unittest test.test_pxea
or
>>> python setup.py test

It should give an output similar to below
..
----------------------------------------------------------------------
Ran 2 tests in 1.220s

OK


## Usage

### PXEA

>>> from pxea.utilities.set_enrichment import get_enrichment_score_and_pval

get_enrichment_score_and_pval(ranked_list, candidates, n_random=n_random, alternative="greater", seed=51234)
    """
    KS based score (~max difference between cumulative distributions of the sample and expected random walk)
    """

Input parameters:
    ranked_list: a list with the ranking of the elements (e.g., pathways proximal to a drug)
    candidates: set of elements (e.g., pathways common to two diseases)
    N: number of pathways in the candidates set (if None, len(candidates) will be used
    n_random: number of shufflings to the ranked list for permutation test based P-value calculation
    (if none, no pvalue is calculated and None is returned instead)
    alternative: greater | less | two-sided
    seed: number to be used initializing random generator (for reproducibility)

Output:
    Returns enrichment score and pvalue


### Proximity

>>> from pxea.utilities.network import calculate_proximity

calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=None, seed=51234, lengths=None)
    """
    Calculate proximity (average distance to the closest node from the 
    first to second)from nodes_from to nodes_to (if degree binning or 
    random nodes are not given, they are generated)
    """

Input parameters:
    network: networkx Graph object
    nodes_from: set of nodes from which proximity is calculated
    nodes_to: set of nodes proximity to which is calculated
    nodes_from_random: random from nodes to check proximity
    nodes_to_random: random to nodes to check proximity
    bins: degree equivalence bins
    n_random: number of randomizations for background closest distance calculation
    min_bin_size: minimum size of the bins for degree binning if None, len(network) // 100 is used
    seed: integer for initializing the state of the random generator
    lengths: precalculated shortest path length dictionary

Output:
    Returns proximity z-score and average path length to nearest nodes in nodes_to


## Data sets

The data sets used in the analysis for the autoimmune diseases are 
under data/ folder at [pxea](https://github.com/emreg00/pxea).


## See also

See [toolbox](https://github.com/emreg00/toolbox) package for various related code.

