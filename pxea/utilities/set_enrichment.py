#######################################################################
# Functions for calculating running sum-like set enrichment score
# e.g. 04/2018
#######################################################################
from __future__ import division, print_function

import numpy
import random 
from . import stat

def get_enrichment_score_and_pval(ranked_list, candidates, N=None, n_random=None, alternative="greater", seed=51234):
    """
    KS based score (~max difference between cumulative distributions of the sample and expected random walk)
    ranked_list: a list with the ranking of the elements (e.g., pathways proximal to a drug)
    candidates: set of elements (e.g., pathways common to two diseases)
    N: number of pathways in the candidates set (if None, len(candidates) will be used
    n_random: number of shufflings to the ranked list for permutation test based P-value calculation
    (if none, no pvalue is calculated and None is returned instead)
    alternative: greater | less | two-sided
    seed: number to be used initializing random generator (for reproducibility)
    RETURNS enrichment score and pvalue
    """
    ks, pval = None, None
    if len(candidates) == 0:
        return ks, pval 
    if seed is not None:
        random.seed(seed)
    ks = stat.ks_score(ranked_list, candidates, N)
    if n_random is not None:
        values = []
        for i in range(n_random):
            random.shuffle(ranked_list)
            ks_random = stat.ks_score(ranked_list, candidates, N)
            values.append(ks_random)
    if n_random is not None:
        values = numpy.array(values)
        #z = (ks - numpy.mean(values)) / numpy.std(values, ddof=1)
        #print(z, numpy.mean(values), numpy.median(values)) #, values
        if alternative == "greater":
            pval = sum(values >= ks) / len(values)
        elif alternative == "less":
            pval = sum(values <= ks) / len(values)
        elif alternative == "two-sided":
            pval = sum(abs(values) >= abs(ks)) / len(values)
        else:
            raise ValueError("Unknown alternative: %s" % alternative)
    return ks, pval


