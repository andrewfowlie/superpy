#########################################################################
#                                                                       #
#    S t a t s                                                          #
#                                                                       #
#########################################################################

import numpy as NP

# Useful statistical functions.


def PosteriorMean(posterior, param):
    """ Calculate the posterior mean.

    Arguments:
    posterior -- Data column of posterior weight.
    param -- Data column of parameter.

    Returns:
    postmean - The posterior mean.

    """
    # Calculate posterior mean - dot product weights with parameter
    # values and normalize.
    postmean = NP.dot(posterior, param) / \
        sum(posterior)
    return postmean


def BestFit(chisq, param):
    """ Calculate the best-fit.

    Arguments:
    chisq -- Data column of chi-squared.
    param -- Data column of parameter.

    Returns:
    bestfit -- The best-fit point.

    """
    # Calculate the best-fit - find the point that corresponds
    # to the smallest chi-squared.
    bestfit = param[chisq.argmin()]
    return bestfit


def PValue(chisq, dof):
    """ Calculate the pvalue.

    Arguments:
    chisq -- Data column of chi-squared.
    dof -- Number of degrees of freedom.

    Returns:
    pvalue -- A p-value for the given chisq, dof.

    """
    from scipy import stats
    # Find the associated p-value. The survival function, sf,
    # is 1 - the CDF.
    pvalue = stats.chi2.sf(chisq.min(), dof)
    return pvalue
