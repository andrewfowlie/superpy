#########################################################################
#                                                                       #
#    2 D    S t a t s                                                   #
#                                                                       #
#########################################################################

# This class contains all the functions for analysing a chain (*.txt file)
# and calculating the 2D stats for a particular pair of variables.

# Internal modules.
import Stats

# External modules.
import numpy as NP
from pylab import *
from scipy import stats


class PosteriorPDF:

    """ Class for calculating two-dimensional posterior PDF. """

    def __init__(self, paramx, paramy, posterior, nbins=50, bin_limits=None):
        """ Histograms the chosen parameters to obtain two-dimensional PDF.

        Arguments:
        paramx -- Data column of parameter x.
        paramy -- Data column of parameter y.
        posterior -- Data column of posterior weight.
        nbins -- Number of bins for histogram.
        bin_limits -- Bin limits for histogram.

        """
        # Outliers sometimes mess up bins. So you might want to
        # specify the bin ranges.
        # 2D Histogram the data - pdf is a matrix.
        pdf, bin_edgesx, bin_edgesy = NP.histogram2d(
            paramx, paramy, nbins,
            range=bin_limits, weights=posterior)
        # Normalize the pdf, so that its maximum value is one.
        # NB will later also normalize so that area is one.
        self.pdf = pdf / pdf.max()
        # Find centres of bins.
        self.centerx = (bin_edgesx[:-1] + bin_edgesx[1:]) * 0.5
        self.centery = (bin_edgesy[:-1] + bin_edgesy[1:]) * 0.5


class ProfileLike:

    """ Class for calculating two-dimensional profile likelihood. """

    def __init__(self, paramx, paramy, chisq, nbins, bin_limits=None):
        """ Maximizes the chisquared to obtain two-dimensional profile likelihood.

        Arguments:
        paramx -- Data column of parameter x.
        paramy -- Data column of parameter y.
        chisq -- Data column of chisq.
        nbins -- Number of bins for histogram.
        bin_limits -- Bin limits for histogram.

        """
        # Bin the data - digitialize will return a column vector
        # containing the bin number for each point in the chain.
        # NB, we discard the pdf.
        pdf, bin_edgesx, bin_edgesy = NP.histogram2d(
            paramx, paramy, nbins,
            range=bin_limits, weights=None)
        bin_numbersx = NP.digitize(paramx, bin_edgesx)
        bin_numbersy = NP.digitize(paramy, bin_edgesy)

        # Subract one from the bin numbers, so that bin numbers,
        # initially (0, nbins+1) because numpy uses extra bins for outliers,
        # match array indices (0, nbins-1).

        # First deal with outliers.
        for i in range(bin_numbersx.size):
            if bin_numbersx[i] == 0:
                bin_numbersx[i] = 1
            if bin_numbersx[i] == nbins + 1:
                bin_numbersx[i] = nbins
            if bin_numbersy[i] == 0:
                bin_numbersy[i] = 1
            if bin_numbersy[i] == nbins + 1:
                bin_numbersy[i] = nbins

        # Now subtract one.
        bin_numbersx = bin_numbersx - 1
        bin_numbersy = bin_numbersy - 1

        # Initialize the profiled chi-squared to something massive.
        profchisq = NP.zeros((nbins, nbins)) + 1e90

        # Min the chi-squared in each bin.
        for i in range(chisq.size):
            # If the chi-squared is less than the current chi-squared.
            if chisq[i] < profchisq[bin_numbersx[i], bin_numbersy[i]]:
                profchisq[bin_numbersx[i], bin_numbersy[i]] = chisq[i]

        # Now exponentiate to obtain likelihood, and normalize.
        self.profchisq = profchisq - profchisq.min()
        self.proflike = NP.exp(- self.profchisq * 0.5)

        # Find centres of bins.
        self.centerx = (bin_edgesx[:-1] + bin_edgesx[1:]) * 0.5
        self.centery = (bin_edgesy[:-1] + bin_edgesy[1:]) * 0.5


class CredibleRegions:

    """ Class for calculating two-dimensional credible regions. """

    def __init__(self, pdf, epsilon=NP.array([0.05, 0.32])):
        """ Calculate the credible regions from marginalised pdf.

        Arguments:
        pdf -- Marginalised two-dimensional posterior PDF.
        epsilon -- Probability levels.

        """
        # Smallest possible regions that contain given amound of PDF.
        # Sum pdf over pdf > crediblelevel, increasing crediblelevel
        # from 0 until sum pdf < 1 - epsilon.
        # We only need the crediblelevel to plot a contour.

        # Normalize pdf so that area is one.
        pdf = pdf / pdf.sum()

        # Set the crediblelevels to zero.
        self.crediblelevel = NP.zeros(epsilon.size)

        # Loop over contours required.
        for i in range(epsilon.size):
            # While the sum of the PDF over the regions with density
            # > credible level is greater than 1 - epsilon,
            # increase the credible level slowly.
            # Once this exits, we have found credible level
            while ma.masked_where(pdf < self.crediblelevel[i],
                                  pdf).sum() > 1 - epsilon[i]:
                self.crediblelevel[i] += 0.0001
        # NB that you could increase accuracy by decreasing the increment.
        # If contours don't appear on a plot, try altering that.


class ConfidenceIntervals:

    """ Class for calculating two-dimensional confidence limits. """

    def __init__(self, epsilon=NP.array([0.05, 0.32])):
        """ Calculate the confidence intervals.

        Arguments:
        chisq -- Data column of profiled chisq.
        epsilon -- Confidence levels.

        """
        # First invert the epsilons to delta chi-squareds with inverse
        # cumalative chi-squared distribution with 2 dof.
        deltachisq = stats.chi2.ppf(1 - epsilon, 2)
        # Convert these into PL values.
        self.deltaPL = NP.exp(- deltachisq / 2)

        # That's all we need! - we will simply plot contours of
        # deltaPL.
