#########################################################################
#                                                                       #
#    P r i o r s                                                        #
#    Priors are set here.                                               #
#                                                                       #
#########################################################################

# External modules.
import math
import scipy
from scipy import special

# Superpy modules.
import numpy as NP
import Cube


def myprior(cube, ndim, nparams):
    """ MultiNest callback function. Sets the model parameters
    from the unit hypercube.
    Arguments:
    cube -- Unit hypercube from which model parameters are mapped.
    ndim -- Number of model paramers.
    nparams -- Total number of parameters in the cube .

    Returns:

    """
    # Set-up the paramters - CMSSM model.
    Model = CMSSMModelTracker()

    # Set the model parameters from the cube.
    Model.SetParams(cube)

    # Copy them to the cube so that they are printed, in alphabetical order.
    for name in sorted(Model.param.keys(), key=str.lower):
        Cube.AddCube(cube, Model.param[name].value, Cube.label, name)

    # Print-out cube model parameters for debugging.
    print '*****************************************************'
    print 'Parameters:'
    for i in range(Cube.AddCube.count() + 1):  # Count begins at 0 - need +1.
        print Cube.label[i], cube[i]


#########################################################################

# Model class - this is where everything is stored.
# This class contains all the parameters for the model,
# and the functions to evalulate them.

class CMSSMModelTracker:

    """ Contains CMSSM model parameters, priors and functions to
    evaluate the parameters. """

    def __init__(self):
        """ Sets-up the model.
        Specify the priors for the model and nuisance parameters.
        """

        self.param = {}  # Holds model's parameters.

        # Set the priors for each parameter. NB Gaussian priors are in a way
        # infinite, so you don't need to specify upper and lower bounds.
        # NB If you use a log parameter, be sure its always strictly > 0.
        # The parameters are in GeV.

        # CMSSM parameters.

        # Universal scalar soft-breaking mass at GUT scale.
        self.param['m0'] = LogParameter(100, 4000)

        # Universal gaugino soft-breaking mass at GUT scale.
        self.param['m12'] = LogParameter(100, 2000)

        # Universal trilinear soft-breaking at GUT scale.
        self.param['a0'] = LinearParameter(-4000, 4000)

        # Ratio of Higgs VEVs, defined at EW scale.
        self.param['tanbeta'] = LinearParameter(3, 62)

        # Sign of Higgs parameter.
        self.param['signmu'] = FixedParameter(1)

        # Nusisance parameters.

        # Bottom MS-bar mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/DataBlock.action?node=Q005M
        self.param['mb'] = GaussParameter(4.18, 0.03)

        # Top pole mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/DataBlock.action?node=Q007T7
        self.param['mt'] = GaussParameter(
            173.21,
            (0.51 ** 2 + 0.71 ** 2) ** 0.5)

        # Strong coupling.
        # PDG
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-qcd.pdf
        self.param['alphas'] = GaussParameter(0.1185, 0.0006)

        # Reciprocal of EM coupling at MZ.
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-standard-model.pd
        self.param['invalpha'] = GaussParameter(127.940, 0.014)

    def SetParams(self, x):
        """ Sets the parameters from the unit hypercube.
        Arguments:
        x -- The whole unit hypercube vector.

        Returns:

        """
        for i, name in enumerate(self.param.keys()):
            self.param[name].SetValue(x[i])

#########################################################################

# These classes are all the types of prior one might use, they all have a
# SetValue function, which sets the value of parameter from MultiNest's unit
# hypercube. i.e. all priors are mapped from a Uniform(0,1) distribution.


@Cube.memoize
class LogParameter:

    """ A parameter with a log prior. """

    def __init__(self, min, max):
        """ Initializes a parameter with a log prior.
        Arguments:
        min -- Minimum of prior.
        max -- Maximum of prior.

        Returns:

        """
        self.minimum = min
        self.maximum = max
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Calculates value of parameter with log prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Log the ranges, map from 0-1 to min-max in log space,
        # then exponentiate.
        self.value = math.pow(10, math.log10(self.minimum) + x *
                              (math.log10(self.maximum) -
                               math.log10(self.minimum)))


@Cube.memoize
class LinearParameter:

    """ A parameter with a linear prior. """

    def __init__(self, min, max):
        """ Initializes a parameter with a linear prior .
        Arguments:
        min -- Minimum of prior.
        max -- Maximum of prior.

        Returns:

        """
        self.minimum = min
        self.maximum = max
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Calculates value of parameter with a linear prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Map from 0-1 to min-max.
        self.value = self.minimum + x * (self.maximum - self.minimum)


@Cube.memoize
class FixedParameter:

    """ A parameter with a fixed value. """

    def __init__(self, fixed):
        """ Initializes a parameter with a fixed value.
        You might, for example, want to fix sgn mu in the CMSSM.
        Arguments:
        fixed -- Fixed value of parameters.

        Returns:

        """
        self.fix = fixed
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter with a fixed prior, ignore unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube, ignored.

        Returns:

        """
        # NB that although x isn't used, it should stil be here
        # so that this class can be used like the others.
        self.value = self.fix
        self.arg = self.args


@Cube.memoize
class SignParameter:

    """ A parameter that is either +/- 1. """

    def __init__(self):
        """ Initializes a parameter that is either +/-1.
        For example, sign mu in the CMSSM.
        Arguments:

        Returns:

        """
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter to +/- 1, with equal probability.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # You could change this to favour a particular sign.
        if x > 0.5:
            self.value = +1
        else:
            self.value = -1


@Cube.memoize
class GaussParameter:

    """ A parameter with a Gaussian prior """

    def __init__(self, mu, sigma):
        """ Initializes a parameter with a Gaussian prior.
        Arguments:
        mu -- Mean of Gaussian.
        sigma -- Variance of Gaussian.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter with a Gaussian prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Map from 0-1 to a Gaussian! erfcinv is a the inverse of the
        # complementary Gaussian error function.
        self.value = self.mu + self.sigma * math.sqrt(2) *\
            scipy.special.erfcinv(2 * (1 - x))

@Cube.memoize
class Update:

    """ Use one-dimensional  PDF from previous scan as prior for
    a new scan! """

    def __init__(self, filename, column, weight=2, nbins=50):
        """ Load a *.txt file, import the data, and draw samples.

        Arguments:
        filename -- Name of PDF file.
        weight -- Column of PDF weight.
        column -- Column to generate point for.
        nbins -- Number of bins for histogram.

        """
        self.nbins = nbins
        import sys
        sys.path.append("../SuperPlot")
        import OneDim
        import PlotMod
        from scipy.interpolate import interp1d

        # Make one-dimensional histogram. First open chain.
        data = PlotMod.OpenChain(filename)
        # Make one-dimensional histogram of relevant data, i.e.
        # one-dimensional PDF.
        self.hist = OneDim.PosteriorPDF(
            data[column],
            data[weight],
            nbins=self.nbins)
        # Normalize PDF to that sum is one.
        self.hist.pdf = self.hist.pdf / sum(self.hist.pdf)

        # Find CDF of PDF.
        self.cdf = {}
        # Set zeroth bin.
        self.cdf[0] = self.hist.pdf[0]
        # Ignore zeroth bin. Recursively find CDF.
        for bin in range(1,self.nbins):
            self.cdf[bin] = self.cdf[bin - 1] + self.hist.pdf[bin]

        # Find inverse CDF.
        # Make grid of inverse CDF data.
        probability = NP.linspace(0, 1, 100)
        inverse_CDF = NP.zeros(len(probability))
        for i, p in enumerate(probability):
            inverse_CDF[i] = self.inverse_CDF(p)

        # Make interpolation function for inverse CDF.
        bin_width = self.hist.bins[1] - self.hist.bins[0]
        self.inverse_CDF_interp = interp1d(probability, inverse_CDF, kind='cubic')

    def inverse_CDF(self, x):
        """ Find inverse CDF.

        Arguments:
        x -- Cumulative probablity value, from 0 to 1.

        Returns:
        bin -- Bin center.

        """
        # Once specified probability exceeded, stop and
        # return center of bin.
        for bin in range(self.nbins):
            if self.cdf[bin] >= x:
                # Return bin center.
                return self.hist.bins[bin]

        # Return final bin, if no other option.
        return self.hist.bins[self.nbins - 1]

    def SetValue(self, x):
        """ Sample point from one-dimensional PDF.

        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Interpolate bin number from CDF.
        self.value = self.inverse_CDF_interp(x)