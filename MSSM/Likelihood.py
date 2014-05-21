#########################################################################
#                                                                       #
#    L i k e l i h o o d                                                #
#                                                                       #
#########################################################################

# External modules.
import os
import math
import scipy
import subprocess
import pyslha
import tempfile
import numpy as NP
import sys
import StringIO

# SuperPy modules.
import Priors
import Cube
import fastlim_superpy


def myloglike(cube, ndim, nparams):
    """ MultiNest callback function.
    Calculate the model's predictions (by calling external programs),
    saves the extra parameters in the cube, and returns log likelihood to Multinest.

    Arguments:
    cube -- Unit hypercube from which model parameters are mapped.
    ndim -- Number of model paramers.
    nparams -- Total number of parameters in the cube.

    Returns: The total log likelihood for the model parameter
    point under consideration.

    """
    # Set up constraints class.
    Constraints = CMSSMConstraintTracker()

    # Copy cube to constraints, so it can work out predictions etc.
    for i, name in enumerate(sorted(Priors.CMSSMModelTracker().param.keys(), key=str.lower)):
        Constraints.param[name] = cube[i]

    # Set predictions and loglikes.
    Constraints.SetPredictions()
    Constraints.SetLogLike()

    # Copy constraints to cube.
    for name in sorted(Constraints.constraint.keys(), key=str.lower):
        Cube.AddCube(
            cube,
            Constraints.constraint[name].theory,
            Cube.label,
            name)

    # Copy associated chi2s to cube. Better to print chi2 than loglike,
    # beceause MultiNest prints chi2 and they can be treated in the
    # same way when plotting.
    for name in sorted(Constraints.constraint.keys(), key=str.lower):
        Cube.AddCube(cube, -
                     2 *
                     Constraints.constraint[name].loglike, Cube.label, 'chi2:' +
                     name)

    # Copy SLHA masses to the cube.
    for key in Constraints.masses:
        Cube.AddCube(
            cube,
            Constraints.masses[key],
            Cube.label,
            'Mass:' +
            str(key))

    # Copy mu-parameter to the cube.
    Cube.AddCube(cube, Constraints.mu, Cube.label, 'mu')

    # Copy neutralino mixing to the cube.
    for key in Constraints.neutralino:
        Cube.AddCube(
            cube,
            Constraints.neutralino[key],
            Cube.label,
            'NMIX:' +
            str(key))

    # Print-out cube for debugging.
    print 'Predictions:'
    for label, param in zip(Cube.label.itervalues(), cube):
        print label, param
    print 'Total loglike', Constraints.loglike

    # Start cube count from 0 again.
    Cube.AddCube.reset()

    # Print an info file for the cube.
    # Decorator insures this is printed only once.
    # Point must be physical, else the info will be incomplete.
    if Constraints.physical:
        Cube.PrintInfo()

    # Return the log likelihood to MultiNest.
    return Constraints.loglike


#########################################################################

# This class evalulates likelihoods and the model's predictions etc.

class CMSSMConstraintTracker:

    """ Contains CMSSM constraints and functions to
    evaluate the model's predictions. """

    def __init__(self):
        """ Sets-up the model. """
        # Dictionary of model's parameters.
        self.param = {}
        # Dictionary of constraint information.
        self.constraint = {}
        # Whether the model has an acceptable mass spectrum, EWSB etc.
        self.physical = True
        # The total loglike associated with model.
        self.loglike = -1e101
        # List SLHA masses.
        self.masses = []

        # You set the constraint values here. You must make sure that
        # you have computed all the relevant theoretical values etc.

        # Naturalness priors.
        # Implemented in modified SOFTSUSY.
        # Although it is a prior, easier to put here in the likelihood.
        self.constraint['Natural'] = ExternalConstraint()

        # LEP, Higgs and Tevatron data.
        # Implemented in HiggsSignals.
        # http://higgsbounds.hepforge.org/
        self.constraint['Higgs'] = ExternalConstraint()

        # LHC direct searches.
        # Implemented in Fast-Lim.
        # http://fastlim.web.cern.ch/fastlim/
        self.constraint['LHC'] = ExternalConstraint()

        # Interpolate lower bound on (m0,m12) plane.
        # https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CombinedSummaryPlots/SUSY/ATLAS_SUSY_MSUGRA/ATLAS_SUSY_MSUGRA.png
        # ATLAS-CONF-2013-047
        self.constraint['LHC_interp'] = InterpolateLowerConstraint(
            'atlas_m0m12.dat',
            0.01)

        # Relic density of neutralinos.
        # Planck.
        # http://arxiv.org/pdf/1303.5076v2.pdf
        self.constraint['oh2'] = GaussConstraintFractionalTau(
            0.1199,
            0.0027,
            0.1)

        # Anomalous magnetic moment of muon.
        # PDG.
        # http://pdg.lbl.gov/2013/reviews/rpp2013-rev-g-2-muon-anom-mag-moment.pdf
        self.constraint['gm2'] = GaussConstraint(28.8e-10, 8.0e-10, 1e-10)

        # BR(b -> s gamma).
        # HFAG.
        # http://www.slac.stanford.edu/xorg/hfag2/rare/2013/radll/OUTPUT/HTML/radll_table3.html
        self.constraint['bsg'] = GaussConstraint(3.43e-4, 0.22e-4, 0.21e-4)

        # BR(Bs -> mu mu).
        # PDG.
        # http://pdg8.lbl.gov/rpp2013v2/pdgLive/Particle.action?node=S086#decays
        self.constraint['bsmumu'] = GaussConstraintFractionalTau(
            3.2e-9,
            1.5e-9,
            0.14)

        # BR(b -> tau nu).
        # HFAG.
        # http://www.slac.stanford.edu/xorg/hfag2/rare/2013/radll/OUTPUT/HTML/radll_table7.html
        self.constraint['btaunu'] = GaussConstraint(1.14e-4, 0.22e-4, 0.38e-4)

        # Spin-indenpendent WIMP-nucleon scattering cross section.
        # LUX.
        # http://arxiv.org/abs/1310.8214
        self.constraint['sigsip'] = InterpolateUpperConstraint('lux.dat', 10)

        # W-boson mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2013v2/pdgLive/Particle.action?node=S043
        self.constraint['mw'] = GaussConstraint(80.385, 0.015, 0.015)

        # Leptonic sin eff theta, effective weak-mixing angle.
        # PDG.
        # http://pdg.lbl.gov/2013/reviews/rpp2013-rev-standard-model.pdf
        self.constraint['sineff'] = GaussConstraint(0.23146, 0.00012, 15e-5)

        # delta MBs.
        # HFAG.
        # http://www.slac.stanford.edu/xorg/hfag/osc/PDG_2014/#DMS
        self.constraint['deltaMb'] = GaussConstraint(17.761, 0.022, 2.4)

    def SetPredictions(self):
        """ Run the auxilliary programs for a particular model
        point to find the model's preddictions.
        Arguments:

        Returns:

        """

        # Call SOFTSUSY and read the mass spectrum.
        print "Calling SOFTSUSY..."
        self.softsusy()
        self.readslha()

        # Set LHC interpolation parameters.
        self.constraint['LHC_interp'].theory = self.param['m0']
        self.constraint['LHC_interp'].theory = self.param['m12']

        # Call auxillary programs if physical.
        if self.physical:
            print "Calling SUSY-HIT..."
            self.susyhit()
        if self.physical:
            print "Calling Fast-Lim..."
            self.fastlim()
        if self.physical:
            print "Calling FeynHiggs..."
            self.feynhiggs()
        if self.physical:
            print "Calling micrOMEGAs..."
            self.micromegas()
        if self.physical:
            print "Calling SuperISO..."
            self.superiso()
        if self.physical:
            print "Calling HiggsSignals..."
            self.higgssignals()
        if self.physical:
            print "Finding naturalness priors..."
            self.naturalness()

    def softsusy(self):
        """Call SoftSUSY to obtain predictions for model.
        Sets SLHA - the name of an SLHA mass spectrum for the
        model point.
        Arguments:

        Returns:

        """
        # Create SLHA input file.
        self.SLHAIN = self.writeslha(self.param)

        # Need to run through shell for < redirect.
        self.SLHA = RunProgram(
            '../softsusy-*/softpoint.x leshouches <' +
            self.SLHAIN,
            './',
            '',
            shell=True)
        self.physical = CheckProgram(
            self.SLHA, ["problem", "invalid", "warning", "Incorrect"])

    def writeslha(self, param, MZ=9.11876000e+01):
        """ Write an SLHA input file, SLHAIN, for a given
        parameter point.

        Arguments:
        param -- Dictionary of model parameters.
        MZ -- Value of Z-boson mass.

        Return:
        Name of SLHA input file.

        """
        SLHAIN = tempfile.NamedTemporaryFile(delete=False)
        # Pass information in SLHA format. NB that you could easily add more
        # parameters to the EXTPAR section, to relax the CMSSM.
        SLHAIN.write(
            """
Block MODSEL                    # Select model
 1   1                          # sugra
Block SMINPUTS                  # Standard Model inputs
 1   %s                         # alpha^(-1) SM MSbar(MZ)
 2   1.16637000e-05             # G_Fermi
 3   %s                         # alpha_s(MZ) SM MSbar
 4   %s                         # MZ(pole)
 5   %s                         # mb(mb) SM MSbar
 6   %s                         # mtop(pole)
 7   1.77700000e+00             # mtau(pole)
Block MINPAR    # Input parameters
 1   %s                         # m0
 2   %s                         # m12
 3   %s                         # tanb
 4   %s                         # sign(mu)
 5   %s                         # A0
Block SOFTSUSY                  # SOFTSUSY specific inputs
 1   1.000000000e-03            # tolerance
 2   2.000000000e+00            # up-quark mixing (=1) or down (=2)
 5   1.000000000E+00            # 2-loop running
 3   0.000000000E+00            # printout
 """ %
            (str(
                param['invalpha']), str(
                param['alphas']), str(
                MZ), str(
                param['mb']), str(
                param['mt']), str(
                param['m0']), str(
                param['m12']), str(
                param['tanbeta']), str(
                param['signmu']), str(
                param['a0'])))

        # Close file so that output is flushed.
        SLHAIN.close()
        print "SLHA input file:", SLHAIN.name
        return SLHAIN.name

    def readslha(self):
        """ Read an SLHA file with PySLHA.
        Populates masses, mu and neutralino mixings.
        If the point is unphysical,
        populate the mass blocks with zeros.
        Arguments:

        Returns:

        """
        blocks = {}
        try:
            # Read the blocks in the SLHA file.
            self.blocks, self.decays = pyslha.readSLHAFile(self.SLHA)
        except:
            # With expected running, shouldn't get any problems. But best
            # to be defensive. A missing mass block would cause an
            # exception, but e.g. stau LSP would not.
            self.physical = False
            print 'Caught trouble in the SLHA file.'

            # Still need to return data of correct length.
            self.masses = [0] * 33
            self.mu = 0.
            self.neutralino = [0] * 16
        else:
            # Save the mass block.
            self.masses = self.blocks['MASS'].entries
            # Pick out mu-parameter.
            self.mu = self.blocks['HMIX'].entries[1]
            # Pick out neutralino mixing.
            self.neutralino = self.blocks['NMIX'].entries

    def naturalness(self, MZ=9.11876000e+01):
        """ Calculate the naturalness priors.
        We have modifed SOFTSUSY to return the
        relevant derivatives and parameters to the SLHA file.

        We are calculating the Jacobain for (mu, b) -> (tanbeta, MZ).
        J = | d mu / d MZ | * | d b / d tanbeta |

        We also add 1/x factors for logarithmic priors for the
        mu and b parameters. We have:

        Effective prior = J * 1/mu * 1/b

        Arguments:
        MZ -- Value of the Z-boson mass.

        """

        # Read relevant derivatives and parameters.
        # Multiply by 2*MZ - we have dMZ^2, want dMZ.
        mu_MZ = 2. * MZ / ReadParameter(self.SLHA, '# dMZ^2/dmu=')
        b_TanBeta = 1. / ReadParameter(self.SLHA, '# dtanbeta/dm3sq=')
        mu = ReadParameter(self.SLHA, '# mu=')
        b = ReadParameter(self.SLHA, '# m3sq=')

        # Set the "likelihood" associated with naturalness.
        # Of course, this is NOT a true chi-squared/likelihood from an experiment.
        # But it is convenient to treat it as such.
        try:
            self.constraint['Natural'].loglike = NP.log10(
                abs(mu_MZ / mu *
                    b_TanBeta / b))
        except:
            self.constraint['Natural'].loglike = -1e101

    def fastlim(self):
        """ Call Fast-Lim and note whether point rejected.
        Fast-lim requires decay tables from SUSY HIT.
        Arguments:

        Returns:

        """

        # Wrap around call to save "print" statements in a file.
        stdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        # Find whether rejected, and if so by which analysis.
        reject, name = fastlim_superpy.excluded(self.SLHA_DECAY)
        saved = sys.stdout.getvalue()
        sys.stdout = stdout

        outputfile = tempfile.NamedTemporaryFile(delete=False)
        outputfile.write(saved)
        print "Output file name", outputfile.name

        if not reject:
            self.constraint['LHC'].loglike = 0.
            print "Point accepted by Fast-Lim."
        else:
            print "Point rejected by Fast-Lim analysis:", name
            self.constraint['LHC'].loglike = -1e101

    def micromegas(self):
        """Call micrOMEGAs to obtain DM predictions
        for model.
        Arguments:

        Returns:

        """
        filename = RunProgram(
            './main',
            '../micromegas_3.6.9.2/MSSM',
            self.SLHA)
        self.physical = CheckProgram(filename, ["error"])
        if self.physical:
            self.constraint['oh2'].theory = ReadParameter(
                filename,
                'Xf=',
                split='Omega=')
            # NB converted from pb to cm2.
            self.constraint['sigsip'].theory = ReadParameter(
                filename,
                'proton  SI',
                split='SD') * 1E-36
            self.constraint['sigsip'].theoryx = ReadParameter(
                filename,
                'Dark matter candidate',
                split='mass=')

    def superiso(self):
        """Call SuperIso to obtain B-physics and g-2 predictions for model.
        Arguments:

        Returns:

        """
        filename = RunProgram('./slha.x', '../superiso_v3.3', self.SLHA)
        self.physical = CheckProgram(filename, ["error"])
        if self.physical:
            self.constraint['bsg'].theory = ReadParameter(
                filename,
                'BR(b->s gamma)')
            self.constraint['btaunu'].theory = ReadParameter(
                filename,
                'BR(B->tau nu)')
            self.constraint['bsmumu'].theory = ReadParameter(
                filename,
                'BR(Bs->mu mu)')
            self.constraint['gm2'].theory = ReadParameter(filename, 'a_muon')

    def feynhiggs(self):
        """Call FeynHiggs to obtain EWPO predictions
        for model.
        Arguments:

        Returns:

        """
        filename = RunProgram('./SuperPyFH', '../FeynHiggs-2.10.0', self.SLHA)
        self.physical = CheckProgram(filename, ["error"])
        if self.physical:
            self.constraint['mw'].theory = ReadParameter(filename, 'MWMSSM=')
            self.constraint['sineff'].theory = ReadParameter(
                filename,
                'SW2MSSM=')
            self.constraint['deltaMb'].theory = ReadParameter(
                filename,
                'deltaMsMSSM=')

    def higgssignals(self):
        """Call HiggsSignals to find whether points is excluded by
        LEP, Tevatron, LHC Higgs searches.
        Arguments:

        Returns:

        """
        filename = RunProgram(
            './SuperPy',
            '../HiggsSignals-1.2.0/example_programs',
            self.SLHA)
        self.physical = CheckProgram(filename, ["error"])
        if self.physical:
            self.constraint['Higgs'].loglike = -0.5 * \
                ReadParameter(filename, 'chi^2 (total) =')

    def susyhit(self):
        """Call SUSY-HIT to calculate decay tables.
        Arguments:

        Returns:

        """
        filename = RunProgram('./run', '../susyhit', self.SLHA)
        self.physical = CheckProgram(filename, ["error"])
        if self.physical:
            self.SLHA_DECAY = filename

    def SetLogLike(self):
        """ Loops over the constraints, calculates each log likelihood
        and sums the log likelihoods.
        Arguments:

        Returns:

        """
        if self.physical:
            # Set the loglikelhoods from the constraints's SetLogLike
            # functions.
            self.loglike = 0.
            for name in self.constraint.keys():
                if self.constraint[name].apply:
                    self.loglike = self.loglike + \
                        self.constraint[name].SetLogLike()
        else:
            # If the point is unphysical, all loglikelihoods = biggest possible.
            # Logzero is -1e101.
            for name in self.constraint.keys():
                self.constraint[name].loglike = -1e101
            self.loglike = -1e101

#########################################################################

# These classes are various likelihoods for experimental measurements.
# I have included Gaussians, error functions etc. They all have a SetLogLike
# function, from which the loglikelhood is calculated.


@Cube.memoize
class GaussConstraint:

    """ A Gaussian likelihood constraint. """

    def __init__(self, mu, sigma, tau=0, apply=True):
        """ Initializes a Gaussian constraint.
        Arguments:
        mu -- Experimental mean of measurement.
        sigma -- Experimental variance of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Caclulate Gaussian likelihood from model's predictions and
        data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        self.loglike = - math.pow(self.mu - self.theory, 2) /\
            (2 * (math.pow(self.sigma, 2) +
                  math.pow(self.tau, 2)))
        return self.loglike


@Cube.memoize
class GaussConstraintFractionalTau:

    """ A Gaussian likelihood constraint, with a fractional
    theoretical error.
    """

    def __init__(self, mu, sigma, tau, apply=True):
        """ Initializes a Gaussian constraint, with a fractional
        theoretical error.
        Arguments:
        mu -- Experimental mean of measurement.
        sigma -- Experimental variance of measurement.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Caclulates Gaussian likelihood from model's predictions and
        data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that here theory error = tau * theoretical value.
        self.loglike = - math.pow(self.mu - self.theory, 2) /\
            (2 * (math.pow(self.sigma, 2) +
                  math.pow(self.tau * self.theory, 2)))
        return self.loglike


@Cube.memoize
class UpperConstraint:

    """ Upper bound, Gaussian error function constraint. """

    def __init__(self, limit, tau, apply=True):
        """ Initializes an upper bound constraint .
        Arguments:
        limit -- Upper bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that erfc is a complementary error function.
        try:
            self.loglike = math.log(
                0.5 * scipy.special.erfc((self.theory - self.limit) / (math.sqrt(2) * self.tau)))
        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101
        return self.loglike


@Cube.memoize
class UpperConstraintFractionalTau:

    """ Upper bound, Gaussian error function constraint, with fractional
    theory error.
    """

    def __init__(self, limit, tau, apply=True):
        """ Initializes an upper bound constraint .
        Arguments:
        limit -- Upper bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that erfc is a complementary error function.
        # NB that here theory error = tau * theoretical value.
        try:
            self.loglike = math.log(0.5 *
                                    scipy.special.erfc((self.theory -
                                                        self.limit) /
                                                       (math.sqrt(2) *
                                                        self.tau *
                                                        self.theory)))
        # Sometimes we are trying to take log zero.
        except ValueError as ZeroDivisionError:
            self.loglike = -1e101
        return self.loglike


@Cube.memoize
class LowerConstraint:

    """ Lower bound, Gaussian error function constraint. """

    def __init__(self, limit, tau, apply=True):
        """ Initializes a lower bound constraint.
        Arguments:
        limit -- Lower bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that erf is an error function.
        try:
            self.loglike = math.log(0.5 * (1 + scipy.special.erf((
                self.theory - self.limit) /
                (math.sqrt(2) * self.tau))))
        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101
        return self.loglike


@Cube.memoize
class InterpolateUpperConstraint:

    """ Interpolate am upper 2D limit from a data file. """

    def __init__(self, file, tau, apply=True):
        """ Initializes a 2D upper bound constraint.
        The constraint is on the y-co-ordinate, and is calculated as
        a function of x.

        Arguments:
        file -- Name of the data file.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.tau = tau
        self.apply = apply

        # The x and y theory values.
        self.theoryx = 0.
        self.theory = 999.

        # The derived limit on the y-parameter, which is a function of x.
        self.limit = 0.

        self.loglike = -1e101
        self.arg = self.args

        # Import the data from file. Do this here so
        # we only do it once.
        # The data should be x, y.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # Interpolate the y-value of the limit
        # corresponding to the theory x-value.
        self.limit = NP.interp(self.theoryx, self.data[0],
                               self.data[1], left=None, right=None)

        # Now calcualte likelihood with Gaussian error function - erf.
        like = 0.5 - 0.5 * \
            scipy.special.erf((self.theory - self.limit) / (2. ** 0.5 * self.tau * self.theory))
        try:
            self.loglike = math.log(like)
        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101

        return self.loglike


@Cube.memoize
class InterpolateLowerConstraint:

    """ Interpolate a lower 2D limit from a data file. """

    def __init__(self, file, tau, apply=True):
        """ Initializes a 2D lower bound constraint.
        The constraint is on the y-co-ordinate, and is calculated as
        a function of x.

        Arguments:
        file -- Name of the data file.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.tau = tau
        self.apply = apply

        # The x and y theory values.
        self.theoryx = 0.
        self.theory = 999.

        # The derived limit on the y-parameter, which is a function of x.
        self.limit = 0.

        self.loglike = -1e101
        self.arg = self.args

        # Import the data from file. Do this here so
        # we only do it once.
        # The data should be x, y.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # Interpolate the y-value of the limit
        # corresponding to the theory x-value.
        self.limit = NP.interp(self.theoryx, self.data[0],
                               self.data[1], left=None, right=None)

        # Now calcualte likelihood with Gaussian error function - erf.
        like = 0.5 + 0.5 * \
            scipy.special.erf((self.theory - self.limit) / (2. ** 0.5 * self.tau * self.theory))
        try:
            self.loglike = math.log(like)
        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101

        return self.loglike


@Cube.memoize
class LikeMapConstraint:

    """ Interpolate likelihood data file. """

    def __init__(self, file, apply=True):
        """ Interpolates a likelihood from a data file
        Arguments:
        file -- Name of the data file.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.apply = apply
        # The x parameter.
        self.theory = 0.
        # The y parameter.
        self.theoryy = 0.
        self.loglike = -1e101
        self.arg = self.args
        # Load file now so only loaded once.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        try:
            # Interpolate the likelihood from the datafile.
            self.loglike = math.log(mlab.griddata(
                self.data[0], self.data[1], self.data[2],
                NP.array([self.theory]), NP.array([self.theoryy]),
                interp='nn'))

        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101

        return self.loglike


@Cube.memoize
class ExternalConstraint:

    """ External contribution to likelihood - the likelihood is read from an external
    program.
    """

    def __init__(self, apply=True):
        """ Dummy constraint with no contribution to the likelhoood.
        Arguments:
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ External constraint - nothing to compute.
        Arguments:

        Returns: The loglike from this constraint.

        """
        return self.loglike

#########################################################################
# These are functions for manipulating programs and datafiles.


def RunProgram(executable, path, arguments, shell=False):
    """ Call a program
    and save its output to a file.

    Arguments:
    executable -- Name of executable file.
    path -- Path to executable.
    arguments -- Command-line arguments for executable.
    shell -- Whether to run through a shell, sometimes useufl.

    Returns:
    Name of file containing program output.

    """
    command = [executable, arguments]
    outputfile = tempfile.NamedTemporaryFile(delete=False)
    errorfile = tempfile.NamedTemporaryFile(delete=True)
    try:
        # Run program.
        program = subprocess.Popen(
            command,
            stdout=outputfile.fileno(),
            stderr=errorfile.fileno(), shell=shell,
            cwd=os.path.abspath(path))
        # Wait until it has finished!
        program.wait()

        # Return filename.
        print "Output file name:", outputfile.name
        return outputfile.name
    except:
        import sys
        print 'Error running program.', sys.exc_info()[0]
        return None


def ReadParameter(filename, start, split=None):
    """ Read a parameter from a file.
    Arguments:
    filename -- File to be read from.
    start -- String of beginning of line to read.
    split -- string on which to split the line.

    Returns:
    Parameter retrieved from the file.

    """
    if split is None:
        split = start[-1]

    start = start.strip()
    try:
        for line in open(filename, 'r'):
            if line.lstrip().startswith(start):
                return float(line.split(split)[-1])
    except:
        return 999.
    else:
        return 999.


def CheckProgram(filename, errors):
    """ Check whether a program contained errors.
    Arguments:
    filename -- File to be read from.
    errors -- Keywords that indicate an error.

    Returns:
    Logical, whether program indicates point is physical.

    """
    if filename is None:
        return False

    for line in open(filename, 'r'):
        for word in errors:
            if word in line:
                print "Error in program output."
                print line
                return False
    return True
