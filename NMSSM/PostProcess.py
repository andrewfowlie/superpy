#########################################################################
#                                                                       #
#  P o s t - P r o c e s s                                              #
#                                                                       #
#########################################################################

# Re-calculates likelihood function for a *.txt file.
# But, careful!, doesn't renormalise the posterior weights, though this
# shouldn't be a problem.

# Internal modules.
import Debug as DB  # Debug options.
import Likelihood  # Constraints and likelihood functions.
import PlotMod as PM  # To access this, export PYTHONPATH=$PWD/SuperPlot

# External modules.
import os
import numpy as NP

# Open chain.
chain = os.path.abspath('./example/SB_MO_log_allpost.txt')
cols = len(open(chain, 'rb').readline().split())
rows = len(open(chain, 'rb').readlines())
data = PM.OpenChain(chain)

# Post-process chain.
# Subract 2 for chi-squared and posterior weight.
# + 100 to protect from bugs with the cube not being big enough for the
# likelihood function's cube.
nparams = cols - 2 + 100
ndim = 9  # For CMSSM.

# For row in an existing *.txt file, re-run the likelihood function,
# making alterations to the cube.
for row in range(rows):

    # Initialise the cube to an empty list - it must be a list, not a
    # dictionary.
    cube = [0] * nparams

    # Copy model parameters from the *.txt file to the cube.
    for i in range(ndim):
        # Plus 2 for chi-squared and posterior weight.
        cube[i] = data[i + 2][row]

    # Re-calculate the cube etc with our likelihood function.
    loglike = Likelihood.myloglike(cube, ndim, nparams)

    # Calculate adjusted posterior weight and copy new chi-squared and new cube
    # to post-processed chain.
    pp = {}
    for i in range(len(cube)):
        pp[i + 2] = cube[i]  # We plus 2 for chi-squared and posterior weight.

    # Adjust the posterior weight.
    pp[0] = data[0][row] * NP.exp(- 2 * loglike + 2 * data[1][row])

    # Adjust the chi-squared.
    pp[1] = - 2 * loglike

    # Now write the post-processed cube to a new file.
    name, extension = os.path.splitext(chain)
    ppchain = open(name + '.pp' + extension, 'a')  # Append file.
    for x in pp.itervalues():
        ppchain.write(str(x))  # Data.
        ppchain.write(' ')  # Space.
    # New line.
    ppchain.write('\n')
    ppchain.close()
