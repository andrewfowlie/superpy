#########################################################################
#                                                                       #
# L i n e - P r o c e s s
#                                                                       #
#########################################################################

# Re-calculates likelihood function for single line in a *.txt file.
# You supply the linenumber of interest at the command line.
# Careful! - first line is zeroth line, which might differ from awk and
# Fortran etc.

# Internal modules.
import Debug as DB  # Debug options.
import Likelihood  # Constraints and likelihood functions.
import PlotMod as PM  # To access this, export PYTHONPATH=$PWD/SuperPlot

# External modules.
import sys

# Check that line number was supplied.
if len(sys.argv) != 2:
    sys.exit('You should supply line number.')

# Open the chain with a GUI.
labels, data = PM.OpenData()
ndim = 10  # Number of model parameters - 10 for CNMSSM.
nparams = 100  # Approximate number of items in cube.

# Read line to re-calculate.
# Plust one to ignore first argument is name of file, e.g. python
# LineProcess.py arg1 etc
linenumber = int(sys.argv[1])  # We need to convert from string to integer.

# Initialise the cube to an empty list - it must be a list, not a dictionary.
cube = [0] * nparams

# Copy model parameters from the *.txt file to the cube.
for i in range(ndim):
    # Plus 2 for chi-squared and posterior weight.
    cube[i] = data[i + 2][linenumber]

# Re-calculate the cube etc with this likelihood function.
loglike = Likelihood.myloglike(cube, ndim, nparams)
