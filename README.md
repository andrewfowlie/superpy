# SuperPy
*SuperPy* is a Python code that scans the CMSSM or CNMSSM's parameter space to find the regions that best agree with experimental data with Bayesian statistics.

1. Draws a model from the CMSSM or CNMSSM with your priors, with MultiNest.
2. Calculates predictions for the model with SoftSUSY, SuperISO, microMEGAs, FeynHiggs, SUSY_FLAVOR, HiggsSignals, fastlim etc.
3. Calculates the chi-squared by comparing these predictions with experimental data.
4. Returns this chi-squared to MultiNest, which explores the parameter space with an efficient algorithm.

With SuperPy, it is **easy** to modify the likelihoods, priors or model.

# SuperPlot
*SuperPlot* is a GUI that plots SuperPy, SuperBayeS (with its information file format) or generally MultiNest results. It can calculate and plot:
* One- and two-dimensional marginalised posterior pdf and credible regions.
* One- and two-dimensional marginalised profile likelihood and confidence intervals.
* Best-fit points.
* Posterior means.
* Three-dimensional scatter plots.

# Installation
From within the SuperPy directory,

      make all

or, for only the SuperPlot routines,

      make python

You might need to install the matplotlib Python plotting library.

# Running SuperPy
From within the */SuperPy* directory,

     python SuperPy.py

Alter the settings in
* Likelihoods in *Likelihood.py*
* Priors in *Prior.py*
* Scanning options in *MNOptions.py*

# Running SuperPlot
From within the *SuperPy/SuperPlot* sub-directory,

     python SuperGUI.py

A GUI window will appear, to select a chain. Select e.g. the *.txt* file in the */SuperPy/examples* sub-directory. A second GUI window will appear to select an information file. Select e.g. the *.info file in the */examples* sub-directory. Finally, select the variables and the plot type in the resulting GUI, and click *Make Plot*.

# Altered programs
My changes are marked "SuperPy". I made minor modifications to:

## SUSY_FLAVOR
The default behaviour is to read a fixed file and write to a
particular file. I've changed it to read a filename at the command line and
write to the screen:
* ./susy_flavor_2.5.1/lib/sflav_io.f (altered input/output)
* ./susy_flavor_2.5.1/susy_flavor_prog.f (read SLHA filename from command line)

## Fastlim
Fastlim doesn't like being imported from other directories, because its
file structure isn't relative.
* ./fastlim-1.0/read_data.py (made file names relative)

I also fixed a little bug that it didn't
like reading NMSSM mass spectra because it didn't understand 5 neutralinos etc.
* ./fastlim-1.0/basic_func.py (support NMSSM spectra)

## HiggsSignals
My gripe was that file operations, e.g. call `system('rm -f HS_analyses.txt')` are
present inside the Fortran codes. I removed them, and perform those
operations once only at build (see my makefile).
* ./HiggsSignals-1.2.0/datatable.f90 (removed file operations)
* ./HiggsSignals-1.2.0/HiggsSignals_subroutine.f90 (removed file operations)

I fixed some compilation/library linking issues:
* ./HiggsSignals-1.2.0/configure (set library paths manually)
* ./HiggsSignals-1.2.0/makefile.in (linked my SuperPy program)

I added my program for calling the HiggsSignals functions:
* ./HiggsSignals-1.2.0/example_programs/SuperPy.f90

## MicrOMEGAS
Small changes, so that the file read input files in SLHA format by default.
* ./micromegas_3.6.9.2/MSSM/main.c (read SLHA filename from command line)
* ./micromegas_3.6.9.2/NMSSM/main.c (read SLHA filename from command line)

## SUSY-HIT
Out of the box not gfortran compliant, I fixed that
* ./susyhit/*f (made gfortran compliant)
* ./susyhit/makefile (replaced g77 with gfortran)

Also, made it read SLHA file name from command line
* ./susyhit/sdcay.f  (read SLHA filename from comand line)

## SOFTSUSY
Wanted to add Jacobian for naturalness priors, with derivatives
of (MZ, tanb) wrt (mu, b). Some of the code was already present.
* ./softsusy-3.4.0/src/softpoint.cpp (added fine-tuning functions)
* ./softsusy-3.4.0/src/softsusy.cpp (print derivatives for naturalness priors)

## NMSSMTools
Write fine-tuning parameterswq
* ./NMSSMTools_4.2.1/sources/NS_output.f

Added program to call NMSSM functions and build that code
* ./NMSSMTools_4.2.1/main/pyspec.f
* ./NMSSMTools_4.2.1/main/makefile

## FeynHiggs:
My program to call relevant FeynHiggs libraries:
* ./FeynHiggs-2.10.0/example/SuperPy.F

Fixed a couple of bugs in the SLHA writer, particuarly regarding
NMSSM mass spectra
* ./FeynHiggs-2.10.0/src/SLHA/SLHAWrite.f