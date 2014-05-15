#########################################################################
#                                                                       #
#    M a s s     P l o t                                                #
#                                                                       #
#########################################################################

# SuperPy modules.
import PlotMod as PM
import Stats
import Appearance as AP
import OneDim

# External modules.
import matplotlib.pyplot as plt
import numpy as NP
import sys

# Open the chain with a GUI.
labels, data = PM.OpenData()

# Number of bins.
bins = 50

# Title of plot.
title = "Masses"

# Initalise plot.
fig, ax = PM.NewPlot()
#PM.PlotTicks(AP.xticks,9, ax)
PM.PlotLabels('Particle', 'Mass (GeV)', '')
PM.Appearance()

# Particles to be plotted - label numbers in chain, with 0-based array.
particles = range(12, 44)

# Defensive - make sure we have data for those labels.
if max(particles) > len(labels):
    sys.exit(
    '''No data in chain for those particle numbers.''')

# Alter some default settings.
plt.rcParams['lines.linewidth'] = 12  # Thick lines.
plt.grid(True, which='minor', axis='y')  # Minor grid lines for y-axis.
ax.tick_params(
    axis='x',
    which='minor',
    bottom='off')  # No minor ticks for x-axis.
ax.tick_params(axis='x', which='major', labelsize=10)  # Big x-axis labels.

# Plot information.
dummy = 0  # Counter for x-axis label.
names = []  # List of labels for x-axis.

# Loop of particles to be plotted.
for key in particles:

    # Find relevant quantities from chain.
    name = labels[key]  # Name of particle, mass of which is to be plotted.
    x = data[key]  # Mass to be plotted.
    chisq = data[1]  # Total chi-squared.
    row = NP.argmin(chisq)  # Row number of minimum chi-squared.

    # Calculate statistics.
    profchisq = OneDim.ProfileLike(
        x,
        chisq,
        nbins=bins).profchisq  # Profiled chi-squared.
    bestfit = data[key][row]  # Best-fitting mass.

    xc = OneDim.ProfileLike(x, chisq, nbins=bins).bins  # Bin centers.
    ci = OneDim.ConfidenceIntervals(
        profchisq,
        xc,
        epsilon=AP.epsilon).confint  # Confidence intervals.
    one = ci[1, :]  # One-sigma confidence interval.
    two = ci[0, :]  # Two-sigma confidence interval.

    # Print info to screen - handy.
    print name, bestfit, one, two

    # Make data for x-axis - series of integers.
    xdummy = NP.zeros((one.size)) + dummy  # "Data" for x-axis.

    # Make list of names for x-axis labels.
    # We don't want (GeV) and we don't need m_ either.
    names = names + [name.replace('(GeV)', '').replace('m_', '')]

    # Plot one and two sigma intervals, and bestfit point.
    plt.plot(xdummy, two, '-', color='RoyalBlue', alpha=0.8)
    plt.plot(xdummy, one, '-', color='SeaGreen', alpha=1)
    plt.plot(dummy, bestfit, '*', color='Crimson', ms=15)

    # Increment counter for x-axis "data".
    dummy += 1

# Plot the particle names on the x-axis.
plt.xticks(range(len(names)), names)
locs, xlabels = plt.xticks()
plt.setp(xlabels, rotation=-90)  # Rotate the particle names.
# Might be necessary to make more room.
plt.subplots_adjust(bottom=0.2, left=0.2)

# Proxy for legend - make best-fit, and confidence intervals appear on legend.
plt.plot(-100, 1, '-', color='RoyalBlue', label='$2\sigma$', alpha=0.8)
plt.plot(-100, 1, '-', color='SeaGreen', label='$1\sigma$', alpha=1)
plt.plot(-100, 1, '*', color='Crimson', label='Best-fit', ms=15)

# Set plot limits.
PM.PlotLimits(ax, NP.array([-1, 20, 0, 20000]))

# Make plot.
PM.Legend(title)
PM.SavePlot('Masses')
