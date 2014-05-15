#########################################################################
#                                                                       #
#    P u l l / c h i - s q u a r e d     P l o t                        #
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
import math

# Open the chain with a GUI.
labels, data = PM.OpenData()

# Find minimum chi-squared.
chisq = data[1]
row = NP.argmin(chisq)
mchisq = chisq[row]
print "Min chi-squared is,", mchisq

# Labels for plot.
chilabel = 'Total'
title = ''

# Chi-squared contributions in chain.
vars = range(10, 20)

# Defensive - make sure we have data for those labels.
if max(vars) > len(labels):
    sys.exit(
        '''No data in chain for those chi-squared numbers.''')

# Initialise plot.
fig, ax = PM.NewPlot()
PM.PlotTicks(10, len(vars), ax)
PM.PlotLabels('$\chi^2$', 'Contribution', '')
PM.Appearance()

# Alter some default settings.
# plt.rcParams['lines.linewidth'] = 18 # Thick lines.
plt.grid(True, which='minor', axis='x')  # Minor grid lines for x-axis.
ax.tick_params(
    axis='y',
    which='minor',
    left='off')  # No minor ticks for y-axis.

# Plot limit for chi-squared
limchisq = int(math.floor(mchisq)) + 3  # Minimum chi-squared + padding.

# Set the plot ticks for chi-squared axis.
ax.set_xticks(range(limchisq))  # 0, 1..., ~ minimum chi-squared.
ax.set_xlim([0, limchisq])  # 0 - minimum chi-squared + padding.
# Number of items, with padding i.e. -1 and +1.
ax.set_ylim([-1, len(vars) + 1])

# List of axis labels.
names = []

# Loop over the experimental measurements.
for i, var in enumerate(vars):

    # Plot bar of chi-squared.
    ax.plot([0, data[var][row]], [i, i], '-',
            color='FireBrick', alpha=0.8, solid_capstyle="butt")

    # Append label.
    names.append(labels[var])

# Plot vertical line for total chi-squared.
plt.rcParams['lines.linewidth'] = 4  # Regular lines.
plt.plot([chisq[row],
          chisq[row]],
         [-1,
          len(vars) + 1],
         '--',
         color='FireBrick',
         alpha=0.8,
         label=chilabel,
         solid_capstyle="butt")

# Plot labels.
plt.yticks(range(len(names)), names)
# Might be necessary to make extra space.
plt.subplots_adjust(left=0.4, bottom=0.1)

# Make plot.
PM.Legend(title)
PM.SavePlot('ChiSq')
