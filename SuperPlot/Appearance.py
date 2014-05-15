#########################################################################
#                                                                       #
#    Set colour schemes etc here.                                       #
#                                                                       #
#########################################################################

# External modules.
import numpy as NP
from pylab import get_cmap


class Scheme:

    """ Holds information for how a piece of data should be plotted. """

    def __init__(
            self,
            Colour=None,
            Symbol=None,
            Label=None,
            ColourMap=None,
            Size=5,
            plot_limits=None,
            Colours=None):
        """
        Define an appearance scheme.

        Arguments:
        Colour -- Colour for a line/point.
        Symbol -- Indicates point style e.g. cirlce 'o' or line style e.g '--'.
        Label -- Label for legend.
        ColourMap -- Colour map for 2D plots.
        Size -- Size of points.
        plot_limits -- Axes limits.
        Colours -- List of colours to be iterated, for, e.g., filled contours.

        """
        self.Colour = Colour
        self.Symbol = Symbol
        self.Label = Label
        self.ColourMap = ColourMap
        self.Size = Size
        self.plot_limits = plot_limits
        self.Colours = Colours

#########################################################################

# Schemes for various data types.

PosteriorMean = Scheme(Colour='SeaGreen', Symbol='o', Label='Posterior Mean')

BestFit = Scheme(Colour='Brown', Symbol='*', Label=r'Best-fit point', Size=12)

Posterior = Scheme(
    Colour='RoyalBlue',
    Symbol='-',
    Label=r'Posterior pdf',
    ColourMap=get_cmap('GnBu'),
    Colours=[
        'RoyalBlue',
        'SeaGreen'])

ProfLike = Scheme(
    Colour='DarkOrange',
    Symbol='--',
    Label=r'Profile likelihood',
    ColourMap=get_cmap('Reds'),
    Colours=[
        'DarkOrange',
        'Brown'])

ProfChiSq = Scheme(
    Colour='DarkOrange',
    Symbol='--',
    Label=r'$\Delta \chi^2$',
    Colours=[
        'Gold',
        'Peru'])

Scatter = Scheme(Symbol='o', ColourMap=get_cmap('Reds'), Size=15)

ConfInterval = {}
ConfInterval[0] = Scheme(
    Colour='DarkOrange',
    Symbol='o',
    Label=r'$2\sigma$ confidence interval')
ConfInterval[1] = Scheme(
    Colour='Brown',
    Symbol='o',
    Label=r'$1\sigma$ confidence interval')

CredibleRegion = {}
CredibleRegion[0] = Scheme(
    Colour='SeaGreen',
    Symbol='-',
    Label=r'$2\sigma$ credible region')
CredibleRegion[1] = Scheme(
    Colour='RoyalBlue',
    Symbol='-',
    Label=r'$1\sigma$ credible region')

#########################################################################

# Technical plot options.

dof = 10
epsilon = NP.array([0.05, 0.32])  # Values of alpha, in ascending order.
plot_limits = None  # NP.array((0,0.2,0,1.3))
bin_limits = None  # [[0,1000],[0,1000]]
nbins = 70

#########################################################################

# Plot appearance.

xticks = 5  # Numbers of ticks.
yticks = 5
# Labels for two-dimensional regions.
LevelNames = ['$2\sigma$ region', '$1\sigma$ region']
ChiSqLevelNames = [
    '$2\sigma$ excluded',
    '$1\sigma$ excluded']  # Labels for delta chi-squared regions.
plottitle = 'Fowlie (2014)'
PLTitle = 'PL'  # Legend titles for plots.
PDFTitle = 'PDF'
ScatterTitle = 'Scatter'
OneDimTitle = 'PL and PDF'
ChiSqTitle = None

# Size in inches.
size = (8, 8)

# For theoretical error on delta chi-squared plots.
Tau = 2
TauColour = 'RoyalBlue'
TauLabel = 'Theory error'
