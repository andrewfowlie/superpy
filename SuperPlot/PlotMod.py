##########################################################################
#                                                                       #
#     P l o t   M o d                                                   #
#                                                                       #
#########################################################################

# General functions for plotting data, defined once so that they can be used/edited
# in a consistent manner.

# External modules.
from pylab import *
import os
import gobject
import pygtk
pygtk.require('2.0')
import gtk
import sys
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import NullLocator
import matplotlib.pyplot as plt
import re
import numpy as NP
import cPickle as pickle

# SuperPy modules.
import Appearance as AP
import SLHADictionary as SLHA


def OpenData():
    """ Wrapper for opening and labelling a chain.

    Returns:
    data -- Dictionary with chain.
    label -- Dictionary with chain's labels.

    """
    # Select file with a GUI.
    filename = OpenFileGUI()

    if filename is None:
        sys.exit("No data file was selected - did you press cancel?")

    # Open file.
    data = OpenChain(filename)
    # Assign labels.
    label = LabelChain(data)

    return label, data


def OpenFileGUI():
    """ GUI for opening a file with a file browser.

    Return:
    filename -- Name of file selected with GUI.

    """
    # Select the file from a dialog box.
    dialog = gtk.FileChooserDialog("Open..",
                                   None,
                                   gtk.FILE_CHOOSER_ACTION_OPEN,
                                   (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                    gtk.STOCK_OPEN, gtk.RESPONSE_OK))
    dialog.set_default_response(gtk.RESPONSE_OK)

    # Only show particular files.
    filter = gtk.FileFilter()
    filter.set_name("Text files, serial data or info file.")
    filter.add_pattern("*.txt")
    filter.add_pattern("*.pkl")
    filter.add_pattern("*.info")
    dialog.add_filter(filter)
    filter = gtk.FileFilter()
    dialog.add_filter(filter)

    response = dialog.run()
    if response == gtk.RESPONSE_OK:
        print 'File:', dialog.get_filename(), 'selected.'
    elif response == gtk.RESPONSE_CANCEL:
        print 'Error: no file selected.'

    # Save the file name/path
    filename = dialog.get_filename()
    dialog.destroy()

    return filename


def SaveFileGUI():
    """ GUI for saving a file with a file browser.

    Return:
    filename -- Name of file selected with GUI.

    """
    # Select the file from a dialog box.
    dialog = gtk.FileChooserDialog("Open..",
                                   None,
                                   gtk.FILE_CHOOSER_ACTION_SAVE,
                                   (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                    gtk.STOCK_OPEN, gtk.RESPONSE_OK))
    dialog.set_default_response(gtk.RESPONSE_OK)

    # Only show particular files.
    filter = gtk.FileFilter()
    filter.set_name("PDF, PS, EPS or PNG files.")
    filter.add_pattern("*.pdf")
    filter.add_pattern("*.eps")
    filter.add_pattern("*.ps")
    filter.add_pattern("*.png")
    dialog.add_filter(filter)
    filter = gtk.FileFilter()
    dialog.add_filter(filter)

    response = dialog.run()
    if response == gtk.RESPONSE_OK:
        print 'File:', dialog.get_filename(), 'selected.'
    elif response == gtk.RESPONSE_CANCEL:
        print 'Error: no file selected.'

    # Save the file name/path
    filename = dialog.get_filename()
    dialog.destroy()

    return filename


def OpenChain(filename):
    """ Open a text file or serialised data. If opening a text file, serialise data for future plots.

    Arguments:
    filename -- Name of chain file.

    Returns:
    data -- Dictionary of chain indexed with integers.

    """
    # If serialised data, open it. Else open a data file and serialise data too,
    # with the same name but with the sl suffix.
    name, extension = os.path.splitext(filename)

    if 'pkl' in extension:
        # Open the data from serialised data.
        print 'Opening chain...'
        serialf = open(filename, 'rb')
        data = pickle.load(serialf)
        serialf.close()

    else:
        # Open the *txt file.
        data = OpenDataFile(filename)

        # Serialize the data for quicker future reading.
        print 'Dumping chain...'
        serialf = open(name + '.pkl', 'wb')
        # Try to catch memory errors.
        try:
            pickle.dump(data, serialf)
        except:
            memoryerror
        serialf.close()

    print 'Success: returning chain.'
    return data


def OpenDataFile(filename):
    """ Open a text file and return dictonary of numpy arrays of data.

    Arguments:
    filename -- Name of chain file.

    Returns:
    data -- Dictionary of chain indexed with integers.

    """
    # Find size of text file.

    # Find number of columns
    cols = len(open(filename, 'rb').readline().split())

    # Try quick approximate method that might fail for number of rows.

    # Find total length of file with seek.
    dataf = open(filename, 'rb')
    dataf.seek(0, 2)  # Go to final line.
    datalen = dataf.tell()  # Read position.
    dataf.close()

    # Find length of a single row.
    rowlen = len(open(filename, 'rb').readline())

    # If rows are equal in length, rows = datalen / rowlen.
    if datalen % rowlen == 0:
        # Quick method.
        rows = datalen / rowlen
    else:
        # Slow method.
        rows = len(open(filename, 'rb').readlines())

    # Initialise data as dictionary of arrays.
    data = {}
    for key in range(cols):
        data[key] = NP.zeros(rows)

    # Populate data - NB we convert from string to float.
    print 'Opening chain...'
    row = 0
    for line in open(filename, 'rb'):
        for key, word in enumerate(line.split()):
            data[key][row] = float(word)
        row += 1

    return data


def LabelChain(data):
    """ Match labels to the data.
    i) Load an *info file.
    ii) If that doesn't match, use data indicies.

    Arguments:
    data -- Data chain, to match arguments with.

    Return:
    label -- Dictonary of labels, indexed with the same
    indicies as the data chain.

    """
    label = {}

    # Attempt to load a SuperBayeS style info file.
    print "Try an info file."
    label = ReadInfo()
    if len(data.keys()) != len(label):
        print "Warning: length of info file did not match chain."
        print "Missing labels are integers."
        # Label all unlabelled columns with integers.
        for key in data.keys():
            if label.get(key) is None:
                label[key] = str(key)  # NB convert integer index to string

    return label


def ReadInfo():
    """ Read labels from a SuperBayeS style info file.
    The file is chosen with a GUI.
    """
    label = {}

    # Open the info file with a GUI.
    filename = OpenFileGUI()

    if filename is None:
        return label

    info = open(filename, 'rb')
    for line in info:
        # Look for labX=string.
        if line.startswith('lab'):
            string = line.strip().split('=')[::-1][0]
            index = int(
                re.findall(
                    re.escape('lab') +
                    "(.*)" +
                    re.escape('='),
                    line)[0])

            # Correct index - SuperBayeS info files begin at 1, and ignore
            # posterior weight and ChiSq in positions 0 and 1 of the chain.
            index += 1

            # Replace with dictionary of names if possible.
            if SLHA.name.get(string) is not None:
                string = SLHA.name[string]

            label[index] = string

    info.close()

    # Add posterior weight and ChiSq.
    label[0] = '$p_i$'
    label[1] = '$\chi^2$'

    return label


def PlotData(x, y, Scheme):
    """ Plot a point with a particular color scheme.

    Arguments:
    x -- Data to be plotted on x-axis.
    y -- Data to be plotted on y-axis.
    Scheme -- Object containing plot options.

    """
    plt.plot(
        x,
        y,
        Scheme.Symbol,
        color=Scheme.Colour,
        label=Scheme.Label,
        ms=Scheme.Size)


def NewPlot():
    """ Clear the plot so we start afresh."""
    fig = plt.figure(figsize=AP.size)  # Size in inches.
    ax = fig.add_subplot(1, 1, 1)
    return fig, ax


def SavePlot(name):
    """ Save the plot, with a descriptive name.

    Arguments:
    name -- Prefix of filename, without extension.

    """
    plt.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)


def Appearance():
    """ Specify the plots appearance, with e.g. font types etc.
    """

    # Make the lines thicker.
    plt.rcParams['lines.linewidth'] = 4

    # Plot gridlines over the plot - can be useful
    plt.grid(True)

    # Set the fonts - I like mathpazo/Palatino, but
    # that would introduce an unneccesasry dependency.
    rc('text', usetex=True)
    rc('font',
       **{'family': 'serif',
          'serif': ['Computer Modern Roman'],
          'size': '20'})

    # Set size of plot in inches.
    plt.rcParams['figure.figsize'] = [2.5, 2.5]


def Legend(title=None):
    """ Turn on the legend.

    Arguments:
    title -- Title of legend.

    """
    leg = plt.legend(prop={'size': 16}, shadow=False, fancybox=True,
                     title=title, loc='best', borderaxespad=1.,
                     scatterpoints=1, numpoints=1)
    leg.get_frame().set_alpha(0.5)


def PlotLimits(ax, plot_limits=None):
    """ If specified plot limits, set them.

    Arguments:
    ax -- Axis object.
    plot_limits -- Array of plot limits in order xmin, xmax, ymin, ymax.

    """
    if plot_limits is not None:
        ax.set_xlim([plot_limits[0], plot_limits[1]])
        ax.set_ylim([plot_limits[2], plot_limits[3]])


def PlotTicks(xticks, yticks, ax):
    """ Set the numbers of ticks on the axis.

    Arguments:
    ax -- Axis object.
    xticks -- Number of required major x ticks.
    yticks -- Number of required major y ticks.

    """
    # Set major x, y ticks.
    ax.xaxis.set_major_locator(MaxNLocator(xticks))
    ax.yaxis.set_major_locator(MaxNLocator(yticks))
    # Auto minor x and y ticks.
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


def PlotLabels(xlabel, ylabel, plottitle=''):
    """ Plot axis labels.

    Arguments:
    xlabel -- Label for x-axis.
    ylabel -- Label for y-axis.
    plottile -- Title appearing above plot.

    """
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # Plot title.
    plt.title(plottitle)


def PlotImage(xdata, ydata, data, Scheme, zlabel='', bin_limits=None):
    """ Plot data as an image.

    Arguments:
    xdata -- x-axis data.
    ydata -- y-axis data.
    data -- Image data, i.e., z-axis data.
    Scheme -- Object containing appearance options, colours etc.
    zlabel -- Label for colour bar.
    bin_limits -- Bin limits.

    """

    # Set the aspect so that resulting figure is a square.
    extent = NP.zeros((4))
    extent[0] = min(xdata)
    extent[1] = max(xdata)
    extent[2] = min(ydata)
    extent[3] = max(ydata)

    if bin_limits is None:
        bin_limits = extent
    else:
        bin_limits = NP.array(
            (bin_limits[0][0],
             bin_limits[0][1],
             bin_limits[1][0],
             bin_limits[1][1]))

    aspect = (bin_limits[1] - bin_limits[0]) / (bin_limits[3] - bin_limits[2])

    # Interpolating perhaps misleads, if you don't want it set
    # interpolation='nearest'. NB that imshow is annoying - it reads y,x
    # rather than x,y so we take transpose.
    plt.im = plt.imshow(data.T, cmap=Scheme.ColourMap, extent=bin_limits,
                        interpolation='bilinear', label=Scheme.Label,
                        origin='lower', aspect=aspect)
    # Plot a colour bar.
    cb = plt.colorbar(plt.im, orientation='horizontal', shrink=0.5)
    # Colour bar label.
    cb.ax.set_xlabel(zlabel)


def PlotContour(xdata, ydata, data, levels, names, Scheme, bin_limits=None):
    """ Make unfilled contours for a plot.
    Arguments:
    xdata -- x-axis data.
    ydata -- y-axis data.
    data -- Data to be contoured.
    levels -- Levels at which to draw contours.
    names -- Labels for the contour levels.
    Scheme -- Object containing appearance options, colours etc.
    bin_limits -- Bin limits.
    """

    # Set the bin limits correctly.
    extent = NP.zeros((4))
    extent[0] = min(xdata)
    extent[1] = max(xdata)
    extent[2] = min(ydata)
    extent[3] = max(ydata)

    if bin_limits is None:
        bin_limits = extent
    else:
        bin_limits = NP.array(
            (bin_limits[0][0],
             bin_limits[0][1],
             bin_limits[1][0],
             bin_limits[1][1]))

    # Make the contours of the levels.
    cset = plt.contour(
        data.T,
        levels,
        linewidths=2,
        colors=Scheme.Colour,
        hold='on',
        extent=bin_limits,
        interpolation='bilinear',
        origin='lower',
        linestyles=[
            '--',
            '-'])

    # Set the contour labels - they will show labels.
    fmt = {}
    for i, s in zip(cset.levels, names):
        fmt[i] = s
    plt.clabel(cset, inline=True, fmt=fmt, fontsize=12, hold='on')


def PlotFilledContour(
        xdata,
        ydata,
        data,
        levels,
        names,
        Scheme,
        bin_limits=None):
    """ Make filled contours for a plot.

    Arguments:
    xdata -- x-axis data.
    ydata -- y-axis data.
    data -- Data to be contoured.
    levels -- Levels at which to draw contours.
    names -- Labels for the contour levels.
    Scheme -- Object containing appearance options, colours etc.
    bin_limits -- Bin limits.

    """
    # Set the aspect so that resulting figure is a square.
    extent = NP.zeros((4))
    extent[0] = min(xdata)
    extent[1] = max(xdata)
    extent[2] = min(ydata)
    extent[3] = max(ydata)

    # Find the bin limits, and if neccessary change the format.
    if bin_limits is None:
        bin_limits = extent
    else:
        bin_limits = NP.array(
            (bin_limits[0][0],
             bin_limits[0][1],
             bin_limits[1][0],
             bin_limits[1][1]))

    # We need to ensure levels are in ascending order, and append the list with one.
    # This makes 2 intervals (between 3 values) that will be shown with
    # colours.
    levels = NP.append(levels, 1.0)

    # Filled contours.
    cset = plt.contourf(data.T, levels,
                        colors=Scheme.Colours,
                        hold='on', extent=bin_limits,
                        interpolation='bilinear', origin='lower',
                        alpha=0.7)

    # Plot a proxy for the legend - plot spurious data outside plot limits,
    # with legend entry matching colours of filled contours.
    for i, value in enumerate(Scheme.Colours):
        plt.plot(-1.5 * abs(min(xdata)),
                 1.5 * abs(max(ydata)),
                 's',
                 color=Scheme.Colours[i],
                 label=AP.LevelNames[i],
                 alpha=0.7,
                 ms=15)


def PlotBand(x, y, width, ax):
    """ Plot a band around a line.
    This is for a theoretical error. We find the largest and smallest
    y within +/width of the value of x, and fill between these largest and smallest
    x and y.

    Arguments:
    x -- x-data to be plotted.
    y -- y-data to be plotted.
    width -- Width of band - it is this width on the left and right hand-side.
    ax -- An axis object to plot teh band on.

    """
    # Find upper line, and lower line of the shifted data.
    uy = NP.zeros(len(y))
    ly = NP.zeros(len(y)) + 1e90
    for i in range(len(x)):
        for j in range(len(x)):
            # Find lowest/highest point within width of that point.
            if abs(x[i] - x[j]) < width:
                if y[j] < ly[i]:
                    ly[i] = y[j]
                if y[j] > uy[i]:
                    uy[i] = y[j]

    # Finally plot.
    ax.fill_between(x, ly, uy, where=None, facecolor=AP.TauColour, alpha=0.7)
    # Proxy for legend.
    plt.plot(-1, -1, 's', color=AP.TauColour,
             label=AP.TauLabel, alpha=0.7, ms=15)


def PlotLine(
        filename,
        colour="SteelBlue",
        width=2.5,
        style='--',
        label=None,
        xc=0,
        yc=1):
    """ Plot a line on top an existing plot.

    Arguments:
    filename -- Name of data file to be plotted.
    colour -- Colour of line.
    width -- Width of line.
    style -- Style of line.
    label -- Label with which to annotate line.
    xc -- Column of x-data.
    yc -- Column of y-data.

    """

    # Unpack data from file.
    data = NP.genfromtxt(filename, unpack=True)

    # Order points by their x-value
    order = data[xc].argsort()
    x = data[xc][order]
    y = data[yc][order]

    # Plot data.
    plt.plot(x, y, style, alpha=0.8, c=colour, linewidth=width)
    if label is not None:
        plt.annotate(
            label, xy=(
                x[0], y[0]), xytext=(
                5, 5), textcoords='offset points', va='center', arrowprops=None)


def PlotPoints(
        filename,
        colour="SteelBlue",
        size=2.5,
        style='x',
        label=None,
        xc=0,
        yc=1):
    """ Plot scatter points on top an existing plot.

    Arguments:
    filename -- Name of data file to be plotted.
    colour -- Colour of points.
    size -- Size of points.
    style -- Style of points.
    label -- Label with which to annotate line.
    xc -- Column of x-data.
    yc -- Column of y-data.

    """

    data = NP.genfromtxt(filename, unpack=True)
    plt.plot(
        data[xc][:],
        data[yc][:],
        style,
        alpha=0.8,
        c=colour,
        ms=size,
        label=label)
