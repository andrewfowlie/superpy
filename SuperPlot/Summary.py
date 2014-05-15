#########################################################################
#                                                                       #
#    S u m m a r y                                                      #
#                                                                       #
#########################################################################

# External modules.
import numpy as NP

# SuperPy modules.
import PlotMod as PM
import Stats
import Appearance as AP
import OneDim

# Open the chain with a GUI.
labels, data = PM.OpenData()

# Print information for the parameters.
print 'Param | Best-fit | Posterior Mean | 1 sigma Credible region'
for key, name in labels.iteritems():
    if key == 0 or key == 1 or '\chi^2' in name:
        continue
    x = data[key]
    pw = data[0]
    chisq = data[1]
    bestfit = Stats.BestFit(chisq, x)
    postmean = Stats.PosteriorMean(pw, x)
    pdf = OneDim.PosteriorPDF(
        x,
        pw,
        nbins=AP.nbins,
        bin_limits=AP.bin_limits).pdf
    xc = OneDim.PosteriorPDF(
        x,
        pw,
        nbins=AP.nbins,
        bin_limits=AP.bin_limits).bins
    lowercredibleregion = OneDim.CredibleRegions(
        pdf,
        xc,
        epsilon=AP.epsilon).lowercredibleregion
    uppercredibleregion = OneDim.CredibleRegions(
        pdf,
        xc,
        epsilon=AP.epsilon).uppercredibleregion
    print name, bestfit, postmean, lowercredibleregion[0], uppercredibleregion[0]

# Print best-fit information.
print 'Min ChiSq', data[1].min()
print 'p-value', Stats.PValue(data[1], AP.dof)
