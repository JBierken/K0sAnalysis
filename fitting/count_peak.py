########################################################################################
# Count the number of instances in an invariant mass peak after background subtraction #
########################################################################################
# Note: this could potentially be updated to a new fitting library,
#       to remove all ROOT dependency.


import sys
import os
import numpy as np
from array import array
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
import tools.fittools as ft
import plotting.plotfit as pft

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def count_peak(hist, label, extrainfo, gargs, mode='subtract'):
    ### return the integral (and error estimate) under a peak in a histogram, 
    ### subtracting the background contribution from sidebands
    # input arguments:
    # - hist = ROOT histogram to perform fitting on (left unmodified)
    # - label = only for plotting
    # - extrainfo = only for plotting
    # - gargs = dict containgin global program parameters
    # - mode = 'gfit' or 'subtract'
    #           or 'hybrid', which returns subtract results but makes fancy 'gfit' plot anyway

    # initializations
    nbins = hist.GetNbinsX()
    xlow = hist.GetBinLowEdge(1)
    xhigh = hist.GetBinLowEdge(nbins)+hist.GetBinWidth(nbins)
    fitrange = (xlow, xhigh)
    fitcenter = (xlow + xhigh)/2.
    fithalfwidth = (xhigh - xlow)/4.
    histclone = hist.Clone()
    binlow = histclone.FindBin(fitcenter-fithalfwidth)
    binhigh = histclone.FindBin(fitcenter+fithalfwidth)

    # take into account only sidebands, not peak
    for i in range(binlow, binhigh):
        histclone.SetBinContent(i,0)
        histclone.SetBinError(i,0)

    # make background-only fit
    guess = [0.,0.]
    backfit, paramdict, backfitobj = ft.poly_fit(histclone, fitrange, guess, "Q0")
    if gargs['helpdir'] is not None:
        outfile = os.path.join(gargs['helpdir'], hist.GetName().replace(' ','_')+'_bck.png')
        pft.plot_fit(hist, outfile,
                fitfunc=None, backfit=backfit, label=label, paramdict=paramdict,
                xaxtitle='Invariant mass (GeV)',
                yaxtitle='Number of reconstructed vertices',
                extrainfo=extrainfo, lumi=gargs['lumi'])
    # make signal peak fit if requested
    if(mode=='gfit' or mode=='hybrid'):
        avgbck = paramdict['a0']+0.5*fitcenter*paramdict['a1']
        guess = [
                fitcenter, # peak position
                avgbck*10, # peak 1 height
                (xhigh-xlow)/45., # peak 1 width
                avgbck*10, # peak 2 height
                (xhigh-xlow)/15. # peak 2 width
                ]
        guess += [paramdict['a0'],paramdict['a1']] # background estimate
        globfit, paramdict, globfitobj = ft.poly_plus_doublegauss_fit(hist, fitrange, guess)
        if gargs['helpdir'] is not None:
            outfile = os.path.join(gargs['helpdir'], hist.GetName().replace(' ','_')+'_sig.png')
            pft.plot_fit(hist, outfile,
                    fitfunc=globfit, backfit=backfit, label=label, paramdict=paramdict,
                    xaxtitle='Invariant mass (GeV)',
                    yaxtitle='Number of reconstructed vertices',
                    extrainfo=extrainfo, lumi=gargs['lumi'])

    # METHOD 1: subtract background from peak and count remaining instances
    if(mode=='subtract' or mode=='hybrid'):
        histclone2 = hist.Clone()
        histclone2.Add(backfit,-1)
        npeak = 0.
        npeak_error2 = 0.
        for i in range(binlow+1,binhigh+1):
            npeak += histclone2.GetBinContent(i)
            npeak_error2 += np.power(histclone2.GetBinError(i),2)
        return (npeak,np.sqrt(npeak_error2))

    # METHOD 2: do global fit and determine integral of peak with error
    elif mode=='gfit':
        intpeak = np.sqrt(2*np.pi)*(paramdict['A_{1}']*paramdict['#sigma_{1}']
                                        + paramdict['A_{2}']*paramdict['#sigma_{2}'])
        sidexbinwidth = histclone.GetBinWidth(1) # to do: make more robust and general
        npeak = intpeak/sidexbinwidth # calculate number of instances instead of integral
        npeak_error = globfit.IntegralError(fitrange[0],fitrange[1])
        npeak_error /= sidexbinwidth
        return (npeak,npeak_error)

    else:
        raise Exception('ERROR: peak counting mode not regognized!')


def count_peak_unbinned(values, weights, variable, mode='subtract',
                        label=None, lumi=None, extrainfo=None,
                        histname='sideband', plotdir=None):
    ### wrapper around fitting function using more modern paradigm for input parameters
    # input arguments:
    # - values: np array of sideband values
    # - weights: np array of weights corresponding to values
    # - variable: dict with all information about the sideband variable
    # - mode: passed down to called function
    
    # make a ROOT histogram with the values and weights
    counts = np.histogram(values, variable['bins'], weights=weights)[0]
    errors = np.sqrt(np.histogram(values, variable['bins'], weights=np.power(weights,2))[0])
    hist = ROOT.TH1F(histname, histname, len(variable['bins'])-1, array('f', variable['bins']))
    hist.SetDirectory(0)
    for i, (count, error) in enumerate(zip(counts, errors)):
      hist.SetBinContent(i+1, count)
      hist.SetBinError(i+1, error)

    # make directory to store plots
    if( plotdir is not None and not os.path.exists(plotdir) ): os.makedirs(plotdir)

    # make a dictionary with extra fit info
    fitinfo = {}
    fitinfo['lumi'] = lumi
    fitinfo['helpdir'] = plotdir

    # call underlying function
    return count_peak(hist, label, extrainfo, fitinfo, mode=mode)
