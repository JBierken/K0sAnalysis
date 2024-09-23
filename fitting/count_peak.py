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

# -------------------------------------------------------------------------------------
# Determine FWHM of a given Histogram:
# -------------------------------------------------------------------------------------
def Confidence(hist, confLevel = 0.95):
    ix      = hist.GetXaxis().FindBin(hist.GetMean())
    ixlow   = ix
    ixhigh  = ix
    nb      = hist.GetNbinsX()
    ntot    = hist.Integral()
    nsum    = hist.GetBinContent(ix)
    width   = hist.GetBinWidth(ix)
    
    if ntot==0:                                                 return 0
    while (nsum < confLevel * ntot):
        nlow    = hist.GetBinContent(ixlow-1)   if ixlow>0      else 0
        nhigh   = hist.GetBinContent(ixhigh+1)  if ixhigh<nb    else 0
        if (nsum+max(nlow,nhigh) < confLevel * ntot):
            if (nlow>=nhigh and ixlow>0):
                nsum    += nlow
                ixlow   -=1
                width   += hist.GetBinWidth(ixlow)
            elif ixhigh<nb:
                nsum    += nhigh
                ixhigh  +=1
                width   += hist.GetBinWidth(ixhigh)
            else:       raise ValueError('BOOM')
        else:
            if (nlow>nhigh):
                width   += hist.GetBinWidth(ixlow-1) * (confLevel * ntot - nsum) / nlow
            else:
                width   += hist.GetBinWidth(ixhigh+1) * (confLevel * ntot - nsum) / nhigh
            nsum        = ntot
    
    return width

# -------------------------------------------------------------------------------------
# Make background/signal fit and calculate event count:
#
#   Input arguments:  - hist = ROOT histogram to perform fitting on (left unmodified)
#                     - label = only for plotting
#                     - extrainfo = only for plotting
#                     - gargs = dict containgin global program parameters
#                     - mode = 'gfit' or 'subtract'
#                         or 'hybrid', which returns subtract results but makes fancy 'gfit' plot anyway
#   Return:           - return the integral (and error estimate) under a peak in a histogram, 
#                         subtracting the background contribution from sidebands
# -------------------------------------------------------------------------------------
def count_peak(hist, label, extrainfo, gargs, mode='subtract'):
    # initializations
    nbins           = hist.GetNbinsX()
    xlow            = hist.GetBinLowEdge(1)
    xhigh           = hist.GetBinLowEdge(nbins)+hist.GetBinWidth(nbins)
    fitrange        = (xlow, xhigh)
    fitcenter       = (xlow + xhigh)/2.     # Center of fit range
    fithalfwidth    = (xhigh - xlow)/4.     # differnce from center
    histclone       = hist.Clone()
    binlow          = histclone.FindBin(fitcenter-fithalfwidth)
    binhigh         = histclone.FindBin(fitcenter+fithalfwidth)

    # take into account only sidebands, not peak
    for i in range(binlow, binhigh):
        histclone.SetBinContent(i,0)
        histclone.SetBinError(i,0)

    # make background-only fit
    guess = [0.,0.]
    backfit, paramdict, backfitobj = ft.poly_fit(histclone, fitrange, guess, "WLQ0")
    if gargs['helpdir'] is not None:
        lumitext    = '' if gargs['lumi'] is None else '{0:.3g} '.format(float(gargs['lumi'])/1000.) + 'fb^{-1} (13 TeV)'
        outfile     = os.path.join(gargs['helpdir'], hist.GetName().replace(' ','_')+'_bck.png')
        pft.plot_fit(hist, outfile,
                fitfunc     =None, 
                backfit     =backfit, 
                label       =label, 
                paramdict   =paramdict,
                xaxtitle    ='Invariant mass (GeV)',
                yaxtitle    ='Number of reconstructed vertices',
                extrainfo   =extrainfo, 
                lumitext    =lumitext
        )

    # make signal peak fit if requested
    if(mode=='gfit' or mode=='hybrid'):
        # If more than 50 events --> use Double Gauss
        # Else                   --> use Single Gauss
        if(hist.GetEffectiveEntries() <= 100):
            guess = [
                    fitcenter,                              # peak position
                    hist.GetMaximum()/2,                    # peak 1 height
                    hist.GetRMS()                           # peak 1 width
            ]
            guess += [paramdict['a0'],paramdict['a1']]     # background estimate
            #guess += [0., 0.]                              # background estimate
            
            globfit, paramdict, globfitobji, globres = ft.poly_plus_gauss_fit(hist, fitrange, guess)

        else:
            guess = [
                    fitcenter,                              # peak position
                    hist.GetMaximum()/2,                    # peak 1 height
                    hist.GetRMS()/4,                        # peak 1 width
                    hist.GetMaximum()/2,                    # peak 2 height
                    hist.GetRMS()                           # peak 2 width
            ]
            guess += [paramdict['a0'],paramdict['a1']]     # background estimate
            #guess += [0., 0.]                               # background estimate
           
            globfit, paramdict, globfitobj, globres = ft.poly_plus_doublegauss_fit(hist, fitrange, guess)
        
        # separate background component of fit
        backfit2 = ROOT.TF1("fitfunc","pol1(0)", fitrange[0], fitrange[1])
        backfit2.SetParameters(paramdict['a0'], paramdict['a1'])
        
        if gargs['helpdir'] is not None:
            outfile = os.path.join(gargs['helpdir'], hist.GetName().replace(' ','_')+'_sig.png')
            pft.plot_fit(hist, outfile,
                    fitfunc     =globfit, 
                    backfit     =backfit2, 
                    label       =label, 
                    paramdict   =paramdict,
                    xaxtitle    ='Invariant mass (GeV)',
                    yaxtitle    ='Reconstructed vertices',
                    extrainfo   =extrainfo, 
                    lumitext    =lumitext,
                    sideband    =[fitcenter-fithalfwidth, fitcenter+fithalfwidth]
            )
    
    # ---------------------------------------------------------------------------------
    # Make single Gauss fit for L_xy confidence study:
    guess = [
            fitcenter,                              # peak position
            hist.GetMaximum()/2,                    # peak 1 height
            hist.GetRMS()                           # peak 1 width
    ]
    guess += [paramdict['a0'],paramdict['a1']]      # background estimate
    #guess += [0., 0.]                               # background estimate
           
    globfit2, paramdict2, globfitobj2, globres2 = ft.poly_plus_gauss_fit(hist, fitrange, guess)
    
    # separate background component of fit
    backfit3 = ROOT.TF1("fitfunc","pol1(0)", fitrange[0], fitrange[1])
    backfit3.SetParameters(paramdict2['a0'], paramdict2['a1'])
    
    if gargs['helpdir2'] is not None:
            outfile2 = os.path.join(gargs['helpdir2'], hist.GetName().replace(' ','_')+'_sig.png')
            pft.plot_fit(hist, outfile2,
                    fitfunc     =globfit2, 
                    backfit     =backfit3, 
                    label       =label, 
                    paramdict   =paramdict2,
                    xaxtitle    ='Invariant mass (GeV)',
                    yaxtitle    ='Reconstructed vertices',
                    extrainfo   =extrainfo, 
                    lumitext    =lumitext,
                    sideband    =[fitcenter-fithalfwidth, fitcenter+fithalfwidth]
            )
    # ---------------------------------------------------------------------------------
    
    # METHOD 1: subtract background from peak and count remaining instances
    if(mode=='subtract' or mode=='hybrid'):
                
        # Calculate the event count (integral) and statistical error
        if hist.GetEffectiveEntries()                           > 20:                               # Calculate integral
            if hist.GetEffectiveEntries()                       <= 100: 
                sigfit                                          = ROOT.TF1("signal","gaus(0)")
                sigfit.SetParameters(paramdict[r'A'], paramdict[r'#mu'], paramdict[r'#sigma'])
            else:                                
                sigfit                                          = ROOT.TF1("signal","gaus(0)+gaus(3)")
                sigfit.SetParameters(paramdict[r'A_{1}'], paramdict[r'#mu'], paramdict[r'#sigma_{1}'],
                                                                paramdict[r'A_{2}'], paramdict[r'#mu'], paramdict[r'#sigma_{2}'])
            
            # Determine general fit parameters required for integration 
            params                                              = globres.GetParams()
            cov                                                 = globres.GetCovarianceMatrix()
            sidexbinwidth                                       = histclone.GetBinWidth(1)          # TODO: make more robust and general
           
            # calculate number of instances instead of integral
            
            # Subtract background (only if total fit worked but background fit failed)
            if backfit2.GetMinimum(fitrange[0],fitrange[1]) < -5: 
                npeak                                               = globfit.Integral(fitrange[0],fitrange[1]) / sidexbinwidth
                npeak_error                                         = globfit.IntegralError(fitrange[0],fitrange[1], params, cov.GetMatrixArray()) / sidexbinwidth
            else:
                npeak                                               = sigfit.Integral(fitrange[0],fitrange[1]) / sidexbinwidth
                npeak_error                                         = globfit.IntegralError(fitrange[0],fitrange[1], params, cov.GetMatrixArray()) / sidexbinwidth

            # Calculate peak width and error for confidence study
            confidence                                          = float(abs(globfit2.GetParameter(2)))
            conf_error                                          = float(abs(globfit2.GetParError(2)))
        else:                                                                                       # Cut and Count method
            histclone2                                          = hist.Clone()
            if backfit.GetMinimum(fitrange[0],fitrange[1]) > 0: histclone2.Add(backfit,-1)          # Subtract background (only if background fit is physical)
            npeak, npeak_error                                  = 0., 0.                            # Determine event count/statistical error
            
            # Count number of instances in each bin
            for i in range(0, nbins):
                npeak                                           += histclone2.GetBinContent(i)
                npeak_error                                     += np.power(histclone2.GetBinError(i),2)

            npeak_error                                         = np.sqrt(npeak_error)

            # Calculate peak width and error for confidence study
            confidence                                          = float(abs(histclone2.GetRMS()))
            conf_error                                          = float(abs(histclone2.GetRMSError()))
            #confidence                                          = float(abs(globfit2.GetParameter(2)))
            #conf_error                                          = float(abs(globfit2.GetParError(2)))
        
        print(f'Confidence = ', confidence, ' +- ', conf_error )

        return (npeak, npeak_error, confidence, conf_error)

    # METHOD 2: do global fit and determine integral of peak with error
    elif mode=='gfit':
        intpeak                                                 = np.sqrt(2*np.pi)*(paramdict['A_{1}']*paramdict['#sigma_{1}']
                                                                            + paramdict['A_{2}']*paramdict['#sigma_{2}'])
        sidexbinwidth                                           = histclone.GetBinWidth(1)           # TODO: make more robust and general
        npeak                                                   = intpeak / sidexbinwidth            # calculate number of instances instead of integral
        npeak_error                                             = globfit.IntegralError(fitrange[0],fitrange[1]) / sidexbinwidth
        
        return (npeak, npeak_error, confidence, conf_error)

    else:
        raise Exception('ERROR: peak counting mode not regognized!')

# -------------------------------------------------------------------------------------
# Wrap around fit function:
#   wrapper around fitting function using more modern paradigm for input parameters
#   Input:  - values: np array of sideband values
#           - weights: np array of weights corresponding to values
#           - variable: dict with all information about the sideband variable
#           - mode: passed down to called function
# -------------------------------------------------------------------------------------
def count_peak_unbinned(values, weights, variable, mode='subtract',
                        label=None, lumi=None, extrainfo=None,
                        histname='sideband', plotdir=None):
    # make a ROOT histogram with the values and weights
    counts  = np.histogram(values, variable['bins'], weights=weights)[0]
    errors  = np.sqrt(np.histogram(values, variable['bins'], weights=np.power(weights,2))[0])
    
    hist    = ROOT.TH1F(histname, histname, len(variable['bins'])-1, array('f', variable['bins']))
    hist.SetDirectory(0)
    for i, (count, error) in enumerate(zip(counts, errors)):
      hist.SetBinContent(   i+1, count)
      hist.SetBinError(     i+1, error)

    # make directory to store plots
    singleGaussdir = os.path.join(plotdir, 'singleGauss')
    if( plotdir is not None and not os.path.exists(plotdir) ):          os.makedirs(plotdir)
    if( plotdir is not None and not os.path.exists(singleGaussdir) ):   os.makedirs(singleGaussdir)

    # make a dictionary with extra fit info
    fitinfo                 = {}
    fitinfo['lumi']         = lumi
    fitinfo['helpdir']      = plotdir
    fitinfo['helpdir2']     = singleGaussdir

    # call underlying function
    return count_peak(hist, label, extrainfo, fitinfo, mode=mode)
