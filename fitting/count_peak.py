########################################################################################
# Count the number of instances in an invariant mass peak after background subtraction #
########################################################################################
# Note: this could potentially be updated to a new fitting library,
#       to remove all ROOT dependency.


import sys
import os
import math
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
# Helpers for the adaptive fit strategy and for the signal-only peak integral
# -------------------------------------------------------------------------------------

# significance level of the F-test that decides whether a more complex fit model
# (pol2 instead of pol1 background, double instead of single Gaussian) is needed
ADAPTIVE_FTEST_ALPHA = 0.05

# complexity rank of the signal model (None = background-only fit)
_SIGNAL_RANK = {None: 0, 'single': 0, 'double': 1}


def _fit_quality(fitfunc):
    ndf = int(fitfunc.GetNDF())
    chi2 = float(fitfunc.GetChisquare())
    if ndf <= 0 or not np.isfinite(chi2):
        return chi2, ndf, np.inf
    return chi2, ndf, chi2 / float(ndf)


def _fit_status(fitres):
    # decode the TFitResult status code; with the default (TMinuit) minimizer it is
    #   migrad + 10*minos + 100*hesse + 1000*improve
    # improve=4 ("no better minimum found") is normal and minos errors are not used,
    # so only the migrad and hesse digits and the covariance matrix status matter
    # (CovMatrixStatus: 3 = full accurate, 2 = forced positive-definite, 1 = approximate, 0 = none)
    try:
        status = int(fitres.Status())
        covstatus = int(fitres.CovMatrixStatus())
    except Exception:
        return None, None, None, None
    return status, status % 10, (status // 100) % 10, covstatus


def _make_candidate(model, degree, fitfunc, paramdict, fitobj, fitres):
    chi2, ndf, chi2ndf = _fit_quality(fitfunc)
    status, migrad, hesse, covstatus = _fit_status(fitres)
    converged = (migrad == 0) and np.isfinite(chi2ndf)
    usable = converged and (hesse == 0) and (covstatus == 3)
    return {
        'model': model,
        'background_degree': degree,
        'fitfunc': fitfunc,
        'paramdict': paramdict,
        'fitobj': fitobj,
        'fitres': fitres,
        'chi2': chi2,
        'ndf': ndf,
        'chi2_ndf': chi2ndf,
        'status': status,
        'covstatus': covstatus,
        'converged': converged,
        'usable': usable,
        'selection': None,
    }


def _candidate_name(candidate):
    if candidate['model'] is None:
        return 'pol{}'.format(candidate['background_degree'])
    return '{} + pol{}'.format(candidate['model'], candidate['background_degree'])


def _is_extension(simple, complex_):
    # True if complex_ contains simple as a special case (nested models)
    if simple is complex_:
        return False
    rank_s, rank_c = _SIGNAL_RANK[simple['model']], _SIGNAL_RANK[complex_['model']]
    deg_s, deg_c = simple['background_degree'], complex_['background_degree']
    return rank_c >= rank_s and deg_c >= deg_s and (rank_c, deg_c) != (rank_s, deg_s)


def _ftest_pvalue(simple, complex_):
    # F-test for nested chi2 fits: is the chi2 improvement of the more complex model
    # larger than expected from its extra free parameters?
    dk = simple['ndf'] - complex_['ndf']
    if dk <= 0 or complex_['ndf'] <= 0 or not complex_['chi2'] > 0:
        return 1.
    dchi2 = simple['chi2'] - complex_['chi2']
    if not dchi2 > 0:
        return 1.
    fstat = (dchi2 / dk) / (complex_['chi2'] / complex_['ndf'])
    return float(ROOT.Math.fdistribution_cdf_c(fstat, dk, complex_['ndf']))


def _select_model(candidates, alpha=ADAPTIVE_FTEST_ALPHA, tag='fit'):
    # select a fit model among nested candidates:
    #   - only fits that converged with an accurate covariance matrix are considered
    #   - start from the simplest one and move to a more complex model
    #     only if the F-test p-value is below alpha (smallest p-value first)
    # if no candidate is usable, fall back to the lowest finite chi2/ndf
    # (preferring converged fits) and print a warning.
    usable = [c for c in candidates if c['usable']]
    if len(usable) == 0:
        pool = ([c for c in candidates if c['converged']]
                or [c for c in candidates if np.isfinite(c['chi2_ndf'])]
                or list(candidates))
        best = min(pool, key=lambda c: c['chi2_ndf'] if np.isfinite(c['chi2_ndf']) else np.inf)
        best['selection'] = 'fallback'
        print('WARNING in adaptive {}: no candidate converged with an accurate covariance matrix'.format(tag)
              + ' (status/covstatus: {});'.format(', '.join('{}={}/{}'.format(
                  _candidate_name(c), c['status'], c['covstatus']) for c in candidates))
              + ' falling back to {} (lowest chi2/ndf)'.format(_candidate_name(best)))
        return best

    current = max(usable, key=lambda c: (c['ndf'], -_SIGNAL_RANK[c['model']]))
    steps = []
    while True:
        best_next, best_p = None, alpha
        for c in usable:
            if not _is_extension(current, c):
                continue
            p = _ftest_pvalue(current, c)
            if p < best_p:
                best_next, best_p = c, p
        if best_next is None:
            break
        steps.append('{} -> {} (p = {:.3g})'.format(_candidate_name(current), _candidate_name(best_next), best_p))
        current = best_next
    current['selection'] = 'ftest'
    skipped = [_candidate_name(c) for c in candidates if not c['usable']]
    print('Adaptive {} selection: {}{}{}'.format(
        tag, _candidate_name(current),
        '' if len(steps) == 0 else ' via ' + ', '.join(steps),
        '' if len(skipped) == 0 else ' [not usable: {}]'.format(', '.join(skipped))))
    return current


def _gauss_integral(mu, amplitude, sigma, xlow, xhigh):
    sigma = abs(sigma)
    if sigma == 0:
        return 0.
    s2 = sigma * np.sqrt(2.)
    return amplitude * sigma * np.sqrt(np.pi / 2.) * (math.erf((xhigh - mu) / s2) - math.erf((xlow - mu) / s2))


def _signal_integral(params, model, xlow, xhigh):
    # params follow the ordering of tools/fittools.py:
    #   single: [mu, A, sigma, a0, ...], double: [mu, A1, sigma1, A2, sigma2, a0, ...]
    if model == 'single':
        return _gauss_integral(params[0], params[1], params[2], xlow, xhigh)
    if model == 'double':
        return (_gauss_integral(params[0], params[1], params[2], xlow, xhigh)
                + _gauss_integral(params[0], params[3], params[4], xlow, xhigh))
    raise ValueError('Unknown signal model {}'.format(model))


def _signal_integral_and_error(fitres, model, fitrange):
    # integral of the signal component only, with its uncertainty propagated
    # from the full fit covariance matrix (so correlations between signal and
    # background parameters are taken into account): err^2 = g^T C g,
    # with g the gradient of the signal integral w.r.t. the signal parameters
    # (the integral does not depend on the background parameters).
    nsig = 3 if model == 'single' else 5
    params = np.array([float(fitres.Parameter(i)) for i in range(nsig)])
    cov = np.array([[float(fitres.CovMatrix(i, j)) for j in range(nsig)] for i in range(nsig)])
    integral = _signal_integral(params, model, fitrange[0], fitrange[1])
    grad = np.zeros(nsig)
    for i in range(nsig):
        step = 1e-6 * abs(params[i]) if params[i] != 0 else 1e-9
        up, down = params.copy(), params.copy()
        up[i] += step
        down[i] -= step
        grad[i] = (_signal_integral(up, model, fitrange[0], fitrange[1])
                   - _signal_integral(down, model, fitrange[0], fitrange[1])) / (2. * step)
    var = float(grad.dot(cov).dot(grad))
    error = np.sqrt(var) if (np.isfinite(var) and var >= 0) else np.nan
    return integral, error


def _peak_count_from_fit(globfit, globres, backfit, model, fitrange, binwidth):
    # number of instances in the peak (and its uncertainty) from a signal+background fit
    params = globres.GetParams()
    cov = globres.GetCovarianceMatrix()
    total_error = globfit.IntegralError(fitrange[0], fitrange[1], params, cov.GetMatrixArray()) / binwidth
    # background component unphysical: count the full fit integral (with its own error)
    if backfit.GetMinimum(fitrange[0], fitrange[1]) < -5:
        npeak = globfit.Integral(fitrange[0], fitrange[1]) / binwidth
        return npeak, total_error
    integral, error = _signal_integral_and_error(globres, model, fitrange)
    if not np.isfinite(error):
        print('WARNING in count_peak: signal integral error not computable'
              + ' (invalid covariance matrix?), using the error on the full fit integral instead')
        error = total_error * binwidth
    return integral / binwidth, error / binwidth


def _make_signal_fit(hist, fitrange, fitcenter, model, background_params, background_degree, funcname,
                     xcenter=0., xscale=1.):
    peak_height = max(hist.GetMaximum() / 2., 1e-6)
    peak_width = max(hist.GetRMS(), hist.GetBinWidth(1), 1e-6)

    if model == 'single':
        guess = [fitcenter, peak_height, peak_width]
        if background_degree == 1:
            guess += [background_params['a0'], background_params['a1']]
        else:
            guess += [background_params['a0'], background_params['a1'], background_params['a2']]
        return ft.poly_plus_gauss_fit(hist, fitrange, guess, funcname=funcname, xcenter=xcenter, xscale=xscale)

    if model == 'double':
        guess = [fitcenter, peak_height, peak_width / 4., peak_height, peak_width]
        if background_degree == 1:
            guess += [background_params['a0'], background_params['a1']]
        else:
            guess += [background_params['a0'], background_params['a1'], background_params['a2']]
        return ft.poly_plus_doublegauss_fit(hist, fitrange, guess, funcname=funcname, xcenter=xcenter, xscale=xscale)

    raise ValueError('Unknown signal model {}'.format(model))


def _background_params_for_degree(paramdict, degree):
    background_params = {'a0': paramdict['a0'], 'a1': paramdict['a1']}
    if degree >= 2:
        background_params['a2'] = paramdict.get('a2', 0.)
    return background_params


def _make_background_component(paramdict, degree, fitrange=None, xcenter=0., xscale=1.):
    # polynomial in t = (x-xcenter)/xscale, matching the parametrization of the fit
    t = '((x-({:.10g}))/({:.10g}))'.format(xcenter, xscale)
    formula = '[0]+[1]*{0}'.format(t) if degree == 1 else '[0]+[1]*{0}+[2]*{0}*{0}'.format(t)
    if fitrange is None:
        backfit = ROOT.TF1('fitfunc', formula)
    else:
        backfit = ROOT.TF1('fitfunc', formula, fitrange[0], fitrange[1])
    if degree == 1:
        backfit.SetParameters(paramdict['a0'], paramdict['a1'])
    else:
        backfit.SetParameters(paramdict['a0'], paramdict['a1'], paramdict['a2'])
    return backfit

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
    # Threshold for choosing between single and double Gaussian
    singleGauss = 20000

    # make background-only fit on the sidebands (peak region blanked in histclone)
    # If at most 1000 effective entries --> use linear background
    # Else                              --> use quadratic background
    if(hist.GetEffectiveEntries() <= 1000): guess = [0., 0.]
    else:                                   guess = [0.,0.,0.]
    backfit, paramdict, backfitobj = ft.poly_fit(histclone, fitrange, guess, "MSE0")
    
    if gargs['helpdir'] is not None:
        lumitext    = '' if gargs['lumi'] is None else '{0:.3g} '.format(float(gargs['lumi'])/1000.) + 'fb^{-1} (13.6 TeV)'
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
    # subtract mode still needs the signal fit to estimate the peak integral
    if(mode=='gfit' or mode=='hybrid' or mode=='subtract'):
        # If more than 10000 events     --> use Double Gauss
        # Else                          --> use Single Gauss
        if(hist.GetEffectiveEntries() <= singleGauss):
            guess = [
                    fitcenter,                                                                                      # peak position
                    hist.GetMaximum()/2,                                                                            # peak 1 height
                    hist.GetRMS()                                                                                   # peak 1 width
            ]
            if(hist.GetEffectiveEntries() <= 1000): guess += [paramdict['a0'],paramdict['a1']]                      # background estimate
            else:                                   guess += [paramdict['a0'],paramdict['a1'], paramdict['a2']]     # background estimate
            #guess += [paramdict['a0'],paramdict['a1']]                     # background estimate
            
            globfit, paramdict, globfitobji, globres = ft.poly_plus_gauss_fit(hist, fitrange, guess)

        else:
            guess = [
                    fitcenter,                                              # peak position
                    hist.GetMaximum()/2,                                    # peak 1 height
                    hist.GetRMS()/4,                                        # peak 1 width
                    hist.GetMaximum()/2,                                    # peak 2 height
                    hist.GetRMS()                                           # peak 2 width
            ]
            if(hist.GetEffectiveEntries() <= 1000): guess += [paramdict['a0'],paramdict['a1']]                      # background estimate
            else:                                   guess += [paramdict['a0'],paramdict['a1'], paramdict['a2']]     # background estimate
           
            globfit, paramdict, globfitobj, globres = ft.poly_plus_doublegauss_fit(hist, fitrange, guess)
            print(paramdict) 
        
        # separate background component of fit
        if(hist.GetEffectiveEntries() <= 1000):  
            backfit2 = ROOT.TF1("fitfunc","pol1(0)", fitrange[0], fitrange[1])
            backfit2.SetParameters(paramdict['a0'], paramdict['a1'])
        else:
            backfit2 = ROOT.TF1("fitfunc","pol2(0)", fitrange[0], fitrange[1])
            backfit2.SetParameters(paramdict['a0'], paramdict['a1'], paramdict['a2'])
        
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
    if(hist.GetEffectiveEntries() <= 1000): guess += [paramdict['a0'],paramdict['a1']]                      # background estimate
    else:                                   guess += [paramdict['a0'],paramdict['a1'], paramdict['a2']]     # background estimate
           
    globfit2, paramdict2, globfitobj2, globres2 = ft.poly_plus_gauss_fit(hist, fitrange, guess)
    
    # separate background component of fit
    if(hist.GetEffectiveEntries() <= 1000):  
        backfit3 = ROOT.TF1("fitfunc","pol1(0)", fitrange[0], fitrange[1])
        backfit3.SetParameters(paramdict['a0'], paramdict['a1'])
    else:
        backfit3 = ROOT.TF1("fitfunc","pol2(0)", fitrange[0], fitrange[1])
        backfit3.SetParameters(paramdict['a0'], paramdict['a1'], paramdict['a2'])
    
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
        if hist.GetEffectiveEntries()                           > 50:                               # Calculate integral
            if hist.GetEffectiveEntries()                       <= singleGauss: 
                sigfit                                          = ROOT.TF1("signal","gaus(0)")
                sigfit.SetParameters(paramdict[r'A'], paramdict[r'#mu'], paramdict[r'#sigma'])
                print(f'sigma = ', abs(globfit.GetParameter(2)), ' +- ', abs(globfit2.GetParError(2)))
            else:                                
                sigfit                                          = ROOT.TF1("signal","gaus(0)+gaus(3)")
                sigfit.SetParameters(paramdict[r'A_{1}'], paramdict[r'#mu'], paramdict[r'#sigma_{1}'],
                                                                paramdict[r'A_{2}'], paramdict[r'#mu'], paramdict[r'#sigma_{2}'])
                print(f'sigma_1 = ', abs(globfit.GetParameter(2)), ' +- ', abs(globfit2.GetParError(2)))
                print(f'sigma_2 = ', abs(globfit.GetParameter(4)), ' +- ', abs(globfit2.GetParError(4)))
            
            sidexbinwidth                                       = histclone.GetBinWidth(1)          # TODO: make more robust and general

            # calculate number of instances (signal component only, with its own error;
            # full fit integral if the background component is unphysical)
            signal_model                                        = 'single' if hist.GetEffectiveEntries() <= singleGauss else 'double'
            npeak, npeak_error                                  = _peak_count_from_fit(globfit, globres, backfit2, signal_model,
                                                                                       fitrange, sidexbinwidth)

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
        
        print(f'Integral = ',   npeak,      ' +- ', npeak_error)
        print(f'Confidence = ', confidence, ' +- ', conf_error )

        return (npeak, npeak_error, confidence, conf_error)
        #return (npeak, npeak_error)

    # METHOD 2: do global fit and determine integral of peak with error
    elif mode=='gfit':
        if hist.GetEffectiveEntries() <= singleGauss:
            intpeak                                             = np.sqrt(2*np.pi) * paramdict['A'] * paramdict['#sigma']
        else:
            intpeak                                             = np.sqrt(2*np.pi)*(paramdict['A_{1}']*paramdict['#sigma_{1}']
                                                                                + paramdict['A_{2}']*paramdict['#sigma_{2}'])
        sidexbinwidth                                           = histclone.GetBinWidth(1)           # TODO: make more robust and general
        npeak                                                   = intpeak / sidexbinwidth            # calculate number of instances instead of integral
        npeak_error                                             = globfit.IntegralError(fitrange[0],fitrange[1]) / sidexbinwidth
        
        return (npeak, npeak_error, confidence, conf_error)
        #return (npeak, npeak_error)

    else:
        raise Exception('ERROR: peak counting mode not regognized!')


# -------------------------------------------------------------------------------------
# Adaptive (threshold-free) version of count_peak:
#   - background-only fit: pol1 and pol2 on the sidebands
#   - signal+background fit: {single, double} Gaussian x {pol1, pol2}
#   - in both cases, only fits that converged with an accurate covariance matrix are
#     considered, and a more complex model is only chosen if it improves the fit
#     significantly (F-test for nested models, see _select_model and ADAPTIVE_FTEST_ALPHA)
#   - the peak width ("confidence") is taken from the selected single Gaussian + polN fit
#   - the peak count and its error refer to the signal component only
# -------------------------------------------------------------------------------------
def count_peak_adaptive(hist, label, extrainfo, gargs, mode='subtract', return_metadata=False):
    # initializations
    nbins           = hist.GetNbinsX()
    xlow            = hist.GetBinLowEdge(1)
    xhigh           = hist.GetBinLowEdge(nbins)+hist.GetBinWidth(nbins)
    fitrange        = (xlow, xhigh)
    fitcenter       = (xlow + xhigh)/2.
    fithalfwidth    = (xhigh - xlow)/4.
    histclone       = hist.Clone()
    binlow          = histclone.FindBin(fitcenter-fithalfwidth)
    binhigh         = histclone.FindBin(fitcenter+fithalfwidth)
    lumitext        = '' if gargs['lumi'] is None else '{0:.3g} '.format(float(gargs['lumi'])/1000.) + 'fb^{-1} (13.6 TeV)'
    # background polynomials are parametrized in t = (x-fitcenter)/bkgxscale, i.e. t in [-1,1],
    # which keeps the coefficients nearly uncorrelated (in raw x, the pol2 covariance matrix
    # is so ill-conditioned that HESSE fails and the fit would be rejected as not usable)
    bkgxscale       = (xhigh - xlow)/2.

    # take into account only sidebands, not peak
    for i in range(binlow, binhigh):
        histclone.SetBinContent(i,0)
        histclone.SetBinError(i,0)

    # background-only fits on the sidebands: used as seed for the signal fits,
    # for the background plot and for the cut-and-count fallback
    backfit_candidates = []
    for degree, guess, funcname in (
        (1, [0., 0.], 'backfit_pol1'),
        (2, [0., 0., 0.], 'backfit_pol2'),
    ):
        fitfunc, fitparamdict, fitobj, fitres = ft.poly_fit(histclone, fitrange, guess, "MSE0",
                                                            funcname=funcname, return_fitresult=True,
                                                            xcenter=fitcenter, xscale=bkgxscale)
        backfit_candidates.append(_make_candidate(None, degree, fitfunc, fitparamdict, fitobj, fitres))

    best_backfit = _select_model(backfit_candidates, tag='background fit')
    best_backfit_params = best_backfit['paramdict']
    best_backfit_degree = best_backfit['background_degree']
    backfit = best_backfit['fitfunc']
    paramdict = best_backfit_params

    if gargs['helpdir'] is not None:
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

    # signal+background fits
    # (the single Gaussian fits are always needed for the peak width;
    #  the double Gaussian fits only if the fit result is used for the peak count)
    need_signal_fit = (mode == 'gfit' or mode == 'hybrid' or hist.GetEffectiveEntries() > 50)
    models = ('single', 'double') if need_signal_fit else ('single',)
    signal_candidates = []
    for background_degree in (1, 2):
        background_seed = _background_params_for_degree(best_backfit_params, background_degree)
        for model in models:
            funcname = 'adaptive_{}_pol{}'.format(model, background_degree)
            fitfunc, fitparamdict, fitobj, fitres = _make_signal_fit(
                hist, fitrange, fitcenter, model, background_seed, background_degree, funcname,
                xcenter=fitcenter, xscale=bkgxscale
            )
            signal_candidates.append(_make_candidate(model, background_degree, fitfunc, fitparamdict, fitobj, fitres))

    # peak width ("confidence"): selected single Gaussian + polN fit
    best_confidence = _select_model([c for c in signal_candidates if c['model'] == 'single'], tag='confidence fit')
    globfit2 = best_confidence['fitfunc']
    paramdict2 = best_confidence['paramdict']
    confidence_background_degree = best_confidence['background_degree']
    backfit3 = _make_background_component(paramdict2, confidence_background_degree, fitrange,
                                          xcenter=fitcenter, xscale=bkgxscale)

    best_signal = None
    if need_signal_fit:
        best_signal = _select_model(signal_candidates, tag='signal fit')

        globfit = best_signal['fitfunc']
        paramdict = best_signal['paramdict']
        globres = best_signal['fitres']
        signal_model = best_signal['model']
        signal_background_degree = best_signal['background_degree']

        print('Adaptive signal fit: {} + pol{} with chi2/ndf = {:.4g}/{:d} = {:.4g}'.format(
            signal_model, signal_background_degree, best_signal['chi2'], best_signal['ndf'], best_signal['chi2_ndf']))

        backfit2 = _make_background_component(paramdict, signal_background_degree, fitrange,
                                              xcenter=fitcenter, xscale=bkgxscale)

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

    print('Adaptive confidence fit: single + pol{} with chi2/ndf = {:.4g}/{:d} = {:.4g}'.format(
        confidence_background_degree, best_confidence['chi2'], best_confidence['ndf'], best_confidence['chi2_ndf']))

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

    def _metadata(npeak, npeak_error, confidence, conf_error):
        return {
            'peak_count': npeak,
            'peak_error': npeak_error,
            'confidence': confidence,
            'conf_error': conf_error,
            'background_model': 'pol{}'.format(best_backfit_degree),
            'signal_model': signal_model,
            'signal_background_model': 'pol{}'.format(signal_background_degree),
            'signal_selection': best_signal['selection'],
            'signal_status': best_signal['status'],
            'signal_covstatus': best_signal['covstatus'],
            'background_chi2': best_backfit['chi2'],
            'background_ndf': best_backfit['ndf'],
            'background_chi2_ndf': best_backfit['chi2_ndf'],
            'signal_chi2': best_signal['chi2'],
            'signal_ndf': best_signal['ndf'],
            'signal_chi2_ndf': best_signal['chi2_ndf'],
            'confidence_chi2': best_confidence['chi2'],
            'confidence_ndf': best_confidence['ndf'],
            'confidence_chi2_ndf': best_confidence['chi2_ndf'],
        }

    sidexbinwidth = histclone.GetBinWidth(1)    # TODO: make more robust and general

    if mode == 'subtract' or mode == 'hybrid':
        if hist.GetEffectiveEntries() > 50 and best_signal is not None:
            npeak, npeak_error = _peak_count_from_fit(globfit, globres, backfit2, signal_model,
                                                      fitrange, sidexbinwidth)
            confidence = float(abs(globfit2.GetParameter(2)))
            conf_error = float(abs(globfit2.GetParError(2)))
        else:
            histclone2 = hist.Clone()
            if backfit.GetMinimum(fitrange[0],fitrange[1]) > 0:
                histclone2.Add(backfit,-1)
            npeak, npeak_error = 0., 0.

            for i in range(0, nbins):
                npeak += histclone2.GetBinContent(i)
                npeak_error += np.power(histclone2.GetBinError(i),2)

            npeak_error = np.sqrt(npeak_error)
            confidence = float(abs(histclone2.GetRMS()))
            conf_error = float(abs(histclone2.GetRMSError()))

        print(f'Integral = ',   npeak,      ' +- ', npeak_error)
        print(f'Confidence = ', confidence, ' +- ', conf_error )

        if return_metadata and best_signal is not None:
            return _metadata(npeak, npeak_error, confidence, conf_error)
        return (npeak, npeak_error, confidence, conf_error)

    elif mode == 'gfit':
        if signal_model == 'single':
            intpeak = np.sqrt(2*np.pi) * paramdict['A'] * paramdict[r'#sigma']
        else:
            intpeak = np.sqrt(2*np.pi) * (paramdict['A_{1}']*paramdict[r'#sigma_{1}'] + paramdict['A_{2}']*paramdict[r'#sigma_{2}'])

        npeak = intpeak / sidexbinwidth
        npeak_error = globfit.IntegralError(fitrange[0],fitrange[1]) / sidexbinwidth

        confidence = float(abs(globfit2.GetParameter(2)))
        conf_error = float(abs(globfit2.GetParError(2)))

        if return_metadata and best_signal is not None:
            return _metadata(npeak, npeak_error, confidence, conf_error)
        return (npeak, npeak_error, confidence, conf_error)

    else:
        raise Exception('ERROR: peak counting mode not regognized!')


def count_peak_compare(hist, label, extrainfo, gargs, mode='subtract'):
    compare_gargs = dict(gargs)
    compare_gargs['helpdir'] = None
    compare_gargs['helpdir2'] = None
    
    # Get adaptive result with metadata
    adaptive_result = count_peak_adaptive(hist, label, extrainfo, compare_gargs, mode=mode, return_metadata=True)
    
    # Get legacy result
    legacy_tuple = count_peak(hist, label, extrainfo, compare_gargs, mode=mode)
    
    # Infer legacy model choice based on effective entries (same logic as count_peak)
    effective_entries = hist.GetEffectiveEntries()
    legacy_background_model = 'pol1' if effective_entries <= 1000 else 'pol2'
    legacy_signal_model = 'single' if effective_entries <= 20000 else 'double'
    
    legacy_result = {
        'peak_count': legacy_tuple[0],
        'peak_error': legacy_tuple[1],
        'confidence': legacy_tuple[2],
        'conf_error': legacy_tuple[3],
        'background_model': legacy_background_model,
        'signal_model': legacy_signal_model,
        'background_chi2': None,  # Not tracked in legacy
        'background_ndf': None,
        'background_chi2_ndf': None,
        'signal_chi2': None,
        'signal_ndf': None,
        'signal_chi2_ndf': None,
        'confidence_chi2': None,
        'confidence_ndf': None,
        'confidence_chi2_ndf': None,
    }
    
    return {
        'legacy': legacy_result,
        'adaptive': adaptive_result,
        'mode': mode,
        'label': label,
        'extrainfo': extrainfo,
        'histogram_name': hist.GetName(),
        'effective_entries': effective_entries,
    }

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
                        histname='sideband', plotdir=None,
                        fit_strategy='legacy'):
    # make a ROOT histogram with the values and weights
    counts  = np.histogram(values, variable['bins'], weights=weights)[0]
    errors  = np.sqrt(np.histogram(values, variable['bins'], weights=np.power(weights,2))[0])
    
    hist    = ROOT.TH1F(histname, histname, len(variable['bins'])-1, array('f', variable['bins']))
    hist.SetDirectory(0)
    for i, (count, error) in enumerate(zip(counts, errors)):
      hist.SetBinContent(   i+1, count)
      hist.SetBinError(     i+1, error)

    # diagnostic: compare numpy counts and ROOT histogram
    try:
        heff = hist.GetEffectiveEntries()
    except Exception:
        heff = None
    print('COUNT_PEAK_UNBINNED_DEBUG histname={} np_sum_counts={} hist.Integral={} hist.GetEffectiveEntries={}'.format(
        histname, np.sum(counts), hist.Integral(), heff))

    # make directory to store plots
    singleGaussdir = None
    if plotdir is not None:
        singleGaussdir = os.path.join(plotdir, 'singleGauss')
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        if not os.path.exists(singleGaussdir):
            os.makedirs(singleGaussdir)

    # make a dictionary with extra fit info
    fitinfo                 = {}
    fitinfo['lumi']         = lumi
    fitinfo['helpdir']      = plotdir
    fitinfo['helpdir2']     = singleGaussdir
    fitinfo['fitstrategy']  = fit_strategy

    # decorate extrainfo so plots display which fit strategy was used
    if extrainfo is None:
        extrainfo = ''
    if fit_strategy is not None:
        extrainfo = 'Fit strategy: {}<<'.format(fit_strategy) + extrainfo
    
    # call underlying function
    if fit_strategy == 'adaptive':
        return count_peak_adaptive(hist, label, extrainfo, fitinfo, mode=mode)
    if fit_strategy == 'compare':
        return count_peak_compare(hist, label, extrainfo, fitinfo, mode=mode)

    return count_peak(hist, label, extrainfo, fitinfo, mode=mode)
