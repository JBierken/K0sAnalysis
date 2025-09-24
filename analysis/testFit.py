import ROOT
import numpy as np
from functools import partial

# Open input ROOT file
inName          = 'test_2.root'
inFile          = ROOT.TFile.Open(inName, "read")
#hist            = inFile.Get("2017 data_bin2")
keylist         = inFile.GetListOfKeys()
for key in keylist:
    print(key.GetName())
    hist        = inFile.Get(key.GetName())
    print('--')
# get histogram parameters
nbins           = hist.GetNbinsX()
xlow            = hist.GetBinLowEdge(1)
xhigh           = hist.GetBinLowEdge(nbins)+hist.GetBinWidth(nbins)
fitrange        = (xlow, xhigh)

def poly(x, par, degree=0):
    # par[k] = coefficient with x**k
    # note: in new versions of ROOT or CMSSW or python,
    #       one cannot rely on len(par) anymore to get the degree
    #       of the polynomial; instead it is passed here explicitly.
    res = 0.
    for k in range(degree+1):
        res += par[k]*np.power(x[0], k)
    return res

# set up fit parameters
guess           = [0., 0.]
degree          = len(guess)-1
fitobj          = partial(poly, degree=degree)
fitfunc         = ROOT.TF1("fitfunc", fitobj, fitrange[0], fitrange[1], len(guess))
for i,val in enumerate(guess):
    fitfunc.SetParameter(i,val)

# do fit
fitresult       = hist.Fit(fitfunc)
#fitresult       = hist.Fit("gaus", "WLQ0")
#fitresult       = hist.Fit("gaus", "W0")

# write hist to file (first check)
outHistFile     = ROOT.TFile.Open("test_4.root", "recreate")
outHistFile.cd()
hist.Write()
fitfunc.Write()
outHistFile.Close()
