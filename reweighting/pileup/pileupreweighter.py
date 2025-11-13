######################
# Pileup reweighting #
######################

import os
import sys
import numpy as np
import ROOT
#sys.path.append('/user/jbierken/CMSSW_12_4_12/src/K0sAnalysis/tools')
#import pileuptools as pu


def get_pileup_profile(campaign, year):
    if campaign=='run2preul':
        ### help function to get correct pileup profile fil (should be better)
        path                    = '/user/jbierken/CMSSW_14_0_15/src/K0sAnalysis/reweighting/pileup/data'
        path                    = os.path.join(path, campaign)
        
        if year.startswith('2016'):     year = '2016'
        if year not in ['2016','2017','2018']:
            msg                 = 'ERROR: year {} not recogized.'.format(year)
            raise Exception(msg)
        
        pufile                  = os.path.join(path,'dataPuHist_{}Inclusive_central.root'.format(year))
        histname                = 'pileup'
    
    elif campaign=='run2ul':
        
        ### help function to get correct pileup profile fil (should be better)
        path                    = '/user/jbierken/CMSSW_14_0_15/src/K0sAnalysis/reweighting/pileup/'
        
        if year.startswith('2016'):     year = '2016'
        if year not in ['2016','2017','2018']:
            msg                 = 'ERROR: year {} not recogized.'.format(year)
            raise Exception(msg)
        

        pufile                  = os.path.join(path, 'UL{}'.format(year), 'Pileup_ratio.root')
        histname                = 'pileup'

    elif campaign=='run3':
        
        ### help function to get correct pileup profile fil (should be better)
        path                    = '/user/jbierken/CMSSW_14_0_15/src/K0sAnalysis/reweighting/pileup/'
        
        #if year.startswith('2022'):     year = '2022'
        #if year.startswith('2023'):     year = '2023'
        if year not in ['2022', '2022preEE', '2022postEE', '2022combined', '2023', '2023combined', '2023preBPix', '2023postBPix', '2024']:
            msg                 = 'ERROR: year {} not recogized.'.format(year)
            raise Exception(msg)

        year                    = year.rstrip('prepostcombinedABCDEFGHI')
        pufile                  = os.path.join(path, '{}'.format(year), 'Pileup_ratio.root')
        histname                = 'pileup'

    else:
        msg                     = 'ERROR: campaign {} not recognized.'.format(campaign)
        raise Exception(msg)

    if not os.path.exists(pufile):
        msg                     = 'ERROR: pileup file {} does not seem to exist.'.format(pufile)
        raise Exception(msg)
    
    return (pufile, histname)


class PileupReweighter(object):

    def __init__(self, campaign, year, pufile=None, histname=None):
        
        ### initialize a pileup reweighter
        self.campaign           = campaign
        self.year               = year

        # get correct file
        if( pufile is None or histname is None ):   (self.pufile, self.histname) = get_pileup_profile(campaign, year)
        else:                                       (self.pufile, self.histname) = (pufile, histname)

        # open the file and get the histogram
        f                       = ROOT.TFile.Open(self.pufile)

        try:
            self.puhist         = f.Get(self.histname)
            _                   = self.puhist.GetBinContent(0)
        except:
            msg                 = 'ERROR: pileup profile in data could not be loaded.'
            raise Exception(msg)

        self.puhist.SetDirectory(ROOT.gROOT)
        f.Close()

        # initialize scale histogram
        # note: for pre-UL analyses, the scale histogram should be calculated
        #       for each sample separately, given the pileup distribution in data
        #       and the one for the specific sample;
        #       while for UL analyses, the pileup histogram loaded above
        #       gives directly the required ratio.

        self.scalehist                  = None
        if self.campaign=='run2ul' or self.campaign=='run3':   self.scalehist = self.puhist.Clone()

    def initsample(self, sample):
        ### initialize the reweighter for a given sample

        # skip in case of UL sample, no initialization needed
        if self.campaign=='run2ul':     return
        if self.campaign=='run2ul' or self.campaign=='run3':   return 

        # get true interaction profile from sample
        f                       = ROOT.TFile.Open(sample)
        try:
            inthist             = f.Get('nTrueInteractions')
            _                   = inthist.GetBinContent(0)
        except:
            msg                 = 'ERROR: interactions profile in simulation could not be loaded.'
            raise Exception(msg)

        inthist.SetDirectory(ROOT.gROOT)
        inthist.Scale(1./inthist.GetSumOfWeights())
        f.Close()

        # determine reweighting factors
        self.scalehist          = self.puhist.Clone()
        self.scalehist.Scale(1./self.puhist.GetSumOfWeights())
        self.scalehist.Divide(inthist)

    def getreweight(self, ntrueint):
        ### get the pileup reweighting factor for a given number of true interactions
        if self.scalehist is None:
            msg                 = 'ERROR: pileup reweighter not yet initialized with a sample'
            raise Exception(msg)

        isscalar                = False
        if( isinstance(ntrueint, int) or isinstance(ntrueint, float) ):
            isscalar            = True
            ntrueint            = np.array([ntrueint])

        res                     = np.zeros(len(ntrueint))
        for i in range(len(ntrueint)):
            res[i]              = self.scalehist.GetBinContent(self.scalehist.GetXaxis().FindBin(ntrueint[i]))
        if isscalar: res        = res[0]
        return res
