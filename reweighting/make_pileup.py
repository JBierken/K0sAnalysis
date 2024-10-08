#!/usr/bin/env python
from __future__ import print_function
import os 
import sys
import uproot
import uproot3
from uproot3_methods.classes import TH1
import numpy as np
import argparse

# Get BaseDir
baseDir     = '/user/jbierken/CMSSW_12_4_12/src/K0sAnalysis/reweighting'

def get_args():
    # get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--era',          required=True)
    parser.add_argument('-m', '--makeMC',       default=False,      action = 'store_true')
    
    return parser.parse_args()

# -------------------------------------------------------------------------
def make_pileup_mc(era,  base_dir=baseDir):
    # Get correct setup for each era
    if era == 'UL2016':
        from SimGeneral\
                .MixingModule\
                .mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi import mix
    elif era == 'UL2017':
        from SimGeneral\
                .MixingModule\
                .mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi import mix
    elif era == 'UL2018':
        from SimGeneral\
                .MixingModule\
                .mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi import mix
    elif era == '2016':
        from SimGeneral\
                .MixingModule\
                .mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import mix
    elif era == '2017':
        from SimGeneral\
                .MixingModule\
                .mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi import mix
    elif era == '2018':
        from SimGeneral\
                .MixingModule\
                .mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi import mix
    else:
        print('Unrecognized era', args.era)
        sys.exit(0)

    # Calculate pileup for MC
    values  = np.array([float(x) for x in mix.input.nbPileupEvents.probValue])
    edges   = np.arange(len(values)+1)
    hist    = TH1.from_numpy((values, edges))

    # Write pileup to ROOT file
    path = '{dir}/pileup/{era}/'.format(dir=baseDir, era=args.era)
    os.system('mkdir -p {path}'.format(path=path))
    with uproot3.recreate('{path}/mcPileup.root'.format(path=path)) as f:
         f['pileup'] = hist

# -------------------------------------------------------------------------
def get_pileup(era, base_dir=baseDir): # code below taken from spark tnp

    '''
    Get the pileup distribution scalefactors to apply to simulation
    for a given era.
    '''

    # get the pileup
    dataPileup = {
            'UL2016_APV':   base_dir + '/pileup/UL2016/dataPileup_nominal.root',
            'UL2016':       base_dir + '/pileup/UL2016/dataPileup_nominal.root',
            'UL2017':       base_dir + '/pileup/UL2017/dataPileup_nominal.root',
            'UL2018':       base_dir + '/pileup/UL2018/dataPileup_nominal.root'
            }
    mcPileup = {
            'UL2016_APV':   base_dir + '/pileup/UL2016/mcPileup.root',
            'UL2016':       base_dir + '/pileup/UL2016/mcPileup.root',
            'UL2017':       base_dir + '/pileup/UL2017/mcPileup.root',
            'UL2018':       base_dir + '/pileup/UL2018/mcPileup.root'
            }
    
    with uproot.open(dataPileup[era]) as f:
        data_edges  = f['pileup'].axis(0).edges()
        data_pileup = f['pileup'].values()
        data_pileup /= sum(data_pileup)
    
    with uproot.open(mcPileup[era]) as f:
        mc_edges    = f['pileup'].axis(0).edges()
        mc_pileup   = f['pileup'].values()
        mc_pileup   /= sum(mc_pileup)
    
    pileup_edges    = data_edges if len(data_edges) < len(mc_edges) else mc_edges
    pileup_ratio    = [d/m if m else 1.0 for d, m in zip(
                        data_pileup[:len(pileup_edges)-1], mc_pileup[:len(pileup_edges)-1])]
    
    hist    = TH1.from_numpy((np.array(pileup_ratio), np.array(pileup_edges)))

    # Write pileup to ROOT file
    path = '{dir}/pileup/{era}/'.format(dir=base_dir, era=args.era)
    os.system('mkdir -p {path}'.format(path=path))
    with uproot3.recreate('{path}/Pileup_ratio.root'.format(path=path)) as f:
     f['pileup'] = hist


    #return pileup_ratio, pileup_edges

# -------------------------------------------------------------------------
if __name__ == "__main__":

    # get arguments
    args            = get_args()

    sys.stderr.write('###starting###\n')
    
    # make MC pileup histogram
    if args.makeMC: make_pileup_mc(era=args.era)

    # Make pileup ratio histogram
    get_pileup(era=args.era)
    
    sys.stderr.write('###done###\n')
