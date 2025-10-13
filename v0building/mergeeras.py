#########################
# Merge eras into years #
#########################
# WARNING: lots of hard-coding, might need to adapt if running of different eras 
# or with different naming conventions.

import sys
import os
import argparse


if __name__=='__main__':
    
    # command line arguments
    parser = argparse.ArgumentParser( description = 'Merge V0 files' )
    parser.add_argument('-i', '--filedir',      required=True,      type=os.path.abspath)
    parser.add_argument('-o', '--outputfile',   default='selected.root')
    parser.add_argument('-v', '--version',      default='run2ul')
    parser.add_argument('-f', '--force',        default=False,      action='store_true')
    args = parser.parse_args()

    # Define years to be merged together. --> Hardcoded for now, should be improved
    if args.version == 'run2preul':
        yeardict = ({
            'DoubleMuon_Run2016': [
                'DoubleMuon_crab_Run2016B-17Jul2018_ver2-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016C-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016D-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016E-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016F-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016G-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016H-17Jul2018-v1_data_Run2016_DoubleMuon'
            ],
            'DoubleMuon_Run2016BtoF': [
                'DoubleMuon_crab_Run2016B-17Jul2018_ver2-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016C-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016D-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016E-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016F-17Jul2018-v1_data_Run2016_DoubleMuon',
            ],
            'DoubleMuon_Run2016GtoH': [
                'DoubleMuon_crab_Run2016G-17Jul2018-v1_data_Run2016_DoubleMuon',
                'DoubleMuon_crab_Run2016H-17Jul2018-v1_data_Run2016_DoubleMuon'
            ],
            'DoubleMuon_Run2017': [
                'DoubleMuon_crab_Run2017B-31Mar2018-v1_data_Run2017_DoubleMuon',
                'DoubleMuon_crab_Run2017C-31Mar2018-v1_data_Run2017_DoubleMuon',
                'DoubleMuon_crab_Run2017D-31Mar2018-v1_data_Run2017_DoubleMuon',
                'DoubleMuon_crab_Run2017E-31Mar2018-v1_data_Run2017_DoubleMuon',
                'DoubleMuon_crab_Run2017F-31Mar2018-v1_data_Run2017_DoubleMuon'
            ],
            'DoubleMuon_Run2018': [
                'DoubleMuon_crab_Run2018A-17Sep2018-v2_data_Run2018_DoubleMuon',
                'DoubleMuon_crab_Run2018B-17Sep2018-v1_data_Run2018_DoubleMuon',
                'DoubleMuon_crab_Run2018C-17Sep2018-v1_data_Run2018_DoubleMuon',
                'DoubleMuon_crab_Run2018D-PromptReco-v2_data_Run2018_DoubleMuon'
            ]
        })
    elif args.version == 'run2ul':
        yeardict = ({
            'DoubleMuon_Run2016PreVFP': [
                'DoubleMuon_crab_Run2016B-ver2_HIPM_UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016C-HIPM_UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016D-HIPM_UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016E-HIPM_UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016F-HIPM_UL2016_MiniAODv2-v1'
            ],
            'DoubleMuon_Run2016PostVFP': [
                'DoubleMuon_crab_Run2016F-UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016G-UL2016_MiniAODv2-v1',
                'DoubleMuon_crab_Run2016H-UL2016_MiniAODv2-v2'
            ],
            'DoubleMuon_Run2017': [
                'DoubleMuon_crab_Run2017B-UL2017_MiniAODv2-v1',
                'DoubleMuon_crab_Run2017C-UL2017_MiniAODv2-v1',
                'DoubleMuon_crab_Run2017D-UL2017_MiniAODv2-v1',
                'DoubleMuon_crab_Run2017E-UL2017_MiniAODv2-v2',
                'DoubleMuon_crab_Run2017F-UL2017_MiniAODv2-v1'
            ],
            'DoubleMuon_Run2018': [
                'DoubleMuon_crab_Run2018A-UL2018_MiniAODv2-v1_data_doublemuon',
                'DoubleMuon_crab_Run2018B-UL2018_MiniAODv2-v1_data_doublemuon',
                'DoubleMuon_crab_Run2018C-UL2018_MiniAODv2-v1_data_doublemuon',
                'DoubleMuon_crab_Run2018D-UL2018_MiniAODv2-v1_data_doublemuon'
            ]
        })
    elif args.version=='run3':
        yeardict = ({
            'Run2022preEE': [
                'DoubleMuon_Run2022B-22Sep2023-v1',
                'DoubleMuon_Run2022C-22Sep2023-v1',
                'Muon_Run2022C-22Sep2023-v1',
                'Muon_Run2022D-22Sep2023-v1'
            ],
            'Run2022postEE': [
                '_Run2022E-22Sep2023-v1',
                '_Run2022F-22Sep2023-v2',
                '_Run2022G-22Sep2023-v1'
            ],
            'Run2023preBPix': [
                'Muon0_Run2023B-22Sep2023-v1',
                'Muon0_Run2023C-22Sep2023_v1-v1',
                'Muon0_Run2023C-22Sep2023_v2-v1',
                'Muon0_Run2023C-22Sep2023_v3-v1',
                'Muon0_Run2023C-22Sep2023_v4-v1',
                'Muon1_Run2023B-22Sep2023-v1',
                'Muon1_Run2023C-22Sep2023_v1-v1',
                'Muon1_Run2023C-22Sep2023_v2-v1',
                'Muon1_Run2023C-22Sep2023_v3-v1',
                'Muon1_Run2023C-22Sep2023_v4-v2'
            ],
            'Run2023postBPix': [
                'Muon0_Run2023D-22Sep2023_v2-v1',
                'Muon1_Run2023D-22Sep2023_v1-v1',
                'Muon1_Run2023D-22Sep2023_v2-v1'
            ]
            #Other years...
        })

    for yeardir, eradirs in yeardict.items():
        yeardir = os.path.join(args.filedir, yeardir)
        # make output directory
        if not os.path.exists(yeardir): os.makedirs(yeardir)
        # make hadd command
        cmd = 'hadd'
        if args.force: cmd += ' -f'
        cmd += ' {}'.format(os.path.join(yeardir, args.outputfile))
        for eradir in eradirs:
            cmd += ' {}*/*.root'.format(os.path.join(args.filedir, eradir))
        # run hadd command
        os.system(cmd)
