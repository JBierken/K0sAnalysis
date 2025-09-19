# Tools for pileup reweighting

See Twiki for more information: 'https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData'

# Centrally produced ROOT histograms for pileup reweighting

Although the users are in general invited to use the pileupCalc script themselves, to produce the ROOT histograms for pileup reweighting, in most cases they can simply pick the histograms in the following directories on lxplus or from the web:

### For Run3

It is important to realize that these are preliminary results and may change in the future, so make sure you use the latest ones.

2023 - BCD: /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/BCD/ \\
2023 - BC: /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/BC/ \\
2023 - D: /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/D/ \\

2022 - BCDEFG: /eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/BCDEFG/ \\
2022 - BCD: /eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/BCD/ \\
2022 - EFG: /eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/EFG/ \\

### For Run2

When using UltraLegacy samples:

2018: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/ \\
2017: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/ \\
2016: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/ \\ 
                                     
When using the previous samples and the previous luminosity calibration (note: the number of bins in this case is not fixed, because each MC campaign used a different setting):

2018: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/PrelLum13TeV/ \\
2017: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/PrelLum13TeV/ \\
2016: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/PrelLum15And1613TeV/ \\

In all those directories you will find (as indicated in the file name) histograms corresponding to the following values of the pp inelastic cross section: 69200 ub (recommended central value), 66000 ub (central value - 4.6%), 72400 ub (central value + 4.6%), 80000 ub (conventional value used for the public plots, agreed with ATLAS years ago).

Note also that, while all histograms in the UltraLegacy directories have the same number of bins (99, as indicated in the file name), this is not true in other cases. This reflects the number of pileup bins used for the generation of MC samples. 
