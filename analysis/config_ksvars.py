# note: this is supposed to be the same configuration as the pre-UL version,
#       only with a different (newer) naming convention for the variables.

import numpy as np

config = ({
     'rpv': {'variablename':'_RPV',
            'xaxtitle': '#Delta_{2D} (cm)',
            'yaxtitle': 'Reconstructed vertices',
            'treename': 'laurelin',
            'extrainfos': ['K^{0}_{S} candidates'],
            'bkgmodes': {
              'bkgdefault': {'type':None, 'info':'Background not subtracted'},
              'bkgsideband': {'type':'sideband', 'info': 'Background subtracted',
                              'sidevariable': '_mass',
                              'sidebins': np.linspace(0.44, 0.56, num=31, endpoint=True)},
            },
            'bins': {
              'finebins':np.linspace(0,20,num=41,endpoint=True),
              'defaultbins':[0.,0.5,1.5,4.,10.,20.],
            },
            'normalization': {
              'normlumi':{'type':'lumi', 'info':'Normalized to luminosity'},
              'normeventyield':{'type':'eventyield', 'info':'Normalized to data events'},
              'normrange':{'type':'range', 'info': 'Normalized for #Delta_{2D} < 0.5 cm', 
                          'normvariable': '_RPV', 'normrange':[0.,0.5]},
            },
    },   
     'ks_pt': {'variablename':'_pt',
            'xaxtitle':'K^{0}_{S} transverse momentum (GeV)',
            'yaxtitle': 'Events',
            'treename': 'laurelin',
            'extrainfos': ['K^{0}_{S} candidates'],
            'bkgmodes': {
              'bkgdefault': {'type':None, 'info':'Background not subtracted'},
              'bkgsideband': {'type':'sideband', 'info': 'Background subtracted',
                              'sidevariable': '_mass',
                              'sidebins': np.linspace(0.44, 0.56, num=31, endpoint=True)},
            },
            'bins': {
              'finebins':np.linspace(0,100,num=161,endpoint=True),
              'defaultbins':[0.,0.5,1.5,4.,10.,20., 30., 40., 50., 60., 70., 80., 90., 100.],
            },
            'normalization': {
              'normlumi':{'type':'lumi', 'info':'Normalized to luminosity'},
              'normeventyield':{'type':'eventyield', 'info':'Normalized to data events'},
              'normrange':{'type':'range', 'info': 'Normalized for #Delta_{2D} < 0.5 cm', 
                          'normvariable': '_RPV', 'normrange':[0.,0.5]},
            },
    },
    'rpvsig': {'variablename':'_RPVSig',
               'xaxtitle':'#Delta_{2D} significance',
               'yaxtitle': 'Reconstructed vertices',
               'treename': 'laurelin',
               'extrainfos': ['K^{0}_{S} candidates'],
               'bkgmodes': {
                 'bkgsideband': {'type':'sideband', 'info': 'Background subtracted',
                                 'sidevariable': '_mass',
                              'sidebins': np.linspace(0.44, 0.56, num=31, endpoint=True)},
               },
               'bins': {
                 'finebins':np.linspace(0,600,num=61,endpoint=True),
               },
               'normalization': {
                 'normrange':{'type':'range', 'info': 'Normalized for #Delta_{2D} < 0.5 cm',
                             'normvariable': '_RPV', 'normrange':[0.,0.5]},
               },
    }
})
