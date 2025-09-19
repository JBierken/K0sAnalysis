import numpy as np

config = ({
    
  'rpv_vs_pt': {'variablename': '_RPV',
                'xaxtitle': '#Delta_{2D} (cm)',
		'yvariablename': '_pt',
                'yaxtitle': 'p_{T} (GeV)',
                'treename': 'laurelin',
                'extrainfos': ['K^{0}_{S} candidates'],
                'bkgmodes': {
                  'bkgsideband': {'type': 'sideband', 'info': 'Background subtracted',
                                  'sidevariable': '_mass',
                                  'sidebins': np.linspace(0.44, 0.56, num=31, endpoint=True)}
                },
                'bins': {
	               'defaultbins':{ 'xbins':[0.,0.5,1.5,4.5,20.,100.],
			          'ybins':[0.,5.,10.,20.,500.] },
                   'finexbins':{ 'xbins':[0.,0.5,1.,1.5,2.,5.,10.,15.,20.],
                                'ybins':[0.,5.,10.,20.] },
                },
                'normalization': {
	          'normrange':{ 'type':'range', 'info': 'Normalized for #Delta_{2D} < 0.5 cm',
                                'normvariable': '_RPV', 'normrange':[0.,0.5] },
	        }
  },
  'rpvsig_vs_pt': {'variablename': '_RPVSig', 
                'xaxtitle': '#Delta_{2D} significance',
                'yvariablename': '_pt', 
                'yaxtitle': 'p_{T} (GeV)',
                'treename': 'laurelin',
                'extrainfos': ['K^{0}_{S} candidates'],
                'bkgmodes': {
                  'bkgsideband': {'type': 'sideband', 'info': 'Background subtracted',
                                  'sidevariable': '_mass',
                                  'sidebins': np.linspace(0.44, 0.56, num=31, endpoint=True)}
                },
                'bins': {
	              'defaultbins':{ 'xbins':[0., 15., 50., 1500.],
                                'ybins':[0.,5.,10.,20., 500.] },
                  'finexbins':{ 'xbins':[0., 15., 50., 100., 200., 300., 400., 500., 600., 800.,1500.],
                                'ybins':[0.,5.,10.,20., 50., 500.] },
                },
                'normalization': {
                  'normrange':{ 'type':'range', 'info': 'Normalized for #Delta_{2D} < 0.5 cm',
                                'normvariable': '_RPV', 'normrange':[0.,0.5] },
                }
  }
  })

