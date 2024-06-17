############################################
# Read RPV variable and make distributions #
############################################
# With options for background subtraction, normalization, etc.


# import external modules
import os
import sys
import copy
import json
import argparse
import uproot
import numpy as np
import ROOT
from array import array
# import framework modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
from fitting.count_peak import count_peak_unbinned
from reweighting.pileup.pileupreweighter import PileupReweighter


if __name__=='__main__':

  sys.stderr.write('###starting###\n')

  # ------------------------------------------------------------------------
  # Get arguments
  # ------------------------------------------------------------------------
  parser = argparse.ArgumentParser( description = 'Fill histograms' )
  # general arguments
  parser.add_argument('-i', '--inputconfig', required=True, type=os.path.abspath)
  parser.add_argument('-t', '--treename', required=True)
  parser.add_argument('-v', '--variable', required=True, type=os.path.abspath)
  parser.add_argument('-o', '--outputfile', required=True)
  parser.add_argument('-n', '--nprocess', type=int, default=-1)
  # arguments for background subtraction
  parser.add_argument('--bkgmode', default=None, choices=[None, 'sideband'])
  parser.add_argument('--sidevariable', default=None, type=os.path.abspath)
  parser.add_argument('--sideplotdir', default=None)
  # arguments for normalization
  parser.add_argument('--normmode', default=None,
    choices=[None, 'lumi', 'yield', 'range', 'eventyield'])
  parser.add_argument('--normvariable', default=None, type=os.path.abspath)
  parser.add_argument('--eventtreename', default=None)
  # arguments for secondary binning
  parser.add_argument('--yvariable', default=None, type=os.path.abspath)
  args = parser.parse_args()

  # ------------------------------------------------------------------------
  # load input configuration
  # ------------------------------------------------------------------------
  with open(args.inputconfig) as f:
    temp        = json.load(f)
  datain        = temp['datain']
  simin         = temp['mcin']
  ndatafiles    = len(datain)
  nsimfiles     = len(simin)

  # prints for checking
  print('Found following input configuration:')
  for i,datadict in enumerate(datain):
    print('  - Data file {}'.format(i+1))
    for key,val in datadict.items():
      print('    {}: {}'.format(key, val))
  for i,simdict in enumerate(simin):
    print('  - Simulation file {}'.format(i+1))
    for key,val in simdict.items():
      print('    {}: {}'.format(key, val))

  # check if all input files exist
  missing = []
  for sample in simin + datain:
    if not os.path.exists(sample['file']):
      missing.append(sample['file'])
  if len(missing)>0:
    msg = 'ERROR: following input files do not seem to exist:\n'
    for f in missing: msg += '  - {}\n'.format(f)
    raise Exception(msg)

  # check luminosity
  totallumi = sum([sample['luminosity'] for sample in simin])
  lumitest  = sum([sample['luminosity'] for sample in datain])
  if( abs(lumitest-totallumi)/float(totallumi)>0.001 ):
    print('WARNING: total luminosity for data and simulation do not agree!')
    print(' (luminosity values for data are only used for plot labels;')
    print(' the values for simulations are used in event weighting and to calculate the sum)')

  # ------------------------------------------------------------------------
  # load variable
  # ------------------------------------------------------------------------
  # load main variable
  with open(args.variable) as f:
    variable = json.load(f)
  print('Found following main variable:')
  for key,val in variable.items(): print('  {}: {}'.format(key,val))

  # load sideband variable
  sidevariable = None # default case if no background subtraction
  if args.bkgmode in ['sideband']:
    if args.sidevariable is None:
      msg = 'ERROR: requested background subtraction via sideband,'
      msg += ' but sideband variable was not specified.'
      raise Exception(msg)
    with open(args.sidevariable) as f:
      sidevariable = json.load(f)
    print('Found following sideband variable:')
    for key,val in sidevariable.items(): print('  {}: {}'.format(key,val))

  # load normalization variable
  normvariable = None # default case if no normalization in range
  if args.normmode in ['range']:
    if args.normvariable is None:
      msg = 'ERROR: requested normalization in range,'
      msg += ' but normalization variable was not specified.'
      raise Exception(msg)
    with open(args.normvariable) as f:
      normvariable = json.load(f)
    print('Found following normalization variable:')
    for key,val in normvariable.items(): print('  {}: {}'.format(key,val))
    # do check on binning of normalization variable
    if len(normvariable['bins'])>2:
      msg = 'WARNING: found more than one bin for normalization variable'
      msg += ' but only min and max are taken into account for normalization range;'
      msg += ' intermediate bin edges are ignored.'
      print(msg)
      normvariable['bins'] = [normvariable['bins'][0], normvariable['bins'][-1]]

  # load secondary variable
  yvariable = None # default case if no secondary variable
  if args.yvariable is not None:
    with open(args.yvariable) as f:
      yvariable = json.load(f)
    print('Found following secondary variable:')
    for key,val in yvariable.items(): print('  {}: {}'.format(key,val))

  # set luminosity and xsection for simulation to 1 if no lumi scaling is requested
  if args.normmode is None:
    for simdict in simin:
      simdict['luminosity'] = 1
      simdict['xsection']   = 1

  # ------------------------------------------------------------------------
  # define help function to process a single file:
  #
  #     input:  - single input file (+ tree name within that file)
  #             - variable (name in the tree + binning)
  #     output: - histogram in the form of two numpy arrays
  #               holding the counts and corresponding errors 
  #               of variable values in chosen binning,
  #               potentially after background subtraction if requested;
  #             - the resulting arrays are two-dimensional 
  #               if a secondary variable is provided.
  #     note:   - the counts are normalized to the provided xsection and lumi
  #               (if isdata is False, else each value simply gets weight 1)
  #     note:   - if variable is None, return sum of weights and corresponding error
  # ------------------------------------------------------------------------
  def get_histogram(inputfile, 
                    treename, 
                    variable        =None, 
                    isdata          =False,
                    yvariable       =None,
                    xsection        =1, 
                    lumi            =1, 
                    weightvarname   ='_weight',
                    hcountername    ='hCounter',
                    sidevariable    =None,
                    label           =None,
                    nentries        =None,
                    year            =None,          # for reweighter
                    campaign        =None           # for reweighter 
                   ):
    # open the file and read hcounter
    print('Now running on file {}...'.format(inputfile))
    with uproot.open(inputfile) as f:
      sumweights     = 1             # default case for data, overwritten for simulation below
      prescale       = None          # to implement later
      if not isdata:
        try: sumweights = f[hcountername].values()[0]
        except:
          msg = 'WARNING: isdata was set to False, but no valid hCounter found in file'
          msg += ' (for provided key {}),'.format(hcountername)
          msg += ' will use sum of weights = 1 for this sample.'
          msg += ' Valid keys are {}'.format(f.keys())
          print(msg)
      
      # get main tree and manage number of entries
      tree              = f[treename]
      nentries_reweight = 1.
      
      if( nentries is not None and nentries>0 and nentries<tree.num_entries ):
        nentries_reweight = tree.num_entries / nentries
      else: nentries = tree.num_entries
      msg   =   'Tree {} was found to have {} entries,'.format(                     treename, tree.num_entries)
      msg   +=  ' of which {} will be read (using reweighting factor {}).'.format(  nentries, nentries_reweight)
      print(msg)

      # get weights
      if isdata: weights    = np.ones(nentries)
      else:
        weights             = tree[weightvarname].array(library='np', entry_stop=nentries)
        weights             = weights / sumweights * xsection * lumi
      weights               = weights * nentries_reweight

      # do reweighting
      if not isdata:
        print('Doing pileup reweighting...')
        pileupreweighter    = PileupReweighter(campaign, year)
        pileupreweighter.initsample(inputfile)
        ntrueint            = tree['_nTrueInt'].array(library='np', entry_stop=nentries)
        pileupreweight      = pileupreweighter.getreweight(ntrueint)
        weights             = np.multiply(weights, pileupreweight)

      # if no variable was specified, return sum of weights
      if variable is None:
        sumweights          = np.sum(weights)
        error               = np.sqrt(np.sum(np.power(weights, 2)))
        
        return (sumweights, error)

      # get the variable and some masks
      varvalues             = tree[variable['variable']].array(library='np', entry_stop=nentries)
      nanmask               = np.isnan(varvalues)
      rangemask             = ((varvalues > variable['bins'][0]) & (varvalues < variable['bins'][-1]))
      totalmask             = ((~nanmask) & rangemask)

      # initialize a dummy secondary variable if it was not provided
      # (easier than if else statements below)
      dim = 1
      if yvariable is not None: dim = 2
      else:
        dummyvar            = tree.keys()[0]
        dummyvalues         = tree[dummyvar].array(library='np', entry_stop=nentries)
        dummymin            = np.min(dummyvalues)
        dummymax            = np.max(dummyvalues)
        yvariable           = {'variable': dummyvar, 'bins': [dummymin/2., dummymax*2]}

      # get the secondary variable
      yvarvalues            = tree[yvariable['variable']].array(library='np', entry_stop=nentries)
      ynanmask              = np.isnan(yvarvalues)
      yrangemask            = ((yvarvalues >= yvariable['bins'][0]) & (yvarvalues <= yvariable['bins'][-1]))
      totalmask             = (totalmask & (~ynanmask) & yrangemask)

      # case of no background subtraction
      if sidevariable is None:
        varvalues           = varvalues[totalmask]
        yvarvalues          = yvarvalues[totalmask]
        weights             = weights[totalmask]
        counts              = np.histogram2d(varvalues, yvarvalues,
                                bins=(variable['bins'], yvariable['bins']),
                                weights=weights)[0]
        confidence          = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        errors              = np.sqrt(np.histogram2d(varvalues, yvarvalues,
                                bins=(variable['bins'], yvariable['bins']), 
                                weights=np.power(weights,2))[0])
       
      # do background subtraction
      if sidevariable is not None:
        # get values of sideband variable
        sidebandvalues          = tree[sidevariable['variable']].array(library='np', entry_stop=nentries)
        # initialize final histograms
        counts                  = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        confidence              = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        errors                  = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        confidence_error        = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        # loop over main variable bins and secondary variable bins
        for i, (low, high) in enumerate(zip(variable['bins'][:-1], variable['bins'][1:])):
          for j, (ylow, yhigh) in enumerate(zip(yvariable['bins'][:-1], yvariable['bins'][1:])):
            # subselect sideband variable values in this main variable bin
            mask                = ((varvalues > low) & (varvalues < high)
                                    & (yvarvalues > ylow) & (yvarvalues < yhigh))
            thissidebandvalues  = sidebandvalues[mask]
            thisweights         = weights[mask]
            # make extra info
            extrainfo           = '{0:.2f} < '.format(low)
            extrainfo           += variable['label']
            extrainfo           += ' < {0:.2f}'.format(high)
            if dim==2:
              extrainfo         += '<< {0:.2f} < '.format(ylow)
              extrainfo         += yvariable['label']
              extrainfo         += ' < {0:.2f}'.format(yhigh)
            
            # fit background and count what is left in peak
            histlabel = 'Data' if isdata else 'Simulation'
            histname                = '{}_bin{}'.format(label, i)
            if dim==2: histname     += '_ybin{}'.format(j)
            (npeak, nerror, conf, conf_error)   = count_peak_unbinned(
                                        thissidebandvalues,
                                        thisweights,
                                        sidevariable,
                                        mode            ='hybrid',
                                        label           =histlabel,
                                        lumi            =lumi,
                                        extrainfo       =extrainfo,
                                        histname        =histname,
                                        plotdir         =args.sideplotdir
                                      )
            
            counts[i,j]             = npeak
            confidence[i,j]         = conf
            errors[i,j]             = nerror
            confidence_error[i,j]   = conf_error
        
        # Calculate the error on the confidence method
        #conf_error      = np.zeros((len(variable['bins'])-1, len(yvariable['bins'])-1))
        #for i, (low, high) in enumerate(zip(variable['bins'][:-1], variable['bins'][1:])):
        #  for j, (ylow, yhigh) in enumerate(zip(yvariable['bins'][:-1], yvariable['bins'][1:])):
        #      conf_error[i,j] = np.sqrt(np.power(confidence[i,j] - np.average(confidence), 2))

       
    # remove superfluous dimension for one-dimensional arrays
    if dim==1:
      counts = counts[:,0]
      errors = errors[:,0]

    # return histogram with counts and corresponding errors
    #return (counts, errors, confidence)
    return (counts, errors, confidence, confidence_error)
  
  # ------------------------------------------------------------------------
  # loop over input files and fill histograms
  # ------------------------------------------------------------------------
  # Data files
  for datadict in datain:
    print('Now running on data file {}...'.format(datadict['file']))
    lumi                        = datadict['luminosity']
    #counts, errors, confidences = get_histogram(
    counts, errors, confidences, conf_errors = get_histogram(
                                                datadict['file'], 
                                                args.treename,
                                                variable        =variable,
                                                yvariable       =yvariable,
                                                isdata          =True, 
                                                lumi            =lumi,
                                                sidevariable    =sidevariable,
                                                label           =datadict['label'].strip(' .'),
                                                nentries        =args.nprocess
                                            )
    datadict['counts']          = counts
    datadict['errors']          = errors
    datadict['confidences']     = confidences
    datadict['conf_errors']     = conf_errors    
  # Simulation files
  for simdict in simin:
    print('Now running on simulation file {}...'.format(simdict['file']))
    xsection                    = simdict['xsection']
    lumi                        = simdict['luminosity']
    counts, errors, confidences, conf_errors = get_histogram(
                                      simdict['file'], 
                                      args.treename,
                                      variable        =variable,
                                      yvariable       =yvariable,
                                      isdata          =False, 
                                      xsection        =xsection, 
                                      lumi            =lumi,
                                      sidevariable    =sidevariable,
                                      label           =simdict['label'].strip(' .'),
                                      nentries        =args.nprocess,
                                      year            =simdict['year'], 
                                      campaign        =simdict['campaign']
                                  )
    simdict['counts']           = counts
    simdict['errors']           = errors
    simdict['confidences']      = confidences
    simdict['conf_errors']      = conf_errors    

  # clip histograms to minimum zero
  for datadict in datain:
    datadict['counts']          = np.clip(datadict['counts'],           0, None)
    datadict['errors']          = np.clip(datadict['errors'],           0, None)
    datadict['confidences']     = np.clip(datadict['confidences'],      0, None)
    datadict['conf_errors']     = np.clip(datadict['conf_errors'],      0, None)
  for simdict in simin:
    simdict['counts']           = np.clip(simdict['counts'],            0, None)
    simdict['errors']           = np.clip(simdict['errors'],            0, None)
    simdict['confidences']      = np.clip(simdict['confidences'],       0, None)
    simdict['conf_errors']      = np.clip(simdict['conf_errors'],       0, None)

  # ------------------------------------------------------------------------
  # Apply Normalization:
  #     - now need to manage additional normalization of histograms.
  #       for normmode=None and normmode='lumi', no additional steps are needed;
  #     - other cases: treated below
  # ------------------------------------------------------------------------
  # for normmode 'yield', normalize sum of simulation to sum of data
  if args.normmode=='yield':
    print('Normalizing simulation yield to data yield...')
    
    # Get Simulation count
    simsum                          = 0
    for simdict in simin: simsum    += np.sum(simdict['counts'])
    
    # Get Data count
    datasum                         = 0
    for datadict in datain: datasum += np.sum(datadict['counts'])
    
    # Compare data vs. sim
    scale                           = datasum / simsum
    
    # Scale simulation to sum of data
    for simdict in simin:
      simdict['counts']             = simdict['counts']*scale
      simdict['errors']             = simdict['errors']*scale
      simdict['confidences']        = simdict['confidences']*scale
      simdict['conf_errors']        = simdict['conf_errors']*scale

  # ------------------------------------------------------------------------
  # for normmode 'range', normalize sum of simulation to sum of data,
  # but the sum is calculated only for a given variable in given range
  if args.normmode=='range':
    print('Normalizing simulation yield to data yield in range...')
    
    # Get Data count
    datasum     = 0
    for datadict in datain:
      counts, _, confs, c_ = get_histogram(
                       datadict['file'], 
                       args.treename,
                       variable    =normvariable,
                       isdata      =True,
                       label       =datadict['label'].strip(' .')+'_normrange',
                       sidevariable=sidevariable,
                       nentries    =args.nprocess
                     )
      if len(counts)!=1:
        msg = 'ERROR: counts has unexpected length, check the binning of normvariable.'
        raise Exception(msg)
      
      datasum  += counts[0]
    
    # Get Simulation count
    simsum      = 0
    for simdict in simin:
      xsection     = simdict['xsection']
      lumi         = simdict['luminosity']
      counts, _, confsy, c_ = get_histogram(
                       simdict['file'], 
                       args.treename,
                       variable    =normvariable,
                       isdata      =False, 
                       xsection    =xsection, 
                       lumi        =lumi,
                       label       =simdict['label'].strip(' .')+'_normrange',
                       sidevariable=sidevariable,
                       nentries    =args.nprocess,
                       year        =simdict['year'], 
                       campaign    =simdict['campaign']
                     )
      if len(counts)!=1:
        msg = 'ERROR: counts has unexpected length, check the binning of normvariable.'
        raise Exception(msg)
      
      simsum  += counts[0]
    
    # Compare data vs. sim
    scale   = datasum / simsum
    
    # Scale simulation to sum of data
    for simdict in simin:
      simdict['counts']             = simdict['counts']*scale
      simdict['errors']             = simdict['errors']*scale
      simdict['confidences']        = simdict['confidences']*scale
      simdict['conf_errors']        = simdict['conf_errors']*scale

  # ------------------------------------------------------------------------
  # for normmode 'eventyield', scale using event weights
  if args.normmode=='eventyield':
    print('Normalizing simulation event yield to data event yield...')
    if args.eventtreename is None:
        msg = 'ERROR: requested normalization by event yield, but event tree name was not specified.'
        raise Exception(msg)
    
    # Get Data count
    datasum     = 0
    for datadict in datain:
      sumweights, _, confs, c_ = get_histogram(
                             datadict['file'], 
                            args.eventtreename, 
                            isdata    =True, 
                            nentries  =args.nprocess
                        )
       
      datasum  += sumweights
    
    # Get Simulation count
    simsum      = 0
    for simdict in simin:
      xsection      = simdict['xsection']
      lumi          = simdict['luminosity']
      sumweights, _, confs, c_ = get_histogram(
                            simdict['file'], 
                            args.eventtreename,
                            isdata      =False, 
                            xsection    =xsection, 
                            lumi        =lumi,
                            nentries    =args.nprocess,
                            year        =simdict['year'], 
                            campaign    =simdict['campaign']
                         )
      
      simsum  += sumweights
    
    # Compare data vs. sim
    scale = datasum / simsum
    
    # Scale simulation to sum of data
    for simdict in simin:
      simdict['counts']             = simdict['counts']*scale
      simdict['errors']             = simdict['errors']*scale
      simdict['confidences']        = simdict['confidences']*scale
      simdict['conf_errors']        = simdict['conf_errors']*scale

  # ------------------------------------------------------------------------
  # write histograms and meta-info to file:
  #     note:   - for now this is done with PyROOT instead of uproot
  # ------------------------------------------------------------------------
  print('Writing histograms to file...')

  # write histogram for integral count
  f = ROOT.TFile.Open(args.outputfile, "recreate")
  # write histograms
  for ddict in simin + datain:
    counts      = ddict['counts']
    errors      = ddict['errors']
    
    # conversion to ROOT.TH1
    if(len(counts.shape)==1):
      hist = ROOT.TH1F(
              ddict['label'], 
              ddict['label'], 
              len(variable['bins'])-1,
              array('f', variable['bins'])
             )
      for i, (count, error) in enumerate(zip(counts, errors)):
        hist.SetBinContent( i+1, count)
        hist.SetBinError(   i+1, error)
    
    # converstion to ROOT.TH2
    elif(len(counts.shape)==2):
      hist = ROOT.TH2F(
              ddict['label'], 
              ddict['label'],
              len(variable['bins'])-1, 
              array('f', variable['bins']),
              len(yvariable['bins'])-1, 
              array('f', yvariable['bins'])
             )
      for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
          hist.SetBinContent(   i+1, j+1, counts[i,j])
          hist.SetBinError(     i+1, j+1, errors[i,j])
          
    else:
      msg = 'ERROR: shape of counts array could not be converted to TH1 or TH2.'
      raise Exception(msg)
    hist.Write(hist.GetName())
  
  # write variable
  varname_st        = ROOT.TNamed('variable', variable['name'])
  varname_st.Write()
  
  # write secondary variable
  if yvariable is not None:
    yvarname_st     = ROOT.TNamed('yvariable', yvariable['name'])
    yvarname_st.Write()
  
  # write normalization
  normalization_st  = ROOT.TNamed('normalization', str(args.normmode))
  normalization_st.Write()
  
  # write norm range
  if args.normmode in ['range']:
    normrange_st    = ROOT.TVectorD(2)
    normrange_st[0] = normvariable['bins'][0]
    normrange_st[1] = normvariable['bins'][1]
    normrange_st.Write("normrange")
    normvariable_st = ROOT.TNamed('normvariable', normvariable['name'])
    normvariable_st.Write()
  
  # write luminosity
  lumi_st           = ROOT.TVectorD(1)
  lumi_st[0]        = totallumi
  lumi_st.Write("lumi")
  
  # write background mode
  bkgmode_st        = ROOT.TNamed('bkgmode', str(args.bkgmode))
  bkgmode_st.Write()
  
  # write tree name
  treename_st       = ROOT.TNamed('treename', str(args.treename))
  treename_st.Write()
  f.Close()

  # write histogram for confidence count
  outputname = args.outputfile.split(".")
  outputname = outputname[0] + "_confidence.root"
  f2 = ROOT.TFile.Open(outputname, "recreate")
  # write histograms
  for ddict in simin + datain:
    confidences = ddict['confidences']
    errors      = ddict['conf_errors']
    
    # conversion to ROOT.TH1
    if(len(confidences.shape)==1):
      hist2 = ROOT.TH1F(
              ddict['label'] + " confidence", 
              ddict['label'] + " confidence", 
              len(variable['bins'])-1,
              array('f', variable['bins'])
             )
      
      for i, (confidence, error) in enumerate(zip(confidences, errors)):
        hist2.SetBinContent( i+1, confidence)
        hist2.SetBinError(   i+1, error)
    
    # converstion to ROOT.TH2
    elif(len(confidences.shape)==2):
      hist2 = ROOT.TH2F(
              ddict['label'], 
              ddict['label'], 
              len(variable['bins'])-1, 
              array('f', variable['bins']),
              len(yvariable['bins'])-1, 
              array('f', yvariable['bins'])
             )
      for i in range(confidences.shape[0]):
        for j in range(confidences.shape[1]):
          hist2.SetBinContent(   i+1, j+1, confidences[i,j])
          hist2.SetBinError(     i+1, j+1, errors[i,j])
          
    else:
      msg = 'ERROR: shape of counts array could not be converted to TH1 or TH2.'
      raise Exception(msg)
    hist2.Write(hist2.GetName())
  
  # write variable
  varname_st        = ROOT.TNamed('variable', variable['name'])
  varname_st.Write()
  
  # write secondary variable
  if yvariable is not None:
    yvarname_st     = ROOT.TNamed('yvariable', yvariable['name'])
    yvarname_st.Write()
  
  # write normalization
  normalization_st  = ROOT.TNamed('normalization', str(args.normmode))
  normalization_st.Write()
  
  # write norm range
  if args.normmode in ['range']:
    normrange_st    = ROOT.TVectorD(2)
    normrange_st[0] = normvariable['bins'][0]
    normrange_st[1] = normvariable['bins'][1]
    normrange_st.Write("normrange")
    normvariable_st = ROOT.TNamed('normvariable', normvariable['name'])
    normvariable_st.Write()
  
  # write luminosity
  lumi_st           = ROOT.TVectorD(1)
  lumi_st[0]        = totallumi
  lumi_st.Write("lumi")
  
  # write background mode
  bkgmode_st        = ROOT.TNamed('bkgmode', str(args.bkgmode))
  bkgmode_st.Write()
  
  # write tree name
  treename_st       = ROOT.TNamed('treename', str(args.treename))
  treename_st.Write()
  f2.Close()
sys.stderr.write('###done###\n')
