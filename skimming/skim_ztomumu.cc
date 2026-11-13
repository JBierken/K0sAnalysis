/*
Perform Z -> mu+mu- skimming.
Runs on a single input file and produces a single output file.
Apart from event selection, branch reduction is performed.
*/

// include c++ library classes 
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>

// include ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVector2.h"
#include "Math/Vector4D.h"

// definition of constants
const static double     ZMass = 91.19;

// definition of b-tag thresholds
// see e.g. subpages of https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
//              ---> Outdated: New Twiki: https://btv-wiki.docs.cern.ch/ 
// the thresholds correspond to loose working point for Deep-CSV (for Run-2) and DeepJet (for Run-3).
double bTagThreshold(   const std::string&  campaign, 
                        const std::string&  year)
{
    if(campaign     == "run2preul")
    {
        // Run-2: Loose working point for Deep-CSV
	    if(     year== "2016")          return 0.2217;
	    else if(year== "2017")          return 0.1522;
	    else if(year== "2018")          return 0.1241;
    } 
    else if(campaign=="run2ul")
    {
        // Run-2: Loose working point for Deep-CSV
	    if(     year=="2016PreVFP")    return 0.2027;
	    else if(year=="2016PostVFP")   return 0.1918;
        else if(year=="2017")          return 0.1355;
        else if(year=="2018")          return 0.1208;
    }
    else if (campaign=="run3")
    {
        // Run-3: Loose working point for DeepJet (applied to _jetDeepFlavor = b+bb+lepb, AK4 PUPPI jets)
        // 2022-2023BPix: official values, checked against BTV jsonpog-integration "deepJet_wp_values" (Nov 2026)
	    if(     year=="2022")          return 0.0583;
	    else if(year=="2022EE")        return 0.0614;
	    else if(year=="2023")          return 0.0479;
	    else if(year=="2023BPix")      return 0.048;
	    else if(year=="2024")          return 0.048; // PLACEHOLDER (2023BPix value): no official 2024 DeepJet WP (2024_Summer24 only has UParT)
	    else if(year=="2025")          return 0.048; // PLACEHOLDER (2023BPix value): no official 2025 DeepJet WP yet
        //else if(year=="2024")          return 0.0246;
	    //else if(year=="2025")          return 0.0246; // SF not yet available for 2025???
    }
    // unknown campaign, or unknown year within a known campaign:
    // fail loudly instead of silently returning a threshold of 0 (which would tag every jet)
    std::string msg                     = "ERROR in bTagThreshold:";
    msg                                 += " campaign '" + campaign + "' or year '" + year + "' not recognized.";
    throw std::runtime_error(msg);
}

// ---------------------------------------------------------------------------------------------
// b-taggers stored in the V0Analyzer ntuples (AK4 PUPPI jets, see V0Analysis/V0Analyzer JetAnalyzer)
//   score branch = the b-vs-all discriminant on which the BTV working points are defined
//   (-1 in the ntuples means: tagger not available in this production, or jet not evaluated)
// ---------------------------------------------------------------------------------------------
enum Tagger { kDeepJet = 0, kPNet, kRobustParT, kUParT, nTaggers };
const std::string taggerNames[nTaggers]    = { "DeepJet",        "PNet",                   "RobustParT",            "UParT"            };
const std::string taggerBranches[nTaggers] = { "_jetDeepFlavor", "_jetParticleNet_BvsAll", "_jetRobustParT_BvsAll", "_jetUParT_BvsAll" };

Tagger taggerFromName(const std::string& name)
{
    for(int t = 0; t < nTaggers; ++t){
        std::string lower = taggerNames[t];
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if(name == taggerNames[t] || name == lower) return static_cast<Tagger>(t);
    }
    throw std::runtime_error("ERROR: unknown b-tagger '" + name + "' (use DeepJet, PNet, RobustParT or UParT)");
}

// Loose working points per tagger and year. Returns -1 if there is no working point.
// Source: BTV jsonpog-integration (cvmfs: /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/<era>/btagging.json.gz,
//         corrections "<tagger>_wp_values"), checked Nov 2026:
//   - 2022, 2022EE, 2023, 2023BPix: DeepJet, ParticleNet, RobustParT
//   - 2024 (Summer24): UParT only
//   - 2025: nothing yet
// Every value without an official BTV number is marked "PLACEHOLDER" (grep for it to update later):
//   DeepJet 2024/2025, PNet 2024/2025, RobustParT 2024/2025 -> 2023BPix value
//   UParT 2022-2023BPix -> 2024 value, UParT 2025 -> 2024 value
// Note: RobustParT is not stored in 2024 MINIv6/2025 ntuples and UParT not in 2022/2023 ntuples;
//       there the tagger is detected as not evaluated and _nimloth_nBJets_<tagger> = -1 regardless of the WP.
double looseWP(Tagger tagger, const std::string& campaign, const std::string& year)
{
    if(tagger == kDeepJet) return bTagThreshold(campaign, year);
    if(campaign != "run3") return -1.;
    if(tagger == kPNet){
        if(     year == "2022")         return 0.047;
        else if(year == "2022EE")       return 0.0499;
        else if(year == "2023")         return 0.0358;
        else if(year == "2023BPix")     return 0.0359;
        else if(year == "2024")         return 0.0359;   // PLACEHOLDER (2023BPix value): no official 2024 WP yet
        else if(year == "2025")         return 0.0359;   // PLACEHOLDER (2023BPix value): no official 2025 WP yet
    }
    else if(tagger == kRobustParT){
        if(     year == "2022")         return 0.0849;
        else if(year == "2022EE")       return 0.0897;
        else if(year == "2023")         return 0.0681;
        else if(year == "2023BPix")     return 0.0683;
        else if(year == "2024")         return 0.0683;   // PLACEHOLDER (2023BPix value): no official 2024 WP yet
        else if(year == "2025")         return 0.0683;   // PLACEHOLDER (2023BPix value): no official 2025 WP yet
    }
    else if(tagger == kUParT){
        if(     year == "2022")         return 0.0246;   // PLACEHOLDER (2024 value): no official 2022 WP
        else if(year == "2022EE")       return 0.0246;   // PLACEHOLDER (2024 value): no official 2022EE WP
        else if(year == "2023")         return 0.0246;   // PLACEHOLDER (2024 value): no official 2023 WP
        else if(year == "2023BPix")     return 0.0246;   // PLACEHOLDER (2024 value): no official 2023BPix WP
        else if(year == "2024")         return 0.0246;
        else if(year == "2025")         return 0.0246;   // PLACEHOLDER (2024 value): no official 2025 WP yet
    }
    return -1.;
}

bool containsAny(const std::string& s, const std::vector<std::string>& patterns){
    for(const auto& p : patterns) if(s.find(p) != std::string::npos) return true;
    return false;
}

std::pair<std::string, std::string> getYearFromFileName(const std::string& fileName){
    // small help function to exctract the campaign and year from the name of an input file

    // Run-2: Ultra-Legacy (UL) campaign
    if( fileName.find(      "RunIISummer20UL16MiniAODAPVv2")    != std::string::npos 
        || fileName.find(   "HIPM_UL2016_MiniAODv2")            != std::string::npos )            
    { 
        return std::make_pair("run2ul",     "2016PreVFP"); 
    }
    else if( fileName.find( "RunIISummer20UL16MiniAODv2")       != std::string::npos
        || fileName.find(   "UL2016_MiniAODv2")                 != std::string::npos )                 
    { 
        return std::make_pair("run2ul",     "2016PostVFP"); 
    }
    else if( fileName.find( "RunIISummer20UL17MiniAODv2")       != std::string::npos
        || fileName.find(   "UL2017_MiniAODv2")                 != std::string::npos )                   
    { 
        return std::make_pair("run2ul",     "2017"); 
    }
    else if( fileName.find( "RunIISummer20UL18MiniAODv2")       != std::string::npos
        || fileName.find(   "UL2018_MiniAODv2")                 != std::string::npos )                   
    { 
        return std::make_pair("run2ul",     "2018"); 
    }
    // Run-2: pre Ultra-Legacy (pre-UL) campaign
    else if( fileName.find( "Run2016")                          != std::string::npos 
	|| fileName.find(       "RunIISummer16")                    != std::string::npos )
    { 
        return std::make_pair("run2preul",  "2016"); 
    }
    else if( fileName.find( "Run2017")                          != std::string::npos 
	|| fileName.find(       "RunIIFall17")                      != std::string::npos )
    { 
        return std::make_pair("run2preul",  "2017"); 
    }
    else if( fileName.find( "Run2018")                          != std::string::npos 
	|| fileName.find(       "RunIIAutumn18")                    != std::string::npos )
    {	
        return std::make_pair("run2preul",  "2018"); 
    }
    // Run-3
    // note: the more specific sub-campaigns (2022EE, 2023BPix) must be checked
    //       before their parent year, since e.g. "Run3Summer22EE" contains "Run3Summer22".
    //       Data is matched on the era in the dataset name (Run2022E/F/G = postEE, Run2023D = postBPix),
    //       simulation on the MC campaign name in the dataset path.
    else if( containsAny(fileName, {"Run2022E", "Run2022F", "Run2022G", "Run3Summer22EE", "Run3Winter22EE"}) )
    {
        return std::make_pair("run3",       "2022EE");
    }
    else if( containsAny(fileName, {"Run2022", "Run3Summer22", "Run3Winter22"}) )
    {
        return std::make_pair("run3",       "2022");
    }
    else if( containsAny(fileName, {"Run2023D", "Run3Summer23BPix", "Run3Winter23BPix"}) )
    {
        return std::make_pair("run3",       "2023BPix");
    }
    else if( containsAny(fileName, {"Run2023", "Run3Summer23", "Run3Winter23"}) )
    {
        return std::make_pair("run3",       "2023");
    }
    else if( containsAny(fileName, {"Run2024", "RunIII2024Summer24", "Run3Summer24", "Run3Winter24"}) )
    {
        return std::make_pair("run3",       "2024");
    }
    else if( containsAny(fileName, {"Run2025", "RunIII2025", "Run3Summer25", "Run3Winter25"}) )
    {
        return std::make_pair("run3",       "2025");
    }
    else{
	    std::string msg     = "ERROR: year could not be extracted from filename";
	    msg                 += " " + fileName;
	    throw std::runtime_error(msg);
    }
} 

std::vector< std::shared_ptr< TH1 > > getHistogramsFromFile( TFile* f )
{
    // get histograms from root file; 
    // note that the file is assumed to be opened and closed outside this function!

    // vector containing all histograms in current file
    std::vector< std::shared_ptr< TH1 > > histogramVector;
    // loop over keys in blackJackAndHookers directory
    TDirectory* dir                     = (TDirectory*) f->Get("blackJackAndHookers");
    TList* keyList                      = dir->GetListOfKeys();
    for( const auto objectPtr : *keyList )
    {
	    // try if a dynamic_cast to a histogram works to check if object is histogram
	    TH1* histPtr                    = dynamic_cast< TH1* >( dir->Get( objectPtr->GetName() ) );
	    if( histPtr )
        {
            // make sure histograms don't get deleted by root upon deletion of TDirectory above
            histPtr->SetDirectory(      gROOT );
	        histogramVector.emplace_back(histPtr );
	    }
    }
    return histogramVector;
}

void skimFile(const std::string&        inputFilePath,
		        const std::string&      outputFilePath,
		        long unsigned           nEntries,
                Tagger                  vetoTagger )
{

    std::cout <<    "skimming "         <<  inputFilePath   <<  std::endl;

    std::pair<std::string, std::string> temp = getYearFromFileName(inputFilePath);
    std::string campaign                = temp.first;
    std::string year                    = temp.second;

    // load input histograms and tree
    TFile* inputFilePtr                 = TFile::Open( inputFilePath.c_str() );
    std::vector< std::shared_ptr< TH1 > > histVector = getHistogramsFromFile( inputFilePtr );
    TTree* inputTreePtr                 = (TTree*) inputFilePtr->Get("blackJackAndHookers/blackJackAndHookersTree");
    
    // disable unneeded branches to save disk space
    inputTreePtr->SetBranchStatus(      "_gen_*",           0);
    inputTreePtr->SetBranchStatus(      "_Flag_*",          0);
    inputTreePtr->SetBranchStatus(      "_HLT_*",           0);
    inputTreePtr->SetBranchStatus(      "_tau*",            0);
    inputTreePtr->SetBranchStatus(      "_ph*",             0);
    inputTreePtr->SetBranchStatus(      "_jet*",            0);
    // note: every branch read in the event loop below must be (re-)enabled here,
    //       since GetEntry does not fill disabled branches (their buffers keep stale values).
    inputTreePtr->SetBranchStatus(      "_jetPt",           1);
    inputTreePtr->SetBranchStatus(      "_jetEta",          1);
    inputTreePtr->SetBranchStatus(      "_jetPhi",          1);
    // b-tag scores (kept in the output tree as well); a branch may be absent in older ntuples
    bool hasBranch[nTaggers];
    for(int t = 0; t < nTaggers; ++t){
        hasBranch[t]                    = ( inputTreePtr->GetBranch(taggerBranches[t].c_str()) != nullptr );
        if(hasBranch[t])                inputTreePtr->SetBranchStatus(taggerBranches[t].c_str(), 1);
    }
    inputTreePtr->SetBranchStatus(      "_met*",            0);

    // make output ROOT file
    TFile* outputFilePtr                = TFile::Open( outputFilePath.c_str() , "RECREATE" );
    outputFilePtr->mkdir(               "blackJackAndHookers" );
    outputFilePtr->cd(                  "blackJackAndHookers" );

    // write histograms to the new file
    for( const auto& histPtr : histVector ) histPtr->Write();

    // make output tree
    TTree* outputTreePtr                = inputTreePtr->CloneTree(0);
    // extend output tree with event variables
    Double_t _nimloth_Mll;              outputTreePtr->Branch("_nimloth_Mll",       &_nimloth_Mll,      "_nimloth_Mll/D");
    UInt_t _nimloth_nJets;              outputTreePtr->Branch("_nimloth_nJets",     &_nimloth_nJets,    "_nimloth_nJets/i");
    // extend output tree with muon variables (always exactly 2 muons per event)
    Double_t    _celeborn_lPt[2];       outputTreePtr->Branch("_celeborn_lPt",      &_celeborn_lPt,     "_celeborn_lPt[2]/D");
    Double_t    _celeborn_lEta[2];      outputTreePtr->Branch("_celeborn_lEta",     &_celeborn_lEta,    "_celeborn_lEta[2]/D");
    Double_t    _celeborn_lPhi[2];      outputTreePtr->Branch("_celeborn_lPhi",     &_celeborn_lPhi,    "_celeborn_lPhi[2]/D");
    UInt_t      _celeborn_lCharge[2];   outputTreePtr->Branch("_celeborn_lCharge",  &_celeborn_lCharge, "_celeborn_lCharge[2]/I");
    // extend output tree with the number of loose b-tagged jets per tagger (-1: tagger or working point not available)
    Int_t _nimloth_nBJets[nTaggers];
    for(int t = 0; t < nTaggers; ++t){
        const std::string name          = "_nimloth_nBJets_" + taggerNames[t];
        outputTreePtr->Branch(name.c_str(), &_nimloth_nBJets[t], (name + "/I").c_str());
    }

    // link input tree branches used to evaluate the skim condition to variables
    // NOTE: nL_max and nJets_max must correspond to the values used in the ntuplizer,
    //       else the skimmer could crash after the event loop with no clear error message... 
    //       would be better to write this as meta-info to the ntuples, 
    //       but not yet done as far as I know.
    //	     for old iteration of the kshort study: nL_max = 20, nJets_max = 20
    //       for new iteration of the kshort study: nL_max = 20, nJets_max = 100
    static const unsigned nL_max                    = 20;
    static const unsigned nJets_max                 = 100;
    UInt_t _nL;                                     inputTreePtr->SetBranchAddress("_nL",               &_nL);
    UInt_t _nLight;                                 inputTreePtr->SetBranchAddress("_nLight",           &_nLight);
    UInt_t _nMu;                                    inputTreePtr->SetBranchAddress("_nMu",              &_nMu);
    Bool_t _lPOGMedium[         nL_max];            inputTreePtr->SetBranchAddress("_lPOGMedium",       _lPOGMedium);
    Double_t _lPt[              nL_max];            inputTreePtr->SetBranchAddress("_lPt",              _lPt);
    Int_t _lCharge[             nL_max];            inputTreePtr->SetBranchAddress("_lCharge",          _lCharge);
    Bool_t _passTrigger_mm;                         inputTreePtr->SetBranchAddress("_passTrigger_mm",   &_passTrigger_mm);
    Double_t _lEta[             nL_max];            inputTreePtr->SetBranchAddress("_lEta",             _lEta);
    Double_t _lPhi[             nL_max];            inputTreePtr->SetBranchAddress("_lPhi",             _lPhi);
    Double_t _lE[               nL_max];            inputTreePtr->SetBranchAddress("_lE",               _lE);
    UInt_t _nJets;                                  inputTreePtr->SetBranchAddress("_nJets",            &_nJets);
    Double_t _jetPt[            nJets_max];         inputTreePtr->SetBranchAddress("_jetPt",            _jetPt);
    Double_t _jetEta[           nJets_max];         inputTreePtr->SetBranchAddress("_jetEta",           _jetEta);
    Double_t _jetPhi[           nJets_max];         inputTreePtr->SetBranchAddress("_jetPhi",           _jetPhi);
    // b-tag scores: the b-vs-all discriminant of each tagger (DeepJet: b+bb+lepb)
    Double_t _jetBTag[nTaggers][nJets_max];
    for(int t = 0; t < nTaggers; ++t){
        std::fill_n(_jetBTag[t], nJets_max, -1.);
        if(hasBranch[t])                            inputTreePtr->SetBranchAddress(taggerBranches[t].c_str(), _jetBTag[t]);
    }

    if(nEntries<=0)                                 nEntries        = inputTreePtr->GetEntries();

    // determine which taggers are usable in this file:
    //   branch present + evaluated (score >= 0) for at least one jet in the b-tag acceptance + a working point for this year.
    //   RobustParT is absent in 2024 MINIv6/2025 ntuples and UParT in 2022/2023 ntuples (stored as -1).
    bool evaluated[nTaggers]                        = {false, false, false, false};
    bool sawCentralJet                              = false;
    const long unsigned nScan                       = std::min<long unsigned>(nEntries, 2000);
    for(long unsigned entry = 0; entry < nScan; ++entry){
        inputTreePtr->GetEntry(entry);
        for(unsigned j = 0; j < _nJets; ++j){
            if(_jetPt[j] < 25. || std::fabs(_jetEta[j]) >= 2.4) continue;
            sawCentralJet                           = true;
            for(int t = 0; t < nTaggers; ++t) if(hasBranch[t] && _jetBTag[t][j] >= 0.) evaluated[t] = true;
        }
    }
    // tiny file without any jet in the b-tag acceptance: availability cannot be tested (and no jet can be b-tagged anyway)
    if(!sawCentralJet) for(int t = 0; t < nTaggers; ++t) evaluated[t] = hasBranch[t];
    bool   usable[nTaggers];
    double wp[nTaggers];
    std::cout << "b-taggers (veto with " << taggerNames[vetoTagger] << "):" << std::endl;
    for(int t = 0; t < nTaggers; ++t){
        wp[t]                                       = looseWP(static_cast<Tagger>(t), campaign, year);
        usable[t]                                   = hasBranch[t] && evaluated[t] && wp[t] > 0.;
        std::cout << "  " << taggerNames[t] << ": branch=" << (hasBranch[t] ? "yes" : "no")
                  << " evaluated=" << (evaluated[t] ? "yes" : "no") << " looseWP=" << wp[t]
                  << (usable[t] ? "" : "  -> not used, _nimloth_nBJets_" + taggerNames[t] + " = -1") << std::endl;
    }
    if(!usable[vetoTagger]){
        throw std::runtime_error("ERROR: b-veto tagger " + taggerNames[vetoTagger] + " is not usable for " + year
                                 + " (missing in the ntuple or no working point)");
    }

    // start event loop
    unsigned                                        debugcounter    = 0;

    for( long unsigned entry = 0; entry < nEntries; ++entry )
    {

	    inputTreePtr->GetEntry(entry);

        // apply event selection

	    // condition: at least two muons
	    if(_nMu < 2)                                continue;

        // make vector of muon indices, useful for indexing (light) lepton arrays further on
        // note: we assume that muons are filled first in the ntuples so indices run from 0 to _nMu
        std::vector<int>                            muoninds(_nMu);
        for(unsigned i = 0; i < _nMu; ++i)          muoninds[i] = i;
        
	    // condition: the number of muons with lPOGMedium must be exactly 2
        unsigned int                                i = 0;
        while(i < muoninds.size())
        {
            if(_lPOGMedium[muoninds[i]])            ++i;
            else                                    muoninds.erase(muoninds.begin()+i);
        }
        if(muoninds.size()!=2)                      continue;

        // order by _lPt
        if(_lPt[muoninds[0]] < _lPt[muoninds[1]])   std::swap(muoninds[0],muoninds[1]);

        // condition: explicit momentum and eta cuts
        if(_lPt[muoninds[0]] < 30. 
            || _lPt[muoninds[1]] < 25.)             continue;
        if(std::abs(_lEta[muoninds[0]]) > 2.4 
            || std::abs(_lEta[muoninds[1]]) > 2.4)  continue;

        // condition: opposite electric charge
        if(_lCharge[muoninds[0]]*_lCharge[muoninds[1]] > 0) continue;

        // set muon variables for output tree
        _celeborn_lPt[0]        = _lPt[muoninds[0]];        _celeborn_lPt[1]        = _lPt[muoninds[1]];
        _celeborn_lEta[0]       = _lEta[muoninds[0]];       _celeborn_lEta[1]       = _lEta[muoninds[1]];
        _celeborn_lPhi[0]       = _lPhi[muoninds[0]];       _celeborn_lPhi[1]       = _lPhi[muoninds[1]];
        _celeborn_lCharge[0]    = _lCharge[muoninds[0]];    _celeborn_lCharge[1]    = _lCharge[muoninds[1]];

        // condition: pass trigger
        if(!_passTrigger_mm)                        continue;

        // condition: dimuon invariant mass
        ROOT::Math::PtEtaPhiEVector vec1(_lPt[muoninds[0]],_lEta[muoninds[0]],
                        _lPhi[muoninds[0]],_lE[muoninds[0]]);
        ROOT::Math::PtEtaPhiEVector vec2(_lPt[muoninds[1]],_lEta[muoninds[1]],
                        _lPhi[muoninds[1]],_lE[muoninds[1]]);
        _nimloth_Mll                                = (vec1+vec2).M();
        if(std::abs(_nimloth_Mll-ZMass)>10.)        continue;

        // condition: b-jet veto
        // note: need to reconsider jet from lepton cleaning, maybe to strict here (e.g taus?)!
        // note: also need to consider jet quality cuts,
        //       especially eta requirement for b-tagging score to make sense
        int nbjets[nTaggers]                        = {0, 0, 0, 0};
        _nimloth_nJets                              = 0;
        for(unsigned j = 0; j < _nJets; ++j)
        {
            // jet pt
            if( _jetPt[j] < 25. )                   continue;

            // remove jets overlapping with leptons (delta phi wrapped into [-pi, pi])
	        bool                                    isleptonjet = false;
	        for(unsigned i = 0; i < _nL; ++i)
            {
                const double dEta                   = _lEta[i] - _jetEta[j];
                const double dPhi                   = TVector2::Phi_mpi_pi(_lPhi[i] - _jetPhi[j]);
		        if(std::sqrt(dEta*dEta + dPhi*dPhi) < 0.4) {isleptonjet = true; break;}
	        }
	        if(isleptonjet)                         continue;
	        ++_nimloth_nJets;

            // check if it is a b-jet (loose working point), separately for each usable tagger
	        bool inBTagAcceptance                   = ( std::fabs(_jetEta[j]) < 2.4 );
            if( !inBTagAcceptance )                 continue;
            for(int t = 0; t < nTaggers; ++t)
            {
                if( usable[t] && _jetBTag[t][j] > wp[t] ) ++nbjets[t];
            }
	    }
        for(int t = 0; t < nTaggers; ++t)           _nimloth_nBJets[t] = usable[t] ? nbjets[t] : -1;

        // condition: b-jet veto with the chosen tagger
	    if(nbjets[vetoTagger] > 0)                  continue;

        // fill new tree
	    ++debugcounter;
        outputTreePtr->Fill();
    }

    // write new tree
    std::cout << "finished event loop, now writing tree" << std::endl;
    outputTreePtr->Write( "",  BIT(2) );

    // close output file
    outputFilePtr->Close();

    // print some results
    std::cout << "skimmed "     << inputFilePath    << std::endl;
    std::cout << "processed "   << nEntries         << " entries, ";
    std::cout << "of which "    << debugcounter     << " passed the selections" << std::endl;
}


int main( int argc, char* argv[] ){
    std::cerr << "###starting###"                   << std::endl;
    if( argc != 4 && argc != 5 )
    {
        std::cerr               << "### ERROR ###: skimmer requires three or four arguments to run :"   << std::endl;
	    std::cerr               << "input_file_path, output_file_path, nentries [, veto_tagger]"        << std::endl;
	    std::cerr               << "veto_tagger: DeepJet (default), PNet, RobustParT or UParT"          << std::endl;
        return -1;
    }

    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );
    std::string&                                    input_file_path     = argvStr[1];
    std::string&                                    output_directory    = argvStr[2];
    long int                                        nEntries_raw        = std::stol(argvStr[3]);
    if(nEntries_raw < 0)                            nEntries_raw        = 0;
    long unsigned                                   nEntries            = (long unsigned) nEntries_raw;
    Tagger                                          vetoTagger          = (argc == 5) ? taggerFromName(argvStr[4]) : kDeepJet;
    skimFile( input_file_path,                      output_directory, nEntries, vetoTagger );
    std::cerr << "###done###"                       << std::endl;
    return 0;
}
