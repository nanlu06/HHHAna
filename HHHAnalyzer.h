 
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 27 15:49:34 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/GluGluHToMuMu_M-125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/AA4BF847-3785-E811-B130-001E67DFFF5F.root
//////////////////////////////////////////////////////////

#ifndef HHHAnalyzer_h
#define HHHAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <map>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "TRandom.h"
#include "MainEvent.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif
// Header file for the classes stored in the TTree if any.

class HHHAnalyzer : public MainEvent {
 public :
   HHHAnalyzer(const TString &inputFileList="foo.txt", const char *outFileName="histo.root", TString dataset="data",const char *isData="F", TString year_num="2017");
   virtual ~HHHAnalyzer();
   void Analyze(bool isData, int option, string outputFileName, string label);
   
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  float getPileupWeight(int);
  float getPileupWeightUp(int);
  float getPileupWeightDown(int);
  void     EventLoop(string, const char *, const char *);
  void     cal_sumOfgw(string , const char *);
  void     SkimTChain(const char *);
  //declare any specific function required
  
  void clearTreeVectors();
  void BookTreeBranches();
  bool DataIs;
  TString year;
  std::string yearst;
  std::map<std::string,float> btag_cut;
  std::map<std::string, float> luminosity;
  std::map<std::string, float> xs;
  std::map<std::string, std::map<std::string, float>> sumOfgenw;

  TH1D *h_sumOfgw = new TH1D("h_sumOfgenWeight","h_sumOfgenWeight",1,0,1);
  TH1D *h_sumOfgpw = new TH1D("h_sumOfgenpuWeight","h_sumOfgenpuWeight",1,0,1);
  
  TFile *pileupWeightFile;
  TH1F *pileupWeightHist, *pileupWeightSysUpHist, *pileupWeightSysDownHist;
    
  TFile *oFile;
  //TFile *ohistFile;
  TTree* tree;
  uint          t_run;
  uint          t_luminosityBlock;
  ulong       t_event;
  float       t_genWeight;
  float       t_puWeight;
  float       t_puWeightUp;
  float       t_puWeightDown;
  float       t_PrefireWeight;
  float       t_PrefireWeight_Up;
  float       t_PrefireWeight_Down;
  int         t_mu1;
  int         t_mu2;
  int         t_index_trigm_mu;
  std::vector<int>           *t_El_genPartIdx;
  std::vector<UChar_t>       *t_El_genPartFlav;
  std::vector <int>          *t_El_charge;
  std::vector<float>         *t_El_pt;
  std::vector<float>         *t_El_phi;
  std::vector<float>         *t_El_eta;   
  std::vector<float>         *t_El_mass;
  std::vector <int>          *t_El_cutBased;   
  std::vector <int>          *t_El_tightCharge;   
  std::vector<bool>          *t_El_cutBased_HEEP;   
  std::vector<bool>          *t_El_isPFcand;   
  std::vector<float>         *t_El_pfRelIso03_all;   
  std::vector<float>         *t_El_pfRelIso03_chg;      
  std::vector<float>         *t_El_miniPFRelIso_all;
  std::vector<float>         *t_El_miniPFRelIso_chg;
  std::vector<float>         *t_El_dxy;   
  std::vector<float>         *t_El_dxyErr;   
  std::vector<float>         *t_El_dz;   
  std::vector<float>         *t_El_dzErr; 
  std::vector<float>         *t_El_sip3d; 
  std::vector<float>         *t_Electron_mvaFall17Iso; 
  std::vector<bool>         *t_Electron_mvaFall17Iso_WP80;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17Iso_WP90;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17Iso_WPL;   //[nElectron]
  std::vector<float>        *t_Electron_mvaFall17noIso;
  std::vector<bool>         *t_Electron_mvaFall17noIso_WP80;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17noIso_WP90;   //[nElectron]
  std::vector<bool>         *t_Electron_mvaFall17noIso_WPL;   //[nElectron]
  
  std::vector<int>           *t_Mu_genPartIdx;
  std::vector<UChar_t>       *t_Mu_genPartFlav;
  std::vector<int>           *t_Mu_charge;  
  std::vector<float>         *t_Mu_EffSF_TRIG;
  std::vector<float>         *t_Mu_EffSFErr_TRIG;
  std::vector<float>         *t_Mu_EffSF_ID;
  std::vector<float>         *t_Mu_EffSF_ID_stat;
  std::vector<float>         *t_Mu_EffSF_ID_syst;
  std::vector<float>         *t_Mu_EffSF_ISO;
  std::vector<float>         *t_Mu_EffSF_ISO_stat;
  std::vector<float>         *t_Mu_EffSF_ISO_syst;
  std::vector<float>         *t_Mu_pt;   
  std::vector<float>         *t_Mu_ptErr;   
  std::vector<float>         *t_Mu_phi;   
  std::vector<float>         *t_Mu_eta;   
  std::vector<float>         *t_Mu_mass;  
  std::vector<float>         *t_Mu_dxy;   
  std::vector<float>         *t_Mu_dxyErr;   
  std::vector<float>         *t_Mu_dz;   
  std::vector<float>         *t_Mu_dzErr;
  std::vector<float>         *t_Mu_sip3d; 
  std::vector<float>         *t_Mu_pfRelIso03_all;   
  std::vector<float>         *t_Mu_pfRelIso03_chg;   
  std::vector<float>         *t_Mu_pfRelIso04_all;   
  std::vector<float>         *t_Mu_miniPFRelIso_all;
  std::vector<float>         *t_Mu_miniPFRelIso_chg;
  std::vector<int>           *t_Mu_tightCharge;   
  std::vector<bool>          *t_Mu_isPFcand;  
  std::vector<bool>          *t_Mu_istracker;
  std::vector<bool>          *t_Mu_isglobal; 
  std::vector<bool>          *t_Mu_mediumId;   
  std::vector<bool>          *t_Mu_softId;   
  std::vector<bool>          *t_Mu_tightId;    
  std::vector<int>           *t_Mu_nStations;   
  std::vector<int>           *t_Mu_nTrackerLayers;   

  float gen_higg0_pt;
  float gen_higg0_eta;
  float gen_higg0_phi;
  float gen_higg0_m;
    
  float gen_higg1_pt;
  float gen_higg1_eta;
  float gen_higg1_phi;
  float gen_higg1_m;
    
  float gen_higg2_pt;
  float gen_higg2_eta;
  float gen_higg2_phi;
  float gen_higg2_m;
    
  float gen_b0_pt;
  float gen_b0_eta;
  float gen_b0_phi;
  float gen_b0_m;

  float gen_b1_pt;
  float gen_b1_eta;
  float gen_b1_phi;
  float gen_b1_m;
    
  float gen_b2_pt;
  float gen_b2_eta;
  float gen_b2_phi;
  float gen_b2_m;

  float gen_b3_pt;
  float gen_b3_eta;
  float gen_b3_phi;
  float gen_b3_m;
    
  float gen_b4_pt;
  float gen_b4_eta;
  float gen_b4_phi;
  float gen_b4_m;

  float gen_b5_pt;
  float gen_b5_eta;
  float gen_b5_phi;
  float gen_b5_m;
    
  float genHHH_pt;
  float genHHH_eta;
  float genHHH_phi;
  float genHHH_m;
    
  float FatJet1_xbb;
  float FatJet1_pt;
  float FatJet1_eta;
  float FatJet1_phi;
  float FatJet1_msoftdrop;
  float FatJet1_particleNet_mass;
  
  float FatJet2_xbb;
  float FatJet2_pt;
  float FatJet2_eta;
  float FatJet2_phi;
  float FatJet2_msoftdrop;
  float FatJet2_particleNet_mass;
    
  float FatJet3_xbb;
  float FatJet3_pt;
  float FatJet3_eta;
  float FatJet3_phi;
  float FatJet3_msoftdrop;
  float FatJet3_particleNet_mass;
    
  std::vector<float>         *t_FatJet_area;  
  std::vector<float>         *t_FatJet_particleNet_HbbvsQCD;  
  std::vector<float>         *t_FatJet_particleNet_mass;  
  std::vector<float>         *t_FatJet_btagDeepB;  
  std::vector<float>         *t_FatJet_eta;  
  std::vector<float>         *t_FatJet_mass;  
  std::vector<float>         *t_FatJet_msoftdrop;  
  std::vector<float>         *t_FatJet_n2b1;  
  std::vector<float>         *t_FatJet_n3b1;  
  std::vector<float>         *t_FatJet_phi;  
  std::vector<float>         *t_FatJet_pt;  
  std::vector<float>         *t_FatJet_tau1;  
  std::vector<float>         *t_FatJet_tau2;  
  std::vector<float>         *t_FatJet_tau3;  
  std::vector<float>         *t_FatJet_tau4;  
  std::vector<int>           *t_FatJet_jetId;  
  std::vector<int>           *t_FatJet_subJetIdx1;  
  std::vector<int>           *t_FatJet_subJetIdx2;  

  std::vector<float>         *t_SubJet_btagCMVA;   
  std::vector<float>         *t_SubJet_btagCSVV2;   
  std::vector<float>         *t_SubJet_btagDeepB;   
  std::vector<float>         *t_SubJet_eta;   
  std::vector<float>         *t_SubJet_mass;   
  std::vector<float>         *t_SubJet_n2b1;   
  std::vector<float>         *t_SubJet_n3b1;   
  std::vector<float>         *t_SubJet_phi;   
  std::vector<float>         *t_SubJet_pt;   
  std::vector<float>         *t_SubJet_tau1;   
  std::vector<float>         *t_SubJet_tau2;   
  std::vector<float>         *t_SubJet_tau3;   
  std::vector<float>         *t_SubJet_tau4;   
  int t_nJet;
  std::vector<float>         *t_Jet_area;     
  std::vector<float>         *t_Jet_btagDeepFlavB;   
  std::vector<float>         *t_Jet_chEmEF;   
  std::vector<float>         *t_Jet_chHEF;   
  std::vector<float>         *t_Jet_eta;   
  std::vector<float>         *t_Jet_mass;   
  std::vector<float>         *t_Jet_neEmEF;   
  std::vector<float>         *t_Jet_neHEF;   
  std::vector<float>         *t_Jet_phi;   
  std::vector<float>         *t_Jet_pt;   
  std::vector<float>         *t_Jet_qgl;   
  std::vector<int>           *t_Jet_jetId;   
  std::vector<int>           *t_Jet_nConstituents;   
  std::vector<int>           *t_Jet_nElectrons;   
  std::vector<int>           *t_Jet_nMuons;   
  std::vector<int>           *t_Jet_puId;   

  int t_nbJet;
  std::vector<float>         *t_bJet_area;    
  std::vector<float>         *t_bJet_btagDeepFlavB;   
  std::vector<float>         *t_bJet_chEmEF;   
  std::vector<float>         *t_bJet_chHEF;   
  std::vector<float>         *t_bJet_eta;   
  std::vector<float>         *t_bJet_mass;   
  std::vector<float>         *t_bJet_neEmEF;   
  std::vector<float>         *t_bJet_neHEF;   
  std::vector<float>         *t_bJet_phi;   
  std::vector<float>         *t_bJet_pt;   
  std::vector<float>         *t_bJet_qgl;   
  std::vector<int>           *t_bJet_jetId;   
  std::vector<int>           *t_bJet_nConstituents;   
  std::vector<int>           *t_bJet_nElectrons;   
  std::vector<int>           *t_bJet_nMuons;   
  std::vector<int>           *t_bJet_puId;   
  std::vector<double>         *t_bJet_SF;
  std::vector<double>         *t_bJet_SFup;
  std::vector<double>         *t_bJet_SFdown;

  float      t_PV_ndof;
  float      t_PV_x;
  float      t_PV_y;
  float      t_PV_z;
  int        t_PV_npvs;
  int        t_PV_npvsGood;

  UInt_t     t_nLHEPdfWeight;
  UInt_t     t_nLHEScaleWeight;
  UInt_t     t_nPSWeight;
    
  float genweight;
  int cat;
  float HHH_m;
  float HHH_pt;
  float HHH_eta;
  float HHH_phi;
  float HH12_m;
  float HH12_pt;
  float HH12_eta;
  float HH12_phi;
  float HH23_m;
  float HH23_pt;
  float HH23_eta;
  float HH23_phi;
  float HH13_m;
  float HH13_pt;
  float HH13_eta;
  float HH13_phi;
  float H1_m;
  float H1_pt;
  float H1_eta;
  float H1_phi;
  float H2_m;
  float H2_pt;
  float H2_eta;
  float H2_phi;
  float H3_m;
  float H3_pt;
  float H3_eta;
  float H3_phi;  
  float bJet1_btagDeepFlavB;
  float bJet2_btagDeepFlavB;
  float bJet3_btagDeepFlavB;
  float bJet4_btagDeepFlavB;
  float bJet5_btagDeepFlavB;
  float bJet6_btagDeepFlavB;
  float bJet1_pt;
  float bJet2_pt;
  float bJet3_pt;
  float bJet4_pt;
  float bJet5_pt;
  float bJet6_pt;
  float bJet1_eta;
  float bJet2_eta;
  float bJet3_eta;
  float bJet4_eta;
  float bJet5_eta;
  float bJet6_eta;
    
};

#endif

#ifdef HHHAnalyzer_cxx
HHHAnalyzer::HHHAnalyzer(const TString &inputFileList, const char *outFileName, TString dataset, const char *isData, TString year_num) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if(*isData!='T') DataIs = false;
  else DataIs = true;
  year = year_num;
  yearst = std::string(year.Data());

  h_sumOfgw->SetBinContent(1,0.0);
  h_sumOfgpw->SetBinContent(1,0.0);
    
  //b-tag score selection loose
  btag_cut["2016"] = 0.0614; //0.3093 //0.7221
  btag_cut["2017"] = 0.0521; //0.3033 //0.7489
  btag_cut["2018"] = 0.0494; //0.2770 //0.7264
    
  //luminosity for each year
  //https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
  luminosity["2016"] = 35.9;
  luminosity["2017"] = 41.5;
  luminosity["2018"] = 59.8;
    
  //xs in unit of fb
  //e.g. 
  //float xsHH = 31.05; //fb
  float xsHHH = 0.1;
  float BRHbb = 5.824E-01;
  float BRHgg = 2.270E-03;
  float BRHtautau = 6.272E-02;
  //ref https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json   
  //ref https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGHH
  xs["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 31.05*5.824E-01*2.270E-03*2;
  xs["hhh6b"] = xsHHH*BRHbb*BRHbb*BRHbb;
  xs["hhh6b_test"] = xsHHH*BRHbb*BRHbb*BRHbb;
    
  std::map<std::string, float> sumOfgenw_2016, sumOfgenw_2017, sumOfgenw_2018;

  sumOfgenw_2018["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 5402.244961268351-122.629;
  sumOfgenw_2017["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 5402.244961268351-122.629;
  sumOfgenw_2016["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 5402.244961268351-122.629;
  sumOfgenw_2017["hhh6b"] = 236500;  
  sumOfgenw_2017["hhh6b_test"] = 500;
    
  sumOfgenw["2016"] = sumOfgenw_2016;
  sumOfgenw["2017"] = sumOfgenw_2017;
  sumOfgenw["2018"] = sumOfgenw_2018;

  TChain *tree = new TChain("Events");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  MainEvent::Init(tree);
  oFile = new TFile(outFileName, "recreate");
  TString histname(outFileName);
  //ohistFile = new TFile("hist_"+histname, "recreate");
  BookTreeBranches();
}

float HHHAnalyzer::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

float HHHAnalyzer::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

float HHHAnalyzer::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}
    
bool HHHAnalyzer::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;

  while ( getline (infile,buffer) )
  {
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
    /*
    if(!DataIs){
      std::string friend_buffer = "";

      std::string delimiter = "/";

      size_t pos = 0;
      std::string buffer_f = "/storage/user/nlu/Hmm/puJEC/"+yearst+"/17July_v1/"; //"/mnt/hadoop/store/user/nlu/Hmm/ntuple/2016/Nano14Dec2018/MC/puJEC/16June_v1/";
      std::string buffer_tmp = buffer;
      int dataname_pos = 7;
      if(buffer.find("storage")!=std::string::npos) dataname_pos = 8;
      int cout=0;
      while ((pos = buffer_tmp.find(delimiter)) != std::string::npos) {
         cout++;
         friend_buffer = buffer_tmp.substr(0, pos);
         std::cout <<"friend_buffer: "<<friend_buffer << std::endl;
         if(cout==dataname_pos){ //8 for storage
                buffer_f +=friend_buffer;
         }
         buffer_tmp.erase(0, pos + delimiter.length());
      }
      std::cout << friend_buffer << std::endl;
      std::cout << "buffer_tmp: "<<buffer_tmp << std::endl;

      buffer_f += buffer_tmp;
      std::cout << "Adding friend tree from " << buffer_f.c_str() << std::endl;
      chain->AddFriend("Friends",buffer_f.c_str());
    }
    */
  }
  infile.close();

  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

 HHHAnalyzer::~HHHAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   h_sumOfgw->Write();
   h_sumOfgpw->Write();
   oFile->Write();
   oFile->Close();
}

Long64_t HHHAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HHHAnalyzer::clearTreeVectors(){
  t_run=0;
  t_luminosityBlock=0;
  t_event=0;
  t_genWeight = -999.;
  t_puWeight = -999.;
  t_puWeightUp = -999.;
  t_puWeightDown = -999.;
  t_PrefireWeight = 1.0;
  t_PrefireWeight_Up = 1.0;
  t_PrefireWeight_Down = 1.0;
  t_mu1=-999; 
  t_mu2=-999;
  t_index_trigm_mu=-999;
  t_nJet=0;
  t_nbJet=0;
  t_El_genPartIdx->clear();
  t_El_genPartFlav->clear();
  t_El_charge->clear();
  t_El_pt->clear();
  t_El_phi->clear();
  t_El_eta->clear();   
  t_El_mass->clear();
  t_El_cutBased->clear();   
  t_El_tightCharge->clear();   
  t_El_cutBased_HEEP->clear();   
  t_El_isPFcand->clear();   
  t_El_pfRelIso03_all->clear();   
  t_El_pfRelIso03_chg->clear();      
  t_El_miniPFRelIso_all->clear();
  t_El_miniPFRelIso_chg->clear();
  t_El_dxy->clear();   
  t_El_dxyErr->clear();   
  t_El_dz->clear();   
  t_El_dzErr->clear();  
  t_El_sip3d->clear();
  t_Electron_mvaFall17Iso->clear(); 
  t_Electron_mvaFall17Iso_WP80->clear();
  t_Electron_mvaFall17Iso_WP90->clear();
  t_Electron_mvaFall17Iso_WPL->clear();
  t_Electron_mvaFall17noIso->clear();
  t_Electron_mvaFall17noIso_WP80->clear();
  t_Electron_mvaFall17noIso_WP90->clear();
  t_Electron_mvaFall17noIso_WPL->clear();

  t_Mu_genPartIdx->clear();
  t_Mu_genPartFlav->clear();  
  t_Mu_charge->clear();
  t_Mu_EffSF_TRIG->clear();
  t_Mu_EffSFErr_TRIG->clear();
  t_Mu_EffSF_ID->clear();
  t_Mu_EffSF_ID_stat->clear();
  t_Mu_EffSF_ID_syst->clear();
  t_Mu_EffSF_ISO->clear();
  t_Mu_EffSF_ISO_stat->clear();
  t_Mu_EffSF_ISO_syst->clear();
  t_Mu_pt->clear();   
  t_Mu_ptErr->clear();   
  t_Mu_phi->clear();   
  t_Mu_eta->clear();   
  t_Mu_mass->clear();  
  t_Mu_dxy->clear();   
  t_Mu_dxyErr->clear();   
  t_Mu_dz->clear();   
  t_Mu_dzErr->clear(); 
  t_Mu_sip3d->clear();
  t_Mu_pfRelIso03_all->clear();   
  t_Mu_pfRelIso03_chg->clear();   
  t_Mu_pfRelIso04_all->clear();   
  t_Mu_miniPFRelIso_all->clear();
  t_Mu_miniPFRelIso_chg->clear();
  t_Mu_tightCharge->clear();   
  t_Mu_isPFcand->clear();   
  t_Mu_isglobal->clear();
  t_Mu_istracker->clear();
  t_Mu_mediumId->clear();   
  t_Mu_softId->clear();   
  t_Mu_tightId->clear();    
  t_Mu_nStations->clear();   
  t_Mu_nTrackerLayers->clear();   

  gen_higg0_pt=-1000;
  gen_higg0_eta=-1000;
  gen_higg0_phi=-1000;
  gen_higg0_m=-1000;
    
  gen_higg1_pt=-1000;
  gen_higg1_eta=-1000;
  gen_higg1_phi=-1000;
  gen_higg1_m=-1000;

  gen_higg2_pt=-1000;
  gen_higg2_eta=-1000;
  gen_higg2_phi=-1000;
  gen_higg2_m=-1000;

  gen_b0_pt=-1000;
  gen_b0_eta=-1000;
  gen_b0_phi=-1000;
  gen_b0_m=-1000;
    
  gen_b1_pt=-1000;
  gen_b1_eta=-1000;
  gen_b1_phi=-1000;
  gen_b1_m=-1000;
    
  gen_b2_pt=-1000;
  gen_b2_eta=-1000;
  gen_b2_phi=-1000;
  gen_b2_m=-1000;
    
  gen_b3_pt=-1000;
  gen_b3_eta=-1000;
  gen_b3_phi=-1000;
  gen_b3_m=-1000;
    
  gen_b4_pt=-1000;
  gen_b4_eta=-1000;
  gen_b4_phi=-1000;
  gen_b4_m=-1000;
    
  gen_b5_pt=-1000;
  gen_b5_eta=-1000;
  gen_b5_phi=-1000;
  gen_b5_m=-1000;
    
  genHHH_pt=-1000;
  genHHH_eta=-1000;
  genHHH_phi=-1000;
  genHHH_m=-1000;
    
  FatJet1_xbb=-1000;
  FatJet1_pt=-1000;
  FatJet1_eta=-1000;
  FatJet1_phi=-1000;
  FatJet1_msoftdrop=-1000;
  FatJet1_particleNet_mass=-1000;
    
  FatJet2_xbb=-1000;
  FatJet2_pt=-1000;
  FatJet2_eta=-1000;
  FatJet2_phi=-1000;
  FatJet2_msoftdrop=-1000;
  FatJet2_particleNet_mass=-1000;
    
  FatJet3_xbb=-1000;
  FatJet3_pt=-1000;
  FatJet3_eta=-1000;
  FatJet3_phi=-1000;
  FatJet3_msoftdrop=-1000;
  FatJet3_particleNet_mass=-1000;
    
  t_FatJet_area->clear();  
  t_FatJet_particleNet_HbbvsQCD->clear();  
  t_FatJet_particleNet_mass->clear();  
  t_FatJet_btagDeepB->clear();  
  t_FatJet_eta->clear();  
  t_FatJet_mass->clear();  
  t_FatJet_msoftdrop->clear();  
  t_FatJet_n2b1->clear();  
  t_FatJet_n3b1->clear();  
  t_FatJet_phi->clear();  
  t_FatJet_pt->clear();  
  t_FatJet_tau1->clear();  
  t_FatJet_tau2->clear();  
  t_FatJet_tau3->clear();  
  t_FatJet_tau4->clear();  
  t_FatJet_jetId->clear();  
  t_FatJet_subJetIdx1->clear();  
  t_FatJet_subJetIdx2->clear();  

  t_SubJet_btagCMVA->clear();   
  t_SubJet_btagCSVV2->clear();   
  t_SubJet_btagDeepB->clear();   
  t_SubJet_eta->clear();   
  t_SubJet_mass->clear();   
  t_SubJet_n2b1->clear();   
  t_SubJet_n3b1->clear();   
  t_SubJet_phi->clear();   
  t_SubJet_pt->clear();   
  t_SubJet_tau1->clear();   
  t_SubJet_tau2->clear();   
  t_SubJet_tau3->clear();   
  t_SubJet_tau4->clear();   
  
  t_Jet_area->clear();     
  t_Jet_btagDeepFlavB->clear();   
  t_Jet_chEmEF->clear();   
  t_Jet_chHEF->clear();   
  t_Jet_eta->clear();   
  t_Jet_mass->clear();   
  t_Jet_neEmEF->clear();   
  t_Jet_neHEF->clear();   
  t_Jet_phi->clear();   
  t_Jet_pt->clear();   
  t_Jet_qgl->clear();   
  t_Jet_jetId->clear();   
  t_Jet_nConstituents->clear();   
  t_Jet_nElectrons->clear();   
  t_Jet_nMuons->clear();   
  t_Jet_puId->clear();   

  t_bJet_area->clear();     
  t_bJet_btagDeepFlavB->clear();   
  t_bJet_chEmEF->clear();   
  t_bJet_chHEF->clear();   
  t_bJet_eta->clear();   
  t_bJet_mass->clear();   
  t_bJet_neEmEF->clear();   
  t_bJet_neHEF->clear();   
  t_bJet_phi->clear();   
  t_bJet_pt->clear();   
  t_bJet_qgl->clear();   
  t_bJet_jetId->clear();   
  t_bJet_nConstituents->clear();   
  t_bJet_nElectrons->clear();   
  t_bJet_nMuons->clear();   
  t_bJet_puId->clear();   
  t_bJet_SF->clear();
  t_bJet_SFup->clear();
  t_bJet_SFdown->clear();

  //SoftActivityJetHT5
  //SoftActivityJetNjets5
  
  t_PV_ndof=-1000;
  t_PV_x=-1000;
  t_PV_y=-1000;
  t_PV_z=-1000;
  t_PV_npvs=-1000;  
  t_PV_npvsGood=-1000;
    
  genweight=0;
  cat = -1;
  HHH_m = -999.;
  HHH_pt = -999.;
  HHH_eta = -999.;
  HHH_phi = -999.;
  HH12_m = -999.;
  HH12_pt = -999.;
  HH12_eta = -999.;
  HH12_phi = -999.;
  HH23_m = -999.;
  HH23_pt = -999.;
  HH23_eta = -999.;
  HH23_phi = -999.;
  HH13_m = -999.;
  HH13_pt = -999.;
  HH13_eta = -999.;
  HH13_phi = -999.;
  H1_m = -999.;
  H1_pt = -999.;
  H1_eta = -999.;
  H1_phi = -999.;
  H2_m = -999.;
  H2_pt = -999.;
  H2_eta = -999.;
  H2_phi = -999.;
  H3_m = -999.;
  H3_pt = -999.;
  H3_eta = -999.;
  H3_phi = -999.;
  bJet1_btagDeepFlavB = -999.;
  bJet2_btagDeepFlavB = -999.;
  bJet3_btagDeepFlavB = -999.;
  bJet4_btagDeepFlavB = -999.;
  bJet5_btagDeepFlavB = -999.;
  bJet6_btagDeepFlavB = -999.;
  bJet1_pt = -999.;
  bJet2_pt = -999.;
  bJet3_pt = -999.;
  bJet4_pt = -999.;
  bJet5_pt = -999.;
  bJet6_pt = -999.;
  bJet1_eta = -999.;
  bJet2_eta = -999.;
  bJet3_eta = -999.;
  bJet4_eta = -999.;
  bJet5_eta = -999.;
  bJet6_eta = -999.;
}

void HHHAnalyzer::BookTreeBranches(){
  tree = new TTree("tree","tree");
  tree->SetAutoSave(10000);

  tree->Branch("t_run", &t_run,"t_run/i");
  tree->Branch("t_luminosityBlock", &t_luminosityBlock,"t_luminosityBlock/i");
  tree->Branch("t_event", &t_event,"t_event/l");
  tree->Branch("t_genWeight", &t_genWeight,"t_genWeight/F");
  tree->Branch("t_puWeight", &t_puWeight,"t_puWeight/F");
  tree->Branch("t_puWeightUp", &t_puWeightUp,"t_puWeightUp/F");
  tree->Branch("t_puWeightDown", &t_puWeightDown,"t_puWeightDown/F");
  //tree->Branch("t_pileupWeight", &t_pileupWeight,"t_pileupWeight/F");
  //tree->Branch("t_pileupupWeight", &t_pileupupWeight,"t_pileupupWeight/F");
  //tree->Branch("t_pileupdnWeight", &t_pileupdnWeight,"t_pileupdnWeight/F");
  tree->Branch("t_PrefireWeight", &t_PrefireWeight, "t_PrefireWeight/F");
  tree->Branch("t_PrefireWeight_Up", &t_PrefireWeight_Up, "t_PrefireWeight_Up/F");
  tree->Branch("t_PrefireWeight_Down", &t_PrefireWeight_Down, "t_PrefireWeight_Down/F");
  tree->Branch("t_mu1", &t_mu1,"t_mu1/I");
  tree->Branch("t_mu2", &t_mu2,"t_mu2/I");
  tree->Branch("t_index_trigm_mu", &t_index_trigm_mu, "t_index_trigm_mu/I");

  t_El_genPartIdx= new std::vector<int>();
  t_El_genPartFlav= new std::vector<UChar_t>();
  t_El_charge= new std::vector<int>();
  t_El_pt= new std::vector<float>();
  t_El_phi= new std::vector<float>();
  t_El_eta= new std::vector<float>();   
  t_El_mass= new std::vector<float>();
  t_El_cutBased= new std::vector<int>();   
  t_El_tightCharge= new std::vector<int>();   
  t_El_cutBased_HEEP= new std::vector<bool>();   
  t_El_isPFcand= new std::vector<bool>();   
  t_El_pfRelIso03_all= new std::vector<float>();   
  t_El_pfRelIso03_chg= new std::vector<float>();      
  t_El_miniPFRelIso_all= new std::vector<float>();
  t_El_miniPFRelIso_chg= new std::vector<float>();
  t_El_dxy= new std::vector<float>();   
  t_El_dxyErr= new std::vector<float>();   
  t_El_dz= new std::vector<float>();   
  t_El_dzErr= new std::vector<float>();  
  t_El_sip3d= new std::vector<float>();
  t_Electron_mvaFall17Iso= new std::vector<float>(); 
  t_Electron_mvaFall17Iso_WP80= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17Iso_WP90= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17Iso_WPL= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso= new std::vector<float>();
  t_Electron_mvaFall17noIso_WP80= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso_WP90= new std::vector<bool>();   //[nElectron]
  t_Electron_mvaFall17noIso_WPL= new std::vector<bool>();   //[nElectron]
 
  tree->Branch("t_El_genPartIdx",   "vector<int>",        &t_El_genPartIdx);
  tree->Branch("t_El_genPartFlav",   "vector<UChar_t>",       &t_El_genPartFlav); 
  tree->Branch("t_El_charge",        "vector<int>"         ,&t_El_charge);
  tree->Branch("t_El_pt",        "vector<float>"         ,&t_El_pt);
  tree->Branch("t_El_phi",        "vector<float>"         ,&t_El_phi);
  tree->Branch("t_El_eta",        "vector<float>"         ,&t_El_eta);   
  tree->Branch("t_El_mass",        "vector<float>"         ,&t_El_mass);
  tree->Branch("t_El_cutBased",        "vector<int>"         ,&t_El_cutBased);   
  tree->Branch("t_El_tightCharge",        "vector<int>"         ,&t_El_tightCharge);   
  tree->Branch("t_El_cutBased_HEEP",        "vector<bool>"         ,&t_El_cutBased_HEEP);   
  tree->Branch("t_El_isPFcand",        "vector<bool>"         ,&t_El_isPFcand);   
  tree->Branch("t_El_pfRelIso03_all",        "vector<float>"         ,&t_El_pfRelIso03_all);   
  tree->Branch("t_El_pfRelIso03_chg",        "vector<float>"         ,&t_El_pfRelIso03_chg);      
  tree->Branch("t_El_miniPFRelIso03_all",        "vector<float>"         ,&t_El_miniPFRelIso_all);
  tree->Branch("t_El_miniPFRelIso03_chg",        "vector<float>"         ,&t_El_miniPFRelIso_chg);
  tree->Branch("t_El_dxy",        "vector<float>"         ,&t_El_dxy);   
  tree->Branch("t_El_dxyErr",        "vector<float>"         ,&t_El_dxyErr);   
  tree->Branch("t_El_dz",        "vector<float>"         ,&t_El_dz);   
  tree->Branch("t_El_dzErr",        "vector<float>"         ,&t_El_dzErr);   
  tree->Branch("t_El_sip3d",        "vector<float>"         ,&t_El_sip3d);
  tree->Branch("t_Electron_mvaFall17Iso",        "vector<float>",         &t_Electron_mvaFall17Iso); 
  tree->Branch("t_Electron_mvaFall17Iso_WP80",        "vector<bool>",         &t_Electron_mvaFall17Iso_WP80);
  tree->Branch("t_Electron_mvaFall17Iso_WP90",        "vector<bool>",         &t_Electron_mvaFall17Iso_WP90);
  tree->Branch("t_Electron_mvaFall17Iso_WPL",        "vector<bool>",         &t_Electron_mvaFall17Iso_WPL);
  tree->Branch("t_Electron_mvaFall17noIso",        "vector<float>",         &t_Electron_mvaFall17noIso);
  tree->Branch("t_Electron_mvaFall17noIso_WP80",        "vector<bool>",         &t_Electron_mvaFall17noIso_WP80);
  tree->Branch("t_Electron_mvaFall17noIso_WP90",        "vector<bool>",         &t_Electron_mvaFall17noIso_WP90);
  tree->Branch("t_Electron_mvaFall17noIso_WPL",        "vector<bool>",         &t_Electron_mvaFall17noIso_WPL);
 
  t_Mu_genPartIdx= new std::vector<int>();
  t_Mu_genPartFlav= new std::vector<UChar_t>(); 
  t_Mu_charge= new std::vector<int>();   
  t_Mu_EffSF_TRIG= new std::vector<float>();
  t_Mu_EffSFErr_TRIG= new std::vector<float>();  
  t_Mu_EffSF_ID= new std::vector<float>();
  t_Mu_EffSF_ID_stat= new std::vector<float>();
  t_Mu_EffSF_ID_syst= new std::vector<float>();
  t_Mu_EffSF_ISO= new std::vector<float>();
  t_Mu_EffSF_ISO_stat= new std::vector<float>();
  t_Mu_EffSF_ISO_syst= new std::vector<float>();
  t_Mu_pt= new std::vector<float>();   
  t_Mu_ptErr= new std::vector<float>();   
  t_Mu_phi= new std::vector<float>();   
  t_Mu_eta= new std::vector<float>();   
  t_Mu_mass= new std::vector<float>();  
  t_Mu_dxy= new std::vector<float>();   
  t_Mu_dxyErr= new std::vector<float>();   
  t_Mu_dz= new std::vector<float>();   
  t_Mu_dzErr= new std::vector<float>(); 
  t_Mu_sip3d= new std::vector<float>();
  t_Mu_pfRelIso03_all= new std::vector<float>();   
  t_Mu_pfRelIso03_chg= new std::vector<float>();   
  t_Mu_pfRelIso04_all= new std::vector<float>();   
  t_Mu_miniPFRelIso_all= new std::vector<float>();
  t_Mu_miniPFRelIso_chg= new std::vector<float>();
  t_Mu_tightCharge= new std::vector<int>();   
  t_Mu_isPFcand= new std::vector<bool>();  
  t_Mu_isglobal= new std::vector<bool>();
  t_Mu_istracker= new std::vector<bool>(); 
  t_Mu_mediumId= new std::vector<bool>();   
  t_Mu_softId= new std::vector<bool>();   
  t_Mu_tightId= new std::vector<bool>();    
  t_Mu_nStations= new std::vector<int>();   
  t_Mu_nTrackerLayers= new std::vector<int>();   

  tree->Branch("t_Mu_genPartIdx"    , "vector<int>"     ,&t_Mu_genPartIdx);
  tree->Branch("t_Mu_genPartFlav",  "vector<UChar_t>",   &t_Mu_genPartFlav);
  tree->Branch("t_Mu_charge"    , "vector<int>"         ,&t_Mu_charge );
  tree->Branch("t_Mu_EffSF_TRIG",    "vector<float>"   ,&t_Mu_EffSF_TRIG);
  tree->Branch("t_Mu_EffSFErr_TRIG",    "vector<float>"   ,&t_Mu_EffSFErr_TRIG);
  tree->Branch("t_Mu_EffSF_ID",    "vector<float>"   ,&t_Mu_EffSF_ID);
  tree->Branch("t_Mu_EffSF_ID_stat",    "vector<float>"   ,&t_Mu_EffSF_ID_stat);
  tree->Branch("t_Mu_EffSF_ID_syst",    "vector<float>"   ,&t_Mu_EffSF_ID_syst);
  tree->Branch("t_Mu_EffSF_ISO",    "vector<float>"   ,&t_Mu_EffSF_ISO);
  tree->Branch("t_Mu_EffSF_ISO_stat",    "vector<float>"   ,&t_Mu_EffSF_ISO_stat);
  tree->Branch("t_Mu_EffSF_ISO_syst",    "vector<float>"   ,&t_Mu_EffSF_ISO_syst);
  tree->Branch("t_Mu_pt"    , "vector<float>"         ,&t_Mu_pt );   
  tree->Branch("t_Mu_ptErr"    , "vector<float>"         ,&t_Mu_ptErr );   
  tree->Branch("t_Mu_phi"    , "vector<float>"         ,&t_Mu_phi );   
  tree->Branch("t_Mu_eta"    , "vector<float>"         ,&t_Mu_eta );   
  tree->Branch("t_Mu_mass"    , "vector<float>"         ,&t_Mu_mass );  
  tree->Branch("t_Mu_dxy"    , "vector<float>"         ,&t_Mu_dxy );   
  tree->Branch("t_Mu_dxyErr"    , "vector<float>"         ,&t_Mu_dxyErr );   
  tree->Branch("t_Mu_dz"    , "vector<float>"         ,&t_Mu_dz );   
  tree->Branch("t_Mu_dzErr"    , "vector<float>"         ,&t_Mu_dzErr );   
  tree->Branch("t_Mu_sip3d"   , "vector<float>"         ,&t_Mu_sip3d );
  tree->Branch("t_Mu_pfRelIso03_all"    , "vector<float>"         ,&t_Mu_pfRelIso03_all );   
  tree->Branch("t_Mu_pfRelIso03_chg"    , "vector<float>"         ,&t_Mu_pfRelIso03_chg );   
  tree->Branch("t_Mu_pfRelIso04_all"    , "vector<float>"         ,&t_Mu_pfRelIso04_all );   
  tree->Branch("t_Mu_miniPFRelIso_all"    , "vector<float>"         ,&t_Mu_miniPFRelIso_all );
  tree->Branch("t_Mu_miniPFRelIso_chg"    , "vector<float>"         ,&t_Mu_miniPFRelIso_chg );
  tree->Branch("t_Mu_tightCharge"    , "vector<int>"         ,&t_Mu_tightCharge );   
  tree->Branch("t_Mu_isPFcand"    , "vector<bool>"         ,&t_Mu_isPFcand );   
  tree->Branch("t_Mu_isglobal"    , "vector<bool>"         ,&t_Mu_isglobal );
  tree->Branch("t_Mu_istracker"    , "vector<bool>"         ,&t_Mu_istracker );
  tree->Branch("t_Mu_mediumId"    , "vector<bool>"         ,&t_Mu_mediumId );   
  tree->Branch("t_Mu_softId"    , "vector<bool>"         ,&t_Mu_softId );   
  tree->Branch("t_Mu_tightId"    , "vector<bool>"         ,&t_Mu_tightId );    
  tree->Branch("t_Mu_nStations"    , "vector<int>"         ,&t_Mu_nStations );   
  tree->Branch("t_Mu_nTrackerLayers"    , "vector<int>"         ,&t_Mu_nTrackerLayers );
  /*
  gen_higg0_pt=-1000;
  gen_higg0_eta=-1000;
  gen_higg0_phi=-1000;
  gen_higg0_m=-1000;
    
  gen_higg1_pt=-1000;
  gen_higg1_eta=-1000;
  gen_higg1_phi=-1000;
  gen_higg1_m=-1000;
    
  gen_higg2_pt=-1000;
  gen_higg2_eta=-1000;
  gen_higg2_phi=-1000;
  gen_higg2_m=-1000;
    
  gen_b0_pt=-1000;
  gen_b0_eta=-1000;
  gen_b0_phi=-1000;
  gen_b0_m=-1000;
    
  gen_b1_pt=-1000;
  gen_b1_eta=-1000;
  gen_b1_phi=-1000;
  gen_b1_m=-1000;
  
  gen_b2_pt=-1000;
  gen_b2_eta=-1000;
  gen_b2_phi=-1000;
  gen_b2_m=-1000;
    
  gen_b3_pt=-1000;
  gen_b3_eta=-1000;
  gen_b3_phi=-1000;
  gen_b3_m=-1000;
    
  gen_b4_pt=-1000;
  gen_b4_eta=-1000;
  gen_b4_phi=-1000;
  gen_b4_m=-1000;
    
  gen_b5_pt=-1000;
  gen_b5_eta=-1000;
  gen_b5_phi=-1000;
  gen_b5_m=-1000;
    
  genHHH_pt=-1000;
  genHHH_eta=-1000;
  genHHH_phi=-1000;
  genHHH_m=-1000;
  */
  tree->Branch("gen_higg0_pt",   &gen_higg0_pt,"gen_higg0_pt/F");  
  tree->Branch("gen_higg0_eta",   &gen_higg0_eta,"gen_higg0_eta/F");
  tree->Branch("gen_higg0_phi",  &gen_higg0_phi,"gen_higg0_phi/F");
  tree->Branch("gen_higg0_m",   &gen_higg0_m,"gen_higg0_m/F");
    
  tree->Branch("gen_higg1_pt",   &gen_higg1_pt,"gen_higg1_pt/F");  
  tree->Branch("gen_higg1_eta",   &gen_higg1_eta,"gen_higg1_eta/F");
  tree->Branch("gen_higg1_phi",  &gen_higg1_phi,"gen_higg1_phi/F");
  tree->Branch("gen_higg1_m",   &gen_higg1_m,"gen_higg1_m/F");
    
  tree->Branch("gen_higg2_pt",   &gen_higg2_pt,"gen_higg2_pt/F");  
  tree->Branch("gen_higg2_eta",   &gen_higg2_eta,"gen_higg2_eta/F");
  tree->Branch("gen_higg2_phi",  &gen_higg2_phi,"gen_higg2_phi/F");
  tree->Branch("gen_higg2_m",   &gen_higg2_m,"gen_higg2_m/F");
    
  tree->Branch("gen_b0_pt",   &gen_b0_pt,"gen_b0_pt/F");  
  tree->Branch("gen_b0_eta",   &gen_b0_eta,"gen_b0_eta/F");
  tree->Branch("gen_b0_phi",  &gen_b0_phi,"gen_b0_phi/F");
  tree->Branch("gen_b0_m",   &gen_b0_m,"gen_b0_m/F");
    
  tree->Branch("gen_b1_pt",   &gen_b1_pt,"gen_b1_pt/F");  
  tree->Branch("gen_b1_eta",   &gen_b1_eta,"gen_b1_eta/F");
  tree->Branch("gen_b1_phi",  &gen_b1_phi,"gen_b1_phi/F");
  tree->Branch("gen_b1_m",   &gen_b1_m,"gen_b1_m/F");
    
  tree->Branch("gen_b2_pt",   &gen_b2_pt,"gen_b2_pt/F");  
  tree->Branch("gen_b2_eta",   &gen_b2_eta,"gen_b2_eta/F");
  tree->Branch("gen_b2_phi",  &gen_b2_phi,"gen_b2_phi/F");
  tree->Branch("gen_b2_m",   &gen_b2_m,"gen_b2_m/F");
   
  tree->Branch("gen_b3_pt",   &gen_b3_pt,"gen_b3_pt/F");  
  tree->Branch("gen_b3_eta",   &gen_b3_eta,"gen_b3_eta/F");
  tree->Branch("gen_b3_phi",  &gen_b3_phi,"gen_b3_phi/F");
  tree->Branch("gen_b3_m",   &gen_b3_m,"gen_b3_m/F");
    
  tree->Branch("gen_b4_pt",   &gen_b4_pt,"gen_b4_pt/F");  
  tree->Branch("gen_b4_eta",   &gen_b4_eta,"gen_b4_eta/F");
  tree->Branch("gen_b4_phi",  &gen_b4_phi,"gen_b4_phi/F");
  tree->Branch("gen_b4_m",   &gen_b4_m,"gen_b4_m/F");
    
  tree->Branch("gen_b5_pt",   &gen_b5_pt,"gen_b5_pt/F");  
  tree->Branch("gen_b5_eta",   &gen_b5_eta,"gen_b5_eta/F");
  tree->Branch("gen_b5_phi",  &gen_b5_phi,"gen_b5_phi/F");
  tree->Branch("gen_b5_m",   &gen_b5_m,"gen_b5_m/F");
   
  tree->Branch("genHHH_pt",   &genHHH_pt,"genHHH_pt/F");  
  tree->Branch("genHHH_eta",   &genHHH_eta,"genHHH_eta/F");
  tree->Branch("genHHH_phi",  &genHHH_phi,"genHHH_phi/F");
  tree->Branch("genHHH_m",   &genHHH_m,"genHHH_m/F");
    
  tree->Branch("FatJet1_xbb",   &FatJet1_xbb, "FatJet1_xbb/F");
  tree->Branch("FatJet1_pt",   &FatJet1_pt, "FatJet1_pt/F");
  tree->Branch("FatJet1_eta",   &FatJet1_eta, "FatJet1_eta/F");
  tree->Branch("FatJet1_phi",   &FatJet1_phi, "FatJet1_phi/F");
  tree->Branch("FatJet1_msoftdrop",   &FatJet1_msoftdrop, "FatJet1_msoftdrop/F");
  tree->Branch("FatJet1_particleNet_mass",   &FatJet1_particleNet_mass, "FatJet1_particleNet_mass/F");
  
  tree->Branch("FatJet2_xbb",   &FatJet2_xbb, "FatJet2_xbb/F");
  tree->Branch("FatJet2_pt",   &FatJet2_pt, "FatJet2_pt/F");
  tree->Branch("FatJet2_eta",   &FatJet2_eta, "FatJet2_eta/F");
  tree->Branch("FatJet2_phi",   &FatJet2_phi, "FatJet2_phi/F");
  tree->Branch("FatJet2_msoftdrop",   &FatJet2_msoftdrop, "FatJet2_msoftdrop/F");
  tree->Branch("FatJet2_particleNet_mass",   &FatJet2_particleNet_mass, "FatJet2_particleNet_mass/F");
  
  tree->Branch("FatJet3_xbb",   &FatJet3_xbb, "FatJet3_xbb/F");
  tree->Branch("FatJet3_pt",   &FatJet3_pt, "FatJet3_pt/F");
  tree->Branch("FatJet3_eta",   &FatJet3_eta, "FatJet3_eta/F");
  tree->Branch("FatJet3_phi",   &FatJet3_phi, "FatJet3_phi/F");
  tree->Branch("FatJet3_msoftdrop",   &FatJet3_msoftdrop, "FatJet3_msoftdrop/F");
  tree->Branch("FatJet3_particleNet_mass",   &FatJet3_particleNet_mass, "FatJet3_particleNet_mass/F");

  t_FatJet_area= new std::vector<float>();  
  t_FatJet_particleNet_HbbvsQCD= new std::vector<float>();  
  t_FatJet_particleNet_mass= new std::vector<float>();  
  t_FatJet_btagDeepB= new std::vector<float>();  
  t_FatJet_eta= new std::vector<float>();  
  t_FatJet_mass= new std::vector<float>();  
  t_FatJet_msoftdrop= new std::vector<float>();  
  t_FatJet_n2b1= new std::vector<float>();  
  t_FatJet_n3b1= new std::vector<float>();  
  t_FatJet_phi= new std::vector<float>();  
  t_FatJet_pt= new std::vector<float>();  
  t_FatJet_tau1= new std::vector<float>();  
  t_FatJet_tau2= new std::vector<float>();  
  t_FatJet_tau3= new std::vector<float>();  
  t_FatJet_tau4= new std::vector<float>();  
  t_FatJet_jetId= new std::vector<int>();  
  t_FatJet_subJetIdx1= new std::vector<int>();  
  t_FatJet_subJetIdx2= new std::vector<int>();  

  tree->Branch("t_FatJet_area",         "vector<float>", &t_FatJet_area);  
  tree->Branch("t_FatJet_particleNet_HbbvsQCD",         "vector<float>", &t_FatJet_particleNet_HbbvsQCD);  
  tree->Branch("t_FatJet_particleNet_mass",         "vector<float>", &t_FatJet_particleNet_mass);  
  tree->Branch("t_FatJet_btagDeepB",         "vector<float>", &t_FatJet_btagDeepB);  
  tree->Branch("t_FatJet_eta",         "vector<float>", &t_FatJet_eta);  
  tree->Branch("t_FatJet_mass",         "vector<float>", &t_FatJet_mass);  
  tree->Branch("t_FatJet_msoftdrop",         "vector<float>", &t_FatJet_msoftdrop);  
  tree->Branch("t_FatJet_n2b1",         "vector<float>", &t_FatJet_n2b1);  
  tree->Branch("t_FatJet_n3b1",         "vector<float>", &t_FatJet_n3b1);  
  tree->Branch("t_FatJet_phi",         "vector<float>", &t_FatJet_phi);  
  tree->Branch("t_FatJet_pt",         "vector<float>", &t_FatJet_pt);  
  tree->Branch("t_FatJet_tau1",         "vector<float>", &t_FatJet_tau1);  
  tree->Branch("t_FatJet_tau2",         "vector<float>", &t_FatJet_tau2);  
  tree->Branch("t_FatJet_tau3",         "vector<float>", &t_FatJet_tau3);  
  tree->Branch("t_FatJet_tau4",         "vector<float>", &t_FatJet_tau4);  
  tree->Branch("t_FatJet_jetId",         "vector<int>", &t_FatJet_jetId);  
  tree->Branch("t_FatJet_subJetIdx1",         "vector<int>", &t_FatJet_subJetIdx1);  
  tree->Branch("t_FatJet_subJetIdx2",         "vector<int>", &t_FatJet_subJetIdx2);  
  
  t_SubJet_btagCMVA= new std::vector<float>();   
  t_SubJet_btagCSVV2= new std::vector<float>();   
  t_SubJet_btagDeepB= new std::vector<float>();   
  t_SubJet_eta= new std::vector<float>();   
  t_SubJet_mass= new std::vector<float>();   
  t_SubJet_n2b1= new std::vector<float>();   
  t_SubJet_n3b1= new std::vector<float>();   
  t_SubJet_phi= new std::vector<float>();   
  t_SubJet_pt= new std::vector<float>();   
  t_SubJet_tau1= new std::vector<float>();   
  t_SubJet_tau2= new std::vector<float>();   
  t_SubJet_tau3= new std::vector<float>();   
  t_SubJet_tau4= new std::vector<float>();   
  
  tree->Branch("t_SubJet_btagCMVA",        "vector<float>", &t_SubJet_btagCMVA);   
  tree->Branch("t_SubJet_btagCSVV2",        "vector<float>", &t_SubJet_btagCSVV2);   
  tree->Branch("t_SubJet_btagDeepB",        "vector<float>", &t_SubJet_btagDeepB);   
  tree->Branch("t_SubJet_eta",        "vector<float>", &t_SubJet_eta);   
  tree->Branch("t_SubJet_mass",        "vector<float>", &t_SubJet_mass);   
  tree->Branch("t_SubJet_n2b1",        "vector<float>", &t_SubJet_n2b1);   
  tree->Branch("t_SubJet_n3b1",        "vector<float>", &t_SubJet_n3b1);   
  tree->Branch("t_SubJet_phi",        "vector<float>", &t_SubJet_phi);   
  tree->Branch("t_SubJet_pt",        "vector<float>", &t_SubJet_pt);   
  tree->Branch("t_SubJet_tau1",        "vector<float>", &t_SubJet_tau1);   
  tree->Branch("t_SubJet_tau2",        "vector<float>", &t_SubJet_tau2);   
  tree->Branch("t_SubJet_tau3",        "vector<float>", &t_SubJet_tau3);   
  tree->Branch("t_SubJet_tau4",        "vector<float>", &t_SubJet_tau4);   
   

  t_Jet_area= new std::vector<float>();   
  t_Jet_btagDeepFlavB= new std::vector<float>();   
  t_Jet_chEmEF= new std::vector<float>();   
  t_Jet_chHEF= new std::vector<float>();   
  t_Jet_eta= new std::vector<float>();   
  t_Jet_mass= new std::vector<float>();   
  t_Jet_neEmEF= new std::vector<float>();   
  t_Jet_neHEF= new std::vector<float>();   
  t_Jet_phi= new std::vector<float>();   
  t_Jet_pt= new std::vector<float>();   
  t_Jet_qgl= new std::vector<float>();   
  t_Jet_jetId= new std::vector<int>();   
  t_Jet_nConstituents= new std::vector<int>();   
  t_Jet_nElectrons= new std::vector<int>();   
  t_Jet_nMuons= new std::vector<int>();   
  t_Jet_puId= new std::vector<int>();   

  tree->Branch("t_nJet",  &t_nJet,"t_nJet/I");
  tree->Branch("t_Jet_area"    , "vector<float>"         ,&t_Jet_area);   
  tree->Branch("t_Jet_btagDeepFlavB"    , "vector<float>"         ,&t_Jet_btagDeepFlavB);   
  tree->Branch("t_Jet_chEmEF"    , "vector<float>"         ,&t_Jet_chEmEF);   
  tree->Branch("t_Jet_chHEF"    , "vector<float>"         ,&t_Jet_chHEF);   
  tree->Branch("t_Jet_eta"    , "vector<float>"         ,&t_Jet_eta);   
  tree->Branch("t_Jet_mass"    , "vector<float>"         ,&t_Jet_mass);   
  tree->Branch("t_Jet_neEmEF"    , "vector<float>"         ,&t_Jet_neEmEF);   
  tree->Branch("t_Jet_neHEF"    , "vector<float>"         ,&t_Jet_neHEF);   
  tree->Branch("t_Jet_phi"    , "vector<float>"         ,&t_Jet_phi);   
  tree->Branch("t_Jet_pt"    , "vector<float>"         ,&t_Jet_pt);   
  tree->Branch("t_Jet_qgl"    , "vector<float>"         ,&t_Jet_qgl);   
  tree->Branch("t_Jet_jetId"    , "vector<int>"         ,&t_Jet_jetId);   
  tree->Branch("t_Jet_nConstituents"    , "vector<int>"         ,&t_Jet_nConstituents);   
  tree->Branch("t_Jet_nElectrons"    , "vector<int>"         ,&t_Jet_nElectrons);   
  tree->Branch("t_Jet_nMuons"    , "vector<int>"         ,&t_Jet_nMuons);   
  tree->Branch("t_Jet_puId"    , "vector<int>"         ,&t_Jet_puId);   

  tree->Branch("t_nbJet",  &t_nbJet,"t_nbJet/I");
  t_bJet_area= new std::vector<float>();     
  t_bJet_btagDeepFlavB= new std::vector<float>();   
  t_bJet_chEmEF= new std::vector<float>();   
  t_bJet_chHEF= new std::vector<float>();   
  t_bJet_eta= new std::vector<float>();   
  t_bJet_mass= new std::vector<float>();   
  t_bJet_neEmEF= new std::vector<float>();   
  t_bJet_neHEF= new std::vector<float>();   
  t_bJet_phi= new std::vector<float>();   
  t_bJet_pt= new std::vector<float>();   
  t_bJet_qgl= new std::vector<float>();   
  t_bJet_jetId= new std::vector<int>();   
  t_bJet_nConstituents= new std::vector<int>();   
  t_bJet_nElectrons= new std::vector<int>();   
  t_bJet_nMuons= new std::vector<int>();   
  t_bJet_puId= new std::vector<int>();   
  t_bJet_SF= new std::vector<double>();
  t_bJet_SFup= new std::vector<double>();
  t_bJet_SFdown= new std::vector<double>();
  
  tree->Branch("t_bJet_area"    , "vector<float>"         ,&t_bJet_area);    
  tree->Branch("t_bJet_btagDeepFlavB"    , "vector<float>"         ,&t_bJet_btagDeepFlavB);   
  tree->Branch("t_bJet_chEmEF"    , "vector<float>"         ,&t_bJet_chEmEF);   
  tree->Branch("t_bJet_chHEF"    , "vector<float>"         ,&t_bJet_chHEF);   
  tree->Branch("t_bJet_eta"    , "vector<float>"         ,&t_bJet_eta);   
  tree->Branch("t_bJet_mass"    , "vector<float>"         ,&t_bJet_mass);   
  tree->Branch("t_bJet_neEmEF"    , "vector<float>"         ,&t_bJet_neEmEF);   
  tree->Branch("t_bJet_neHEF"    , "vector<float>"         ,&t_bJet_neHEF);   
  tree->Branch("t_bJet_phi"    , "vector<float>"         ,&t_bJet_phi);   
  tree->Branch("t_bJet_pt"    , "vector<float>"         ,&t_bJet_pt);   
  tree->Branch("t_bJet_qgl"    , "vector<float>"         ,&t_bJet_qgl);   
  tree->Branch("t_bJet_jetId"    , "vector<int>"         ,&t_bJet_jetId);   
  tree->Branch("t_bJet_nConstituents"    , "vector<int>"         ,&t_bJet_nConstituents);   
  tree->Branch("t_bJet_nElectrons"    , "vector<int>"         ,&t_bJet_nElectrons);   
  tree->Branch("t_bJet_nMuons"    , "vector<int>"         ,&t_bJet_nMuons);   
  tree->Branch("t_bJet_puId"    , "vector<int>"         ,&t_bJet_puId);   
  tree->Branch("t_bJet_SF"    , "vector<double>"         ,&t_bJet_SF);
  tree->Branch("t_bJet_SFup"    , "vector<double>"         ,&t_bJet_SFup);
  tree->Branch("t_bJet_SFdown"    , "vector<double>"     ,&t_bJet_SFdown);


  tree->Branch("t_PV_ndof", &t_PV_ndof, "t_PV_ndof/F");
  tree->Branch("t_PV_x", &t_PV_x,"t_PV_x/F");
  tree->Branch("t_PV_y", &t_PV_y,"t_PV_y/F");
  tree->Branch("t_PV_z", &t_PV_z,"t_PV_z/F");
  tree->Branch("t_PV_npvs",&t_PV_npvs,"t_PV_npvs/I");  
  tree->Branch("t_PV_npvsGood",&t_PV_npvsGood,"t_PV_npvsGood/I"); 
    
  tree->Branch("genweight", &genweight, "genweight/F");
  tree->Branch("cat", &cat, "cat/I");
  tree->Branch("HHH_m", &HHH_m, "HHH_m/F");  
  tree->Branch("HHH_pt", &HHH_pt, "HHH_pt/F");
  tree->Branch("HHH_eta", &HHH_eta, "HHH_eta/F");
  tree->Branch("HHH_phi", &HHH_phi, "HHH_phi/F");
  tree->Branch("HH12_m", &HH12_m, "HH12_m/F");
  tree->Branch("HH12_pt", &HH12_pt, "HH12_pt/F");
  tree->Branch("HH12_eta", &HH12_eta, "HH12_eta/F");
  tree->Branch("HH12_phi", &HH12_phi, "HH12_phi/F");
  tree->Branch("HH23_m", &HH23_m, "HH23_m/F");
  tree->Branch("HH23_pt", &HH23_pt, "HH23_ptt/F");
  tree->Branch("HH23_eta", &HH23_eta, "HH23_eta/F");
  tree->Branch("HH23_phi", &HH23_phi, "HH23_phi/F");
  tree->Branch("HH13_m", &HH13_m, "HH13_m/F");
  tree->Branch("HH13_pt", &HH13_pt, "HH13_pt/F");
  tree->Branch("HH13_eta", &HH13_eta, "HH13_eta/F");
  tree->Branch("HH13_phi", &HH13_phi, "HH13_phi/F");
  tree->Branch("H1_m", &H1_m, "H1_m/F");
  tree->Branch("H1_pt", &H1_pt, "H1_pt/F");
  tree->Branch("H1_eta", &H1_eta, "H1_eta/F");
  tree->Branch("H1_phi", &H1_phi, "H1_phi/F");
  tree->Branch("H2_m", &H2_m, "H2_m/F");
  tree->Branch("H2_pt", &H2_pt, "H2_pt/F");
  tree->Branch("H2_eta", &H2_eta, "H2_eta/F");
  tree->Branch("H2_phi", &H2_phi, "H2_phi/F");
  tree->Branch("H3_m", &H3_m, "H3_m/F");
  tree->Branch("H3_pt", &H3_pt, "H3_pt/F");
  tree->Branch("H3_eta", &H3_eta, "H3_eta/F");
  tree->Branch("H3_phi", &H3_phi, "H3_phi/F");
  tree->Branch("bJet1_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet2_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet3_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet4_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet5_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet6_btagDeepFlavB", &bJet1_btagDeepFlavB, "bJet1_btagDeepFlavB/F");
  tree->Branch("bJet1_pt", &bJet1_pt, "bJet1_pt/F");
  tree->Branch("bJet2_pt", &bJet2_pt, "bJet2_pt/F");
  tree->Branch("bJet3_pt", &bJet3_pt, "bJet3_pt/F");
  tree->Branch("bJet4_pt", &bJet4_pt, "bJet4_pt/F");
  tree->Branch("bJet5_pt", &bJet5_pt, "bJet5_pt/F");
  tree->Branch("bJet6_pt", &bJet6_pt, "bJet6_pt/F");
  tree->Branch("bJet1_eta", &bJet1_eta, "bJet1_eta/F");
  tree->Branch("bJet2_eta", &bJet2_eta, "bJet2_eta/F");
  tree->Branch("bJet3_eta", &bJet3_eta, "bJet3_eta/F");
  tree->Branch("bJet4_eta", &bJet4_eta, "bJet4_eta/F");
  tree->Branch("bJet5_eta", &bJet5_eta, "bJet5_eta/F");
  tree->Branch("bJet6_eta", &bJet6_eta, "bJet6_eta/F");

}
#endif // #ifdef HHHAnalyzer_cxx
