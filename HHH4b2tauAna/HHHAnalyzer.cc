#define HHHAnalyzer_cxx
#include "HHHAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

#ifdef MAKECINT
#pragma link C++ class vector<float>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<int>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<bool>+;
#endif

bool myfunction (int i,int j) { return (i>j); }

int main(int argc, char* argv[])
{

  if(argc < 7) {
    cerr << "Please give 7 arguments: inputFileList / outputFileName / samplename / data type / year / produce gen info or not/ skim cut "<<endl;
    return -1;
  }   
      
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *samplename    = argv[3];
  const char *isData        = argv[4];
  const char *isRunGen      = argv[5];
  TString year_num          = argv[6];
    
  HHHAnalyzer Hmm(inputFileList, outFileName, samplename, isData, year_num);
  cout << "dataset " << samplename << " year " <<year_num<< endl;

  if(argc == 7){
      cout <<"ana events"<<endl;
      //Hmm.cal_sumOfgw(samplename, isData);
      Hmm.EventLoop(samplename, isData, isRunGen);
  }
  else{
      cout <<"skim nanoAOD"<<endl;
      const char *skim_cut = argv[7];
      Hmm.SkimTChain(skim_cut);
  } 

  return 0;
}

void HHHAnalyzer::SkimTChain(const char *skim_cut)
{ 
  //TTree* treeskim = (TTree*)fChain->CopyTree(skim_cut);
  (TTree*)fChain->CopyTree(skim_cut);
}

void HHHAnalyzer::cal_sumOfgw(string samplename, const char *isData)
{   
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout <<"total entries: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    double sum = 0.;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%10000==0) cout <<"entry: "<<jentry<<endl;
        sum = sum + genWeight;
    }
    h_sumOfgw->SetBinContent(1,sum);
}

void HHHAnalyzer::EventLoop(string samplename,const char *isData, const char *isRunGen)
{ 
  if (fChain == 0) return;
  //clearTreeVectors();
  //cout<<"cleared tree vectors\n";
  //BookTreeBranches();
  //cout<<"booked tree branches\n";
    
  float muon_mass = 0.1056583745;
  float lumi = luminosity[yearst]; //fb-1
  cout <<"year "<<yearst<<" lumi "<<lumi<<endl;
  bool datafile = true;
  if(*isData=='F') datafile = false;
  bool debug = false;
    
  double sumOfgenweight = 1.0;
  if(!datafile){
        sumOfgenweight = sumOfgenw[yearst][samplename];
        cout <<"samplename: "<<samplename<< " xs: "<<xs[samplename]<<" sumOfgenweight "<<sumOfgenweight<<endl;
  }
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%1==0) cout <<"entry: "<<jentry<<endl;
      clearTreeVectors();
      //fill event weight
      genweight = genWeight*xs[samplename]*lumi/sumOfgenweight;
     
      //trigger selection
      bool trig_decision = false;
      //if( year=="2016" && (HLT_QuadJet45_TripleBTagCSV_OR_HLT_DoubleJet90_Double30_TripleBTagCSV || HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 || HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20 || HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20 || HLT_AK8PFJet360_TrimMass30 || HLT_AK8PFJet450 || HLT_PFJet450)) trig_decision =true;
      
      //if( year=="2017" && (HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0 || HLT_PFJet450 || HLT_PFJet500 || HLT_AK8PFJet500 || HLT_PFHT1050 || HLT_AK8PFJet360_TrimMass30 || HLT_AK8PFJet380_TrimMass30 || HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFHT800_TrimMass50 || HLT_AK8PFHT750_TrimMass50 || HLT_AK8PFJet330_PFAK8BTagCSV_p17)) trig_decision =true;
      
      //if( year=="2018" && (HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV || HLT_PFHT1050 || HLT_PFJet500 || HLT_AK8PFJet500 || HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFHT800_TrimMass50 || HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)) trig_decision =true;
              
      t_run =run;
      t_luminosityBlock=luminosityBlock;
      t_event=event;
      //cout<<jentry<<" : "<<t_event<<"-------------------\n";
      
      if(((nFatJet>=1 && nJet>=2) || (nJet>=4)) && (nTau>=2)){
          //deepJet https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
          //https://twiki.cern.ch/twiki/bin/view/CMS/DeepJet
          vector<int> JetIndex, bJetIndex, FatJetIndex, TauIndex;
          JetIndex.clear();
          bJetIndex.clear();
          FatJetIndex.clear();
          TauIndex.clear();
          for(int j=0;j<nJet;j++){
            if(Jet_pt[j]>30. && fabs(Jet_eta[j])<2.5 && Jet_jetId[j]>=2 && ((Jet_puId[j]>=4 && Jet_pt[j]<50.) || (Jet_pt[j]>50.))){
               t_nJet++;
               JetIndex.push_back(j);
               if(Jet_btagDeepFlavB[j]>btag_cut[yearst]){
                   t_nbJet++;
                   bJetIndex.push_back(j);
               }
            }
          }//end of njet
    /*
	if(t_Jet_pt->size()>=2){
          for(int k=0;k<t_Jet_pt->size();k++){
            for(int m=k+1; m<t_Jet_pt->size();m++){
              TLorentzVector j1,j2, jj;
              if(k==0 && m==1){
          	  j1.SetPtEtaPhiM((*t_Jet_pt)[0], (*t_Jet_eta)[0],(*t_Jet_phi)[0],(*t_Jet_mass)[0]);
	          j2.SetPtEtaPhiM((*t_Jet_pt)[1], (*t_Jet_eta)[1],(*t_Jet_phi)[1],(*t_Jet_mass)[1]);
	          jj=j1+j2;	
	          t_diJet_pt = jj.Pt();
	          t_diJet_eta=jj.Eta();
	          t_diJet_phi=jj.Phi();
	          t_diJet_mass=jj.M();
                  t_diJet_mass_mo=jj.M();
              }
              else{
                  j1.SetPtEtaPhiM((*t_Jet_pt)[k], (*t_Jet_eta)[k],(*t_Jet_phi)[k],(*t_Jet_mass)[k]);
                  j2.SetPtEtaPhiM((*t_Jet_pt)[m], (*t_Jet_eta)[m],(*t_Jet_phi)[m],(*t_Jet_mass)[m]);
                  jj=j1+j2;
                  if(t_diJet_mass_mo<jj.M()) t_diJet_mass_mo=jj.M();

             }
            }
          }
	}
    */
    //cout <<"t_nbJet and nFat jet before removal: "<<t_nbJet<<" "<<nFatJet<<endl;
    vector<int> bJetIndex_erase;
    bJetIndex_erase.clear();   
	for(int i=0;i<nFatJet;i++){
      //cout <<"fatjet "<<i<<endl;
      double xbb = FatJet_particleNetMD_Xbb[i]/(FatJet_particleNetMD_Xbb[i]+FatJet_particleNetMD_QCD[i]);
      if(FatJet_pt[i]>300. && fabs(FatJet_eta[i])<2.4 && FatJet_jetId[i]>=2 && FatJet_particleNet_mass[i]>50. && xbb>0.8){
          FatJetIndex.push_back(i);   
          //cout <<"fatjet selected "<<i<<endl;
          //remove resoved btag jets overlaping with AK8 jet
          for(int k=0; k<t_nbJet; k++){ 
              //cout <<"bjet "<<k<<" "<<bJetIndex.at(k)<<endl;
              if(DeltaR(FatJet_eta[i], FatJet_phi[i], Jet_eta[bJetIndex.at(k)], Jet_phi[bJetIndex.at(k)])< 0.4){
                  bJetIndex_erase.push_back(k);                  
              }
          }//for each fat jet passing selection, record AK4jet overlap with it to be removed
      }//fat jet selection
      /*
	  t_FatJet_area->push_back(FatJet_area[i]);  
	  t_FatJet_particleNet_HbbvsQCD->push_back(FatJet_particleNetMD_Xbb[i]/(FatJet_particleNetMD_Xbb[i]+FatJet_particleNetMD_QCD[i]));  
	  t_FatJet_eta->push_back(FatJet_eta[i]);  
	  t_FatJet_particleNet_mass->push_back(FatJet_particleNet_mass[i]);  
	  t_FatJet_pt->push_back(FatJet_pt[i]);  
	  t_FatJet_tau1->push_back(FatJet_tau1[i]);  
	  t_FatJet_tau2->push_back(FatJet_tau2[i]);  
	  t_FatJet_tau3->push_back(FatJet_tau3[i]);  
	  t_FatJet_tau4->push_back(FatJet_tau4[i]);  
	  t_FatJet_jetId->push_back(FatJet_jetId[i]);  
	  t_FatJet_subJetIdx1->push_back(FatJet_subJetIdx1[i]);  
	  t_FatJet_subJetIdx2->push_back(FatJet_subJetIdx2[i]); 
      */
	}//loop over nfat jet
    sort(bJetIndex_erase.begin(), bJetIndex_erase.end());
    for(int k=bJetIndex_erase.size()-1; k>=0; k--){ 
        bJetIndex.erase(bJetIndex.begin()+bJetIndex_erase.at(k));
        //cout <<"erase a bjet"<<endl;
    }
    t_nbJet = bJetIndex.size();
    
    if (debug) cout <<"t_nbJet after removal: "<<t_nbJet<<endl;
          
    //selection of taus
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
    //https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorkingLegacyRun2 https://twiki.cern.ch/twiki/bin/view/CMS/Higgs2Tau
    //SVFit https://github.com/SVfit/ClassicSVfit/tree/fastMTT_21_06_2018
    //https://github.com/SVfit/ClassicSVfit/blob/fastMTT_21_06_2018/bin/testClassicSVfit.cc#L121-L131
    //https://cms-nanoaod-integration.web.cern.ch/integration/master/mc102X_doc.html#Tau
    //bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
    for(int i=0;i<nTau;i++){
        if(Tau_pt[i]>20. && fabs(Tau_eta[i])<2.3 && (Tau_idDeepTau2017v2p1VSjet[i]>=8)){
            TauIndex.push_back(i);
        }
    }
          
    //categorization
    if(FatJetIndex.size()>=1){
        FatJet1_xbb = FatJet_particleNetMD_Xbb[FatJetIndex.at(0)]/(FatJet_particleNetMD_Xbb[FatJetIndex.at(0)]+FatJet_particleNetMD_QCD[FatJetIndex.at(0)]);
        FatJet1_pt = FatJet_pt[FatJetIndex.at(0)];
        FatJet1_eta = FatJet_eta[FatJetIndex.at(0)];
        FatJet1_phi = FatJet_phi[FatJetIndex.at(0)];
        FatJet1_msoftdrop = FatJet_msoftdrop[FatJetIndex.at(0)];
        FatJet1_particleNet_mass = FatJet_particleNet_mass[FatJetIndex.at(0)];
        
        if(FatJetIndex.size()>=2){
            FatJet2_xbb = FatJet_particleNetMD_Xbb[FatJetIndex.at(1)]/(FatJet_particleNetMD_Xbb[FatJetIndex.at(1)]+FatJet_particleNetMD_QCD[FatJetIndex.at(1)]);
            FatJet2_pt = FatJet_pt[FatJetIndex.at(1)];
            FatJet2_eta = FatJet_eta[FatJetIndex.at(1)];
            FatJet2_phi = FatJet_phi[FatJetIndex.at(1)];
            FatJet2_msoftdrop = FatJet_msoftdrop[FatJetIndex.at(1)];
            FatJet2_particleNet_mass = FatJet_particleNet_mass[FatJetIndex.at(1)];
            
            if(FatJetIndex.size()>=3){
                FatJet3_xbb = FatJet_particleNetMD_Xbb[FatJetIndex.at(2)]/(FatJet_particleNetMD_Xbb[FatJetIndex.at(2)]+FatJet_particleNetMD_QCD[FatJetIndex.at(2)]);
                FatJet3_pt = FatJet_pt[FatJetIndex.at(2)];
                FatJet3_eta = FatJet_eta[FatJetIndex.at(2)];
                FatJet3_phi = FatJet_phi[FatJetIndex.at(2)];
                FatJet3_msoftdrop = FatJet_msoftdrop[FatJetIndex.at(2)];
                FatJet3_particleNet_mass = FatJet_particleNet_mass[FatJetIndex.at(2)];
            }
        }
    }
    if (debug) cout <<"debugg: before categorization"<<endl;
    //categorization
    TLorentzVector TV_H1, TV_H2, TV_H3, TV_HH12, TV_HH23, TV_HH13, TV_HHH;
    TLorentzVector TV_b1, TV_b2, TV_b3, TV_b4, TV_tau1, TV_tau2;
    
    if(TauIndex.size()>=2){ 
        TV_tau1.SetPtEtaPhiM(Tau_pt[TauIndex.at(0)], Tau_eta[TauIndex.at(0)], Tau_phi[TauIndex.at(0)], Tau_mass[TauIndex.at(0)]);
        TV_tau2.SetPtEtaPhiM(Tau_pt[TauIndex.at(1)], Tau_eta[TauIndex.at(1)], Tau_phi[TauIndex.at(1)], Tau_mass[TauIndex.at(1)]); 
        TV_H3 = TV_tau1 + TV_tau2; 
        if(FatJetIndex.size()>=1 && t_nbJet>=2){
            cat = 0;
            TV_H1.SetPtEtaPhiM(FatJet1_pt, FatJet1_eta, FatJet1_phi, FatJet1_particleNet_mass);  
            TV_b1.SetPtEtaPhiM(Jet_pt[bJetIndex[0]]*Jet_bRegCorr[bJetIndex[0]], Jet_eta[bJetIndex[0]], Jet_phi[bJetIndex[0]], Jet_mass[bJetIndex[0]]);
            TV_b2.SetPtEtaPhiM(Jet_pt[bJetIndex[1]]*Jet_bRegCorr[bJetIndex[1]], Jet_eta[bJetIndex[1]], Jet_phi[bJetIndex[1]], Jet_mass[bJetIndex[1]]);
            TV_H2 = TV_b1 + TV_b2; 
            
        }
        else if(FatJetIndex.size()==0 && t_nbJet>=4){
            cat = 1;           
            TV_b1.SetPtEtaPhiM(Jet_pt[bJetIndex[0]], Jet_eta[bJetIndex[0]], Jet_phi[bJetIndex[0]], Jet_mass[bJetIndex[0]]);
            int bjet2_i = -1, bjet3_i = -1;
            float chi2_min = 9999.; 
            for(int bjet2_index = 1; bjet2_index<4; bjet2_index++){
                TV_b2.SetPtEtaPhiM(Jet_pt[bJetIndex[bjet2_index]]*Jet_bRegCorr[bJetIndex[bjet2_index]], Jet_eta[bJetIndex[bjet2_index]], Jet_phi[bJetIndex[bjet2_index]], Jet_mass[bJetIndex[bjet2_index]]);
                int bjet3_index = -1;
                for(int tmp_index = 1; tmp_index<4; tmp_index++){
                    if(tmp_index!=bjet2_index){
                        bjet3_index = tmp_index;
                        TV_b3.SetPtEtaPhiM(Jet_pt[bJetIndex[bjet3_index]]*Jet_bRegCorr[bJetIndex[bjet3_index]], Jet_eta[bJetIndex[bjet3_index]], Jet_phi[bJetIndex[bjet3_index]], Jet_mass[bJetIndex[bjet3_index]]);
                    }
                }
                //cout <<"debuggg: "<<6-bjet2_index-bjet3_index<<" "<<bjet2_index<<" "<<bjet3_index<<endl;
                TV_b4.SetPtEtaPhiM(Jet_pt[bJetIndex[6-bjet2_index-bjet3_index]]*Jet_bRegCorr[bJetIndex[6-bjet2_index-bjet3_index]], Jet_eta[bJetIndex[6-bjet2_index-bjet3_index]], Jet_phi[bJetIndex[6-bjet2_index-bjet3_index]], Jet_mass[bJetIndex[6-bjet2_index-bjet3_index]]);         
                float chi2_tmp = ((TV_b1 + TV_b2).M() - 125.)*((TV_b1 + TV_b2).M() - 125.) + ((TV_b3 + TV_b4).M() - 125.)*((TV_b3 + TV_b4).M() - 125.);
                if((chi2_tmp<chi2_min) || (bjet2_i==-1)){
                    chi2_min = chi2_tmp;
                    bjet2_i = bjet2_index;
                    bjet3_i = bjet3_index;
                }
            }
            if (debug) cout <<"after pairing 4b "<<bjet2_i<<" "<<bjet3_i<<" "<<6-bjet2_i-bjet3_i<<endl;
            TV_b2.SetPtEtaPhiM(Jet_pt[bJetIndex[bjet2_i]]*Jet_bRegCorr[bJetIndex[bjet2_i]], Jet_eta[bJetIndex[bjet2_i]], Jet_phi[bJetIndex[bjet2_i]], Jet_mass[bJetIndex[bjet2_i]]);
            TV_b3.SetPtEtaPhiM(Jet_pt[bJetIndex[bjet3_i]]*Jet_bRegCorr[bJetIndex[bjet3_i]], Jet_eta[bJetIndex[bjet3_i]], Jet_phi[bJetIndex[bjet3_i]], Jet_mass[bJetIndex[bjet3_i]]);
            TV_b4.SetPtEtaPhiM(Jet_pt[bJetIndex[6-bjet2_i-bjet3_i]]*Jet_bRegCorr[bJetIndex[6-bjet2_i-bjet3_i]], Jet_eta[bJetIndex[6-bjet2_i-bjet3_i]], Jet_phi[bJetIndex[6-bjet2_i-bjet3_i]], Jet_mass[bJetIndex[6-bjet2_i-bjet3_i]]);               
            TV_H1 = TV_b1 + TV_b2;
            TV_H2 = TV_b3 + TV_b4;                       
        }
    }
    vector<double> reco_Higgs_pt;
    vector<TLorentzVector> TV_H;
    reco_Higgs_pt.clear();
    TV_H.clear();
    reco_Higgs_pt.push_back(TV_H1.Pt()); reco_Higgs_pt.push_back(TV_H2.Pt()); reco_Higgs_pt.push_back(TV_H3.Pt());
    TV_H.push_back(TV_H1); TV_H.push_back(TV_H2); TV_H.push_back(TV_H3);
    vector<pair<double, int>> reco_higgs_pt_index = sortArr(reco_Higgs_pt, 3);  
          
    TV_HH12 = TV_H.at(reco_higgs_pt_index.at(2).second)+TV_H.at(reco_higgs_pt_index.at(1).second);
    TV_HH23 = TV_H.at(reco_higgs_pt_index.at(1).second)+TV_H.at(reco_higgs_pt_index.at(0).second); 
    TV_HH13 = TV_H.at(reco_higgs_pt_index.at(0).second)+TV_H.at(reco_higgs_pt_index.at(2).second);  
    TV_HHH = TV_H1 + TV_H2 + TV_H3;
          
    HHH_m = TV_HHH.M();
    HHH_pt = TV_HHH.Pt();
    HHH_eta = TV_HHH.Eta();
    HHH_phi = TV_HHH.Phi();
          
    HH12_m = TV_HH12.M();
    HH12_pt = TV_HH12.Pt();
    HH12_eta = TV_HH12.Eta();
    HH12_phi = TV_HH12.Phi();
    HH23_m = TV_HH23.M();
    HH23_pt = TV_HH23.Pt();
    HH23_eta = TV_HH23.Eta();
    HH23_phi = TV_HH23.Phi();
    HH13_m = TV_HH13.M();
    HH13_pt = TV_HH13.Pt();
    HH13_eta = TV_HH13.Eta();
    HH13_phi = TV_HH13.Phi();
          
    H1_m = TV_H.at(reco_higgs_pt_index.at(2).second).M();
    H1_pt = TV_H.at(reco_higgs_pt_index.at(2).second).Pt();
    H1_eta = TV_H.at(reco_higgs_pt_index.at(2).second).Eta();
    H1_phi = TV_H.at(reco_higgs_pt_index.at(2).second).Phi();
    H2_m = TV_H.at(reco_higgs_pt_index.at(1).second).M();
    H2_pt = TV_H.at(reco_higgs_pt_index.at(1).second).Pt();
    H2_eta = TV_H.at(reco_higgs_pt_index.at(1).second).Eta();
    H2_phi = TV_H.at(reco_higgs_pt_index.at(1).second).Phi();
    H3_m = TV_H.at(reco_higgs_pt_index.at(0).second).M();
    H3_pt = TV_H.at(reco_higgs_pt_index.at(0).second).Pt();
    H3_eta = TV_H.at(reco_higgs_pt_index.at(0).second).Eta();
    H3_phi = TV_H.at(reco_higgs_pt_index.at(0).second).Phi(); 
          
    Hditau_m = TV_H3.M();
    Hditau_pt = TV_H3.Pt();
    Hditau_eta = TV_H3.Eta();
    Hditau_phi = TV_H3.Phi();
          
    if(cat==0){
        Hfatjet_m = TV_H1.M();
        Hfatjet_pt = TV_H1.Pt();
        Hfatjet_eta = TV_H1.Eta();
        Hfatjet_phi = TV_H1.Phi();
        Hdijet1_m = TV_H2.M();
        Hdijet1_pt = TV_H2.Pt();
        Hdijet1_eta = TV_H2.Eta();
        Hdijet1_phi = TV_H2.Phi();
    }
    else if(cat==1){
        if(TV_H2.Pt()>TV_H1.Pt()){
            Hdijet1_m = TV_H2.M();
            Hdijet1_pt = TV_H2.Pt();
            Hdijet1_eta = TV_H2.Eta();
            Hdijet1_phi = TV_H2.Phi();
            Hdijet2_m = TV_H1.M();
            Hdijet2_pt = TV_H1.Pt();
            Hdijet2_eta = TV_H1.Eta();
            Hdijet2_phi = TV_H1.Phi();
        }
        else{
            Hdijet1_m = TV_H1.M();
            Hdijet1_pt = TV_H1.Pt();
            Hdijet1_eta = TV_H1.Eta();
            Hdijet1_phi = TV_H1.Phi();
            Hdijet2_m = TV_H2.M();
            Hdijet2_pt = TV_H2.Pt();
            Hdijet2_eta = TV_H2.Eta();
            Hdijet2_phi = TV_H2.Phi();
        }
    }
          
    /*
	t_PV_ndof = PV_ndof;
	t_PV_x = PV_x;
	t_PV_y = PV_y;
	t_PV_z = PV_z;
	t_PV_npvs = PV_npvs;
	t_PV_npvsGood = PV_npvsGood;
    */          
    }//end of preselection (nFatJets, nJets)

    if (debug) cout <<"debugg: before gen"<<endl;
    if(*isRunGen=='T' && *isData=='F'){
        //gen information: 
        //pdg ID scheme: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD        
        vector<int> index_gen_higgs;
        vector<double> gen_higgs_pt;
        vector<int> index_gen_b;
          
        index_gen_higgs.clear();
        gen_higgs_pt.clear();
        index_gen_b.clear();
          
            if(debug)  cout <<"start to run gen"<<endl;
            for(int igenpart=0; igenpart<nGenPart; igenpart++){   
                //find gen Higgses
                if(fabs(GenPart_pdgId[igenpart])==25 && GenPart_status[igenpart]==62) index_gen_higgs.push_back(igenpart);
            }
            if(index_gen_higgs.size()!=3){
                cout <<"wrong number of Higgses "<<index_gen_higgs.size()<<endl;
                continue;
            }
            for(int jgenHiggs=0; jgenHiggs<3; jgenHiggs++){
                gen_higgs_pt.push_back(GenPart_pt[index_gen_higgs[jgenHiggs]]);
                //find gen b-quarks
                vector<int> tmp;
                for(int igenpart=0; igenpart<nGenPart; igenpart++){
                                
                    if(debug){
                        cout<<"particle Idx ID status "<<igenpart<<" "<<GenPart_pdgId[igenpart] <<" "<<GenPart_status[igenpart]
                        <<" IdxMother "<<GenPart_genPartIdxMother[igenpart]
                        <<" pt "<<GenPart_pt[igenpart]
                        <<" mass "<<GenPart_mass[igenpart]<<endl;
                        cout <<" mother pdgid"<< GenPart_pdgId[GenPart_genPartIdxMother[igenpart]]<< " "<<GenPart_pt[GenPart_genPartIdxMother[igenpart]] <<" "<<GenPart_eta[GenPart_genPartIdxMother[igenpart]]<<endl;
                    }
                    
                    if((fabs(GenPart_pdgId[igenpart])==5 && GenPart_status[igenpart]==23) || (fabs(GenPart_pdgId[igenpart])==15 && (GenPart_status[igenpart]==23 || GenPart_status[igenpart]==24 || GenPart_status[igenpart]==2))){
                        if(GenPart_genPartIdxMother[igenpart] == index_gen_higgs.at(jgenHiggs)){
                            tmp.push_back(igenpart);
                        }
                    }
                }
                if(tmp.size()==2){
                    if(GenPart_pt[tmp.at(0)] > GenPart_pt[tmp.at(1)]){
                      index_gen_b.push_back(tmp.at(0));
                      index_gen_b.push_back(tmp.at(1));
                    }
                    else{
                      index_gen_b.push_back(tmp.at(1));
                      index_gen_b.push_back(tmp.at(0));
                    }
                }
                else{
                    cout <<"error not two b quarks from Higgs decay but "<<tmp.size()<<endl; 
                }
                tmp.clear();
            }
            if(debug) cout <<"index_gen_b.size() "<<index_gen_b.size()<<endl;
        
        if(index_gen_b.size()==6){
            //sort the gen Higgs by increasing pT
            vector<pair<double, int>> gen_higgs_pt_index = sortArr(gen_higgs_pt, 3); 
             
            gen_b0_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[2].second*2)];
            gen_b0_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[2].second*2)];
            gen_b0_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[2].second*2)];
            gen_b0_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[2].second*2)];
            
            gen_b1_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[2].second*2+1)];
            gen_b1_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[2].second*2+1)];
            gen_b1_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[2].second*2+1)];
            gen_b1_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[2].second*2+1)];
                  
            gen_b2_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[1].second*2)];
            gen_b2_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[1].second*2)];
            gen_b2_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[1].second*2)];
            gen_b2_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[1].second*2)];
            
            gen_b3_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[1].second*2+1)];
            gen_b3_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[1].second*2+1)];
            gen_b3_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[1].second*2+1)];
            gen_b3_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[1].second*2+1)];
            
            gen_b4_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[0].second*2)];
            gen_b4_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[0].second*2)];
            gen_b4_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[0].second*2)];
            gen_b4_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[0].second*2)];
            
            gen_b5_pt = GenPart_pt[index_gen_b.at(gen_higgs_pt_index[0].second*2+1)];
            gen_b5_eta = GenPart_eta[index_gen_b.at(gen_higgs_pt_index[0].second*2+1)];
            gen_b5_phi = GenPart_phi[index_gen_b.at(gen_higgs_pt_index[0].second*2+1)];
            gen_b5_m = GenPart_mass[index_gen_b.at(gen_higgs_pt_index[0].second*2+1)];
            
            gen_higg0_pt = GenPart_pt[index_gen_higgs.at(gen_higgs_pt_index[2].second)];
            gen_higg0_eta = GenPart_eta[index_gen_higgs.at(gen_higgs_pt_index[2].second)];
            gen_higg0_phi = GenPart_phi[index_gen_higgs.at(gen_higgs_pt_index[2].second)];
            gen_higg0_m = GenPart_mass[index_gen_higgs.at(gen_higgs_pt_index[2].second)];
            
            gen_higg1_pt = GenPart_pt[index_gen_higgs.at(gen_higgs_pt_index[1].second)];
            gen_higg1_eta = GenPart_eta[index_gen_higgs.at(gen_higgs_pt_index[1].second)];
            gen_higg1_phi = GenPart_phi[index_gen_higgs.at(gen_higgs_pt_index[1].second)];
            gen_higg1_m = GenPart_mass[index_gen_higgs.at(gen_higgs_pt_index[1].second)];
            
            gen_higg2_pt = GenPart_pt[index_gen_higgs.at(gen_higgs_pt_index[0].second)];
            gen_higg2_eta = GenPart_eta[index_gen_higgs.at(gen_higgs_pt_index[0].second)];
            gen_higg2_phi = GenPart_phi[index_gen_higgs.at(gen_higgs_pt_index[0].second)];
            gen_higg2_m = GenPart_mass[index_gen_higgs.at(gen_higgs_pt_index[0].second)];
            
            vector<TLorentzVector> H;
            TLorentzVector HHH;
            for(int jgenHiggs=0; jgenHiggs<3; jgenHiggs++){
                TLorentzVector H_tmp;
                H_tmp.SetPtEtaPhiM(GenPart_pt[index_gen_higgs.at(jgenHiggs)], GenPart_eta[index_gen_higgs.at(jgenHiggs)], GenPart_phi[index_gen_higgs.at(jgenHiggs)], GenPart_mass[index_gen_higgs.at(jgenHiggs)]);
                H.push_back(H_tmp);
            }
            
            HHH = H.at(0) + H.at(1) + H.at(2);
            genHHH_pt = HHH.Pt();
            genHHH_phi = HHH.Phi();
            genHHH_eta = HHH.Eta();
            genHHH_m = HHH.M();  
                
            //cout <<"debug: "<<gen_b0_pt<<endl;
                
            }
        }
        //cout <<"debug: "<<gen_b0_pt<<endl;
        tree->Fill();

      }
      tree->Write();
   }
