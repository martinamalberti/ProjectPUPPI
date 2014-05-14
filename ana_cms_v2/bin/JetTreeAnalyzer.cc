#include "./JetTreeAnalyzer.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"

#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


JetTreeAnalyzer::JetTreeAnalyzer(TTree *tree){
  Init(tree);
}

JetTreeAnalyzer::~JetTreeAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void JetTreeAnalyzer::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize                                                                                                                                        
  // a new tree or chain. Typically here the branch addresses and branch                                                                                                                                        
  // pointers of the tree will be set.                                                                                                                                                                          
  // It is normally not necessary to make changes to the generated                                                                                                                                              
  // code, but the routine can be extended by the user if needed.                                                                                                                                               
  // Init() will be called many times when running on PROOF                                                                                                                                                     
  // (once per file to be processed).                                                                                                                                                                           

  // Set object pointer                                                                                                                                                                                         
  jetUncorrMass = 0;
  jetUncorrPt = 0;
  jetUncorrEta = 0;
  jetUncorrPhi = 0;
  jetTrimmedMass = 0;
  jetMass = 0;
  jetPt = 0;
  jetEta = 0;
  jetPhi = 0;
  jetNparticles = 0;
  jetGenMacthIndex = 0;
  // Set branch addresses and branch pointers                                                                                                                                                                   
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("njetsUncorr", &njetsUncorr, &b_njetsUncorr);
  fChain->SetBranchAddress("njets", &njets, &b_njets);
  fChain->SetBranchAddress("sumEt", &sumEt, &b_sumEt);
  fChain->SetBranchAddress("etMissX", &etMissX, &b_etMissX);
  fChain->SetBranchAddress("etMissY", &etMissY, &b_etMissY);
  fChain->SetBranchAddress("jetUncorrMass", &jetUncorrMass, &b_jetUncorrMass);
  fChain->SetBranchAddress("jetUncorrPt", &jetUncorrPt, &b_jetUncorrPt);
  fChain->SetBranchAddress("jetUncorrEta", &jetUncorrEta, &b_jetUncorrEta);
  fChain->SetBranchAddress("jetUncorrPhi", &jetUncorrPhi, &b_jetUncorrPhi);
  fChain->SetBranchAddress("jetTrimmedMass", &jetTrimmedMass, &b_jetTrimmedMass);
  fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
  fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
  fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
  fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
  fChain->SetBranchAddress("jetNparticles", &jetNparticles, &b_jetNparticles);
  fChain->SetBranchAddress("jetGenMacthIndex", &jetGenMacthIndex, &b_jetGenMacthIndex);
}


Int_t JetTreeAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                                                                                                                                           
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


void JetTreeAnalyzer::bookHistograms(std::string suffix){

  std::cout << "Booking histograms for " << suffix.c_str() << " tree" << std::endl;

  hnjets = new TH1F(("hnjets_"+suffix).c_str(), ("hnjets_"+suffix).c_str(), 50, 0, 50 );

  hrawpt = new TH1F(("hrawpt_"+suffix).c_str(), ("hrawpt_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_pu = new TH1F(("hrawpt_pu_"+suffix).c_str(), ("hrawpt_pu_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_good = new TH1F(("hrawpt_good_"+suffix).c_str(), ("hrawpt_good_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_response = new TH1F(("hrawpt_response_"+suffix).c_str(), ("hrawpt_response_"+suffix).c_str(), 100, -100, 100 );

  hpt = new TH1F(("hpt_"+suffix).c_str(), ("hpt_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu = new TH1F(("hpt_pu_"+suffix).c_str(), ("hpt_pu_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good = new TH1F(("hpt_good_"+suffix).c_str(), ("hpt_good_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response = new TH1F(("hpt_response_"+suffix).c_str(), ("hpt_response_"+suffix).c_str(), 100, -100, 100 );

  hrawpt_leadjet = new TH1F(("hrawpt_leadjet_"+suffix).c_str(), ("hrawpt_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_pu_leadjet = new TH1F(("hrawpt_pu_leadjet_"+suffix).c_str(), ("hrawpt_pu_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_good_leadjet = new TH1F(("hrawpt_good_leadjet_"+suffix).c_str(), ("hrawpt_good_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_response_leadjet = new TH1F(("hrawpt_response_leadjet_"+suffix).c_str(), ("hrawpt_response_leadjet_"+suffix).c_str(), 100, -100, 100 );

  hpt_leadjet = new TH1F(("hpt_leadjet_"+suffix).c_str(), ("hpt_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu_leadjet = new TH1F(("hpt_pu_leadjet_"+suffix).c_str(), ("hpt_pu_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good_leadjet = new TH1F(("hpt_good_leadjet_"+suffix).c_str(), ("hpt_good_leadjet_"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response_leadjet = new TH1F(("hpt_response_leadjet_"+suffix).c_str(), ("hpt_response_leadjet_"+suffix).c_str(), 100, -100, 100 );
  
  heta = new TH1F(("heta_"+suffix).c_str(), ("heta_"+suffix).c_str(), 100, -5, 5 );
  heta_pu = new TH1F(("heta_pu_"+suffix).c_str(), ("heta_pu_"+suffix).c_str(), 100, -5, 5 );
  heta_good = new TH1F(("heta_good_"+suffix).c_str(), ("heta_good_"+suffix).c_str(), 100, -5, 5 );

  heta_leadjet = new TH1F(("heta_leadjet_"+suffix).c_str(), ("heta_leadjet_"+suffix).c_str(), 100, -5, 5 );
  heta_pu_leadjet = new TH1F(("heta_pu_leadjet_"+suffix).c_str(), ("heta_pu_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_good_leadjet = new TH1F(("heta_good_leadjet_"+suffix).c_str(), ("heta_good_leadjet_"+suffix).c_str(), 100, -5, 5 );

  hmass = new TH1F(("hmass_"+suffix).c_str(), ("hmass_"+suffix).c_str(), 100, 0, 100 );
  hmass_response = new TH1F(("hmass_response_"+suffix).c_str(), ("hmass_response_"+suffix).c_str(), 100, -5, 5 );

  hnparticles = new TH1F(("hnparticles_"+suffix).c_str(), ("hnparticles_"+suffix).c_str(), 100, 0, 100 );

}


void JetTreeAnalyzer::fillHistograms(int maxEntries, float minPt){

  std::cout << "Filling histograms..." << std::endl;
  
  if (fChain==0){
    std::cout<<"Error: cannot open " << fChain->GetName() << std::endl;
    exit(0);
  }

  if (maxEntries == -1)
    maxEntries = fChain->GetEntries();

  
  int njets = 0;

  for (int entry = 0; entry < maxEntries; entry++){
    fChain->GetEntry(entry);
    
    if (entry%100==0) std::cout << "Analyzing entry : " << entry << "\r" << std::flush;
    
    // --- Loop over jets in this event                                                                                                                                                                       
    for (int j = 0; j < njetsUncorr; j++){
      
      float pt = jetUncorrPt->at(j);
      
      if (pt < minPt)  continue;

      njets++;
      int imatch = (jetGenMacthIndex)->at(j);
      hrawpt->Fill(jetUncorrPt->at(j));
      if (imatch == -1) hrawpt_pu->Fill(jetUncorrPt->at(j));
      else hrawpt_good->Fill(jetUncorrPt->at(j));
      
    }// end loop over jets 

    hnjets->Fill(float(njets));

  }// end loop over entries

}


void JetTreeAnalyzer::saveHistograms(TFile *outfile){

  std::cout << "Saving histograms ... " << std::endl;
  
  outfile->cd();
  
  hnjets->Write();
  hrawpt->Write();
  hrawpt_pu->Write();
  hrawpt_good->Write();

}
