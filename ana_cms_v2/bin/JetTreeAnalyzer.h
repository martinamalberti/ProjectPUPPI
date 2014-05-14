#ifndef JetTreeAnalyzer_h
#define JetTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>

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

class JetTreeAnalyzer{

 public :

  JetTreeAnalyzer(TTree *tree);
  virtual ~JetTreeAnalyzer();

  virtual void Init(TTree *tree);
  virtual int  GetEntry(Long64_t entry);

  virtual void bookHistograms(std::string suffix);
  virtual void fillHistograms(int maxEntries, float minPt);
  virtual void saveHistograms(TFile *file);

  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain                                                                                                                                          
  Int_t           fCurrent; //!current Tree number in a TChain                                                                                                                                                  

  // Declaration of leaf types                                                                                                                                                                                  
  Int_t           njetsUncorr;
  Int_t           njets;
  Float_t         sumEt;
  Float_t         etMissX;
  Float_t         etMissY;
  std::vector<float>   *jetUncorrMass;
  std::vector<float>   *jetUncorrPt;
  std::vector<float>   *jetUncorrEta;
  std::vector<float>   *jetUncorrPhi;
  std::vector<float>   *jetTrimmedMass;
  std::vector<float>   *jetMass;
  std::vector<float>   *jetPt;
  std::vector<float>   *jetEta;
  std::vector<float>   *jetPhi;
  std::vector<float>   *jetNparticles;
  std::vector<int>     *jetGenMacthIndex;

  // List of branches                                                                                                                                                                                           
  TBranch        *b_njetsUncorr;   //!                                                                                                                                                                          
  TBranch        *b_njets;   //!                                                                                                                                                                                
  TBranch        *b_sumEt;   //!                                                                                                                                                                                
  TBranch        *b_etMissX;   //!                                                                                                                                                                              
  TBranch        *b_etMissY;   //!                                                                                                                                                                              
  TBranch        *b_jetUncorrMass;   //!                                                                                                                                                                        
  TBranch        *b_jetUncorrPt;   //!                                                                                                                                                                          
  TBranch        *b_jetUncorrEta;   //!                                                                                                                                                                         
  TBranch        *b_jetUncorrPhi;   //!                                                                                                                                                                         
  TBranch        *b_jetTrimmedMass;   //!                                                                                                                                                                       
  TBranch        *b_jetMass;   //!                                                                                                                                                                              
  TBranch        *b_jetPt;   //!                                                                                                                                                                                
  TBranch        *b_jetEta;   //!                                                                                                                                                                               
  TBranch        *b_jetPhi;   //!                                                                                                                                                                               
  TBranch        *b_jetNparticles;   //!                                                                                                                                                                        
  TBranch        *b_jetGenMacthIndex;   //!    


  // histograms declaration
  TH1F *hnjets;

  TH1F* hrawpt;
  TH1F* hrawpt_pu;
  TH1F* hrawpt_good;
  TH1F* hrawpt_response;
  TH1F* hpt;
  TH1F* hpt_pu;
  TH1F* hpt_good;
  TH1F* hpt_response;

  TH1F* hrawpt_leadjet;
  TH1F* hrawpt_pu_leadjet;
  TH1F* hrawpt_good_leadjet;
  TH1F* hrawpt_response_leadjet;
  TH1F* hpt_leadjet;
  TH1F* hpt_pu_leadjet;
  TH1F* hpt_good_leadjet;
  TH1F* hpt_response_leadjet;

  TH1F* heta;
  TH1F* heta_pu;
  TH1F* heta_good;
  TH1F* heta_leadjet;
  TH1F* heta_pu_leadjet;
  TH1F* heta_good_leadjet;

  TH1F* hmass;
  TH1F* hmass_response;

  TH1F* hnparticles;

 private:

};
#endif
