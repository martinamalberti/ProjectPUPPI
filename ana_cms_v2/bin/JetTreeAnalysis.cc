#include "./JetTree.h"

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

void bookHistograms(std::vector<TH1F*>* h, std::string suffix);
void fillHistograms(TTree *tree, std::vector<TH1F*>* h, int maxEntries, float minPt);
//void saveHistograms(std::vector<TH1F*> h, TFile *outfile);
void saveHistograms(std::vector<TH1F*> h);

//-------------------------------------------------------
// MAIN
//-------------------------------------------------------

int main( int argc, char **argv ) {

  gROOT->ProcessLine("#include <vector>");

  TFile *inputFile = TFile::Open("../outtree_test.root");
  if (inputFile==0){
    std::cout<<"Error: cannot open " << inputFile->GetName() << std::endl;
    exit(0);
  }


  // puppi tree
  TTree *tree_puppi = (TTree *)inputFile->Get("tree_pf_puppi");

  // histograms booking
  std::vector<TH1F*> h_puppi; 
  bookHistograms(&h_puppi,"puppi");
  
  std::cout<< h_puppi[0]->GetName() << "  " << h_puppi[0]->GetEntries()<<std::endl;

  // fill histograms
  fillHistograms(tree_puppi, &h_puppi, 10, 25.);

  std::cout<<"n histos = " << h_puppi.size() <<std::endl;
  
  std::cout<< h_puppi[0]->GetName() <<std::endl;
  std::cout<< "  " << h_puppi[0]->GetEntries()<<std::endl;

  
  //write histograms in the output file
  //  TFile *outfile = new TFile("out.root","RECREATE");
  //  saveHistograms(h_puppi, outfile);
  //saveHistograms(h_puppi);
}




void bookHistograms(std::vector<TH1F*>* h, std::string suffix){

  std::cout << "Booking histograms for " << suffix.c_str() << " tree" << std::endl;

  h->push_back(new TH1F(("hnjets_"+suffix).c_str(), ("hnjets_"+suffix).c_str(), 50, 0, 50 ));
  h->push_back(new TH1F(("hrawPt_"+suffix).c_str(), ("hrawPt_"+suffix).c_str(), 2000, 0, 2000 ));
  h->push_back(new TH1F(("hrawPt_pu_"+suffix).c_str(), ("hrawPt_pu_"+suffix).c_str(), 2000, 0, 2000 ));
  h->push_back(new TH1F(("hrawPt_good_"+suffix).c_str(), ("hrawPt_good_"+suffix).c_str(), 2000, 0, 2000 ));
  
}


void fillHistograms(TTree *tree, std::vector<TH1F*>* h, int maxEntries, float minPt){

  std::cout << "Filling histograms..." << std::endl;
  
  if (tree==0){
    std::cout<<"Error: cannot open " << tree->GetName() << std::endl;
    exit(0);
  }

  JetTree t(tree);

  if (maxEntries == -1) 
    maxEntries = tree->GetEntries();

  int njets = 0;

  for (int entry = 0; entry < maxEntries; entry++){

    t.GetEntry(entry);

    //    if (entry%1==0) std::cout << "Analyzing entry : " << entry << "\r" << std::flush;
    if (entry%1==0) std::cout << "Analyzing entry : " << entry << std::endl;

    // --- Loop over jets in this event                                                                                                                                                                       
    for (int j = 0; j < t.njetsUncorr; j++){
      
      float pt = t.jetUncorrPt->at(j);
      
      if (pt < minPt)  continue;

      njets++;
      int imatch = (t.jetGenMacthIndex)->at(j);
      
      (*h)[1]->Fill(t.jetUncorrPt->at(j));
      if (imatch == -1) (*h)[2]->Fill(t.jetUncorrPt->at(j));
      else (*h)[3]->Fill(t.jetUncorrPt->at(j));
      
    }// end loop over jets 

    (*h)[0]->Fill(njets);
    
  }// end loop over events

  for (unsigned int i = 0; i < h->size(); i++){
    std::cout << (*h)[i]->GetEntries() <<std::endl;
  }


}



//void saveHistograms(std::vector<TH1F*> h, TFile *outfile){
void saveHistograms(std::vector<TH1F*> h){

  std::cout << "Saving histograms ... " << std::endl;
  
  //outfile->cd();
  std::cout << h.size() <<std::endl;
  std::cout << h[0]->GetEntries() <<std::endl;
  
  //for (unsigned int i = 0; i < h.size(); i++){
  //  if (h[i])    std::cout << h[i]->GetEntries() << std::endl;
  //  //h[i]->Write();
  //  std::cout<< "ciao" << std::endl;
  //}

}
