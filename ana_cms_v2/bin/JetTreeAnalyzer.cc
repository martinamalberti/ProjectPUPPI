#include "./JetTreeAnalyzer.h"

// --- constructor ------------------------------------------------------
JetTreeAnalyzer::JetTreeAnalyzer(TTree *tree){
  Init(tree);
}

// ----------------------------------------------------------------------
JetTreeAnalyzer::~JetTreeAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


// --- Init tree --------------------------------------------------------
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
  jetGenMatchIndex = 0;
  jetGenMass = 0;
  jetGenPt = 0;
  jetGenEta = 0;
  jetGenPhi = 0;
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
  fChain->SetBranchAddress("jetGenMatchIndex", &jetGenMatchIndex, &b_jetGenMatchIndex);
  fChain->SetBranchAddress("jetGenMass", &jetGenMass, &b_jetGenMass);
  fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
  fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
  fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
}

// --- get Tree entry ----------------------------------------------------------------
Int_t JetTreeAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                                                                                                                                           
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


// --- Book histograms ---------------------------------------------------------------
void JetTreeAnalyzer::bookHistograms(std::string suffix = ""){

  std::cout << "Booking histograms for " << suffix.c_str() << " tree" << std::endl;

  hnjets = new TH1F(("hnjets"+suffix).c_str(), ("hnjets"+suffix).c_str(), 50, 0, 50 );

  hrawpt = new TH1F(("hrawpt"+suffix).c_str(), ("hrawpt"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_pu = new TH1F(("hrawpt_pu"+suffix).c_str(), ("hrawpt_pu"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_good = new TH1F(("hrawpt_good"+suffix).c_str(), ("hrawpt_good"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_response = new TH1F(("hrawpt_response"+suffix).c_str(), ("hrawpt_response"+suffix).c_str(), 100, -100, 100 );

  hpt = new TH1F(("hpt"+suffix).c_str(), ("hpt"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu = new TH1F(("hpt_pu"+suffix).c_str(), ("hpt_pu"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good = new TH1F(("hpt_good"+suffix).c_str(), ("hpt_good"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response = new TH1F(("hpt_response"+suffix).c_str(), ("hpt_response"+suffix).c_str(), 100, -100, 100 );

  hrawpt_leadjet = new TH1F(("hrawpt_leadjet"+suffix).c_str(), ("hrawpt_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_pu_leadjet = new TH1F(("hrawpt_pu_leadjet"+suffix).c_str(), ("hrawpt_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_good_leadjet = new TH1F(("hrawpt_good_leadjet"+suffix).c_str(), ("hrawpt_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hrawpt_response_leadjet = new TH1F(("hrawpt_response_leadjet"+suffix).c_str(), ("hrawpt_response_leadjet"+suffix).c_str(), 100, -100, 100 );

  hpt_leadjet = new TH1F(("hpt_leadjet"+suffix).c_str(), ("hpt_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_pu_leadjet = new TH1F(("hpt_pu_leadjet"+suffix).c_str(), ("hpt_pu_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_good_leadjet = new TH1F(("hpt_good_leadjet"+suffix).c_str(), ("hpt_good_leadjet"+suffix).c_str(), 2000, 0, 2000 );
  hpt_response_leadjet = new TH1F(("hpt_response_leadjet"+suffix).c_str(), ("hpt_response_leadjet"+suffix).c_str(), 100, -100, 100 );
  
  heta = new TH1F(("heta"+suffix).c_str(), ("heta"+suffix).c_str(), 100, -5, 5 );
  heta_pu = new TH1F(("heta_pu"+suffix).c_str(), ("heta_pu"+suffix).c_str(), 100, -5, 5 );
  heta_good = new TH1F(("heta_good"+suffix).c_str(), ("heta_good"+suffix).c_str(), 100, -5, 5 );

  heta_leadjet = new TH1F(("heta_leadjet"+suffix).c_str(), ("heta_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_pu_leadjet = new TH1F(("heta_pu_leadjet"+suffix).c_str(), ("heta_pu_leadjet"+suffix).c_str(), 100, -5, 5 );
  heta_good_leadjet = new TH1F(("heta_good_leadjet"+suffix).c_str(), ("heta_good_leadjet"+suffix).c_str(), 100, -5, 5 );

  hmass = new TH1F(("hmass"+suffix).c_str(), ("hmass"+suffix).c_str(), 100, 0, 100 );
  hmass_response = new TH1F(("hmass_response"+suffix).c_str(), ("hmass_response"+suffix).c_str(), 100, -50, 50 );

  hmass_leadjet = new TH1F(("hmass_leadjet"+suffix).c_str(), ("hmass_leadjet"+suffix).c_str(), 100, 0, 100 );
  hmass_response_leadjet = new TH1F(("hmass_response_leadjet"+suffix).c_str(), ("hmass_response_leadjet"+suffix).c_str(), 100, -50, 50 );

  hnparticles = new TH1F(("hnparticles"+suffix).c_str(), ("hnparticles"+suffix).c_str(), 100, 0, 100 );

}


// --- Fill histograms ---------------------------------------------------------------
void JetTreeAnalyzer::fillHistograms(int maxEntries, float minPt){

  std::cout << "Filling histograms..." << std::endl;
  
  if (fChain==0){
    std::cout<<"Error: cannot open " << fChain->GetName() << std::endl;
    exit(0);
  }

  if (maxEntries == -1)
    maxEntries = fChain->GetEntries();

  
  int nj = 0;

  for (int entry = 0; entry < maxEntries; entry++){
    fChain->GetEntry(entry);
    
    if (entry%100==0) std::cout << "Analyzing entry : " << entry << "\r" << std::flush;
        
    // --- Loop over jets in this event                                                                                                                                                                       
    for (int j = 0; j < njetsUncorr; j++){
      
      float pt = jetUncorrPt->at(j);
      
      if (pt < minPt)  continue;

      nj++;
      int imatch = (jetGenMatchIndex)->at(j);

      hrawpt->Fill(jetUncorrPt->at(j));
      hpt->Fill(jetPt->at(j));
      heta->Fill(jetEta->at(j));
      hmass->Fill(jetMass->at(j));
      hnparticles->Fill(jetNparticles->at(j));
      
      if (j == 0){
	hrawpt_leadjet->Fill(jetUncorrPt->at(j));
	hpt_leadjet->Fill(jetPt->at(j));
	heta_leadjet->Fill(jetEta->at(j));
	hmass_leadjet->Fill(jetMass->at(j));
      }

      if (imatch == -1) {
	hrawpt_pu->Fill(jetUncorrPt->at(j));
	hpt_pu->Fill(jetPt->at(j));
	heta_pu->Fill(jetEta->at(j));
	if (j ==0){
	  hrawpt_pu_leadjet->Fill(jetUncorrPt->at(j));
	  hpt_pu_leadjet->Fill(jetPt->at(j));
	  heta_pu_leadjet->Fill(jetEta->at(j));
	}
      }
      else {
	hrawpt_good->Fill(jetUncorrPt->at(j));
	hpt_good->Fill(jetPt->at(j));
	heta_good->Fill(jetEta->at(j));
	if (j == 0){
	  hrawpt_good_leadjet->Fill(jetUncorrPt->at(j));
	  hpt_good_leadjet->Fill(jetPt->at(j));
	  heta_good_leadjet->Fill(jetEta->at(j));
	}
      }

      // -- response plots
      if (imatch > -1){
	hrawpt_response->Fill(jetUncorrPt->at(j)-jetGenPt->at(j));
	hpt_response->Fill(jetPt->at(j)-jetGenPt->at(j));
	hmass_response->Fill(jetMass->at(j)-jetGenMass->at(j));
	if (j == 0){
	  hrawpt_response_leadjet->Fill(jetUncorrPt->at(j)-jetGenPt->at(j));
	  hpt_response_leadjet->Fill(jetPt->at(j)-jetGenPt->at(j));
	  hmass_response_leadjet->Fill(jetMass->at(j)-jetGenMass->at(j));
	}
      }
 
    }// end loop over jets 

    hnjets->Fill(njets);

  }// end loop over entries
}


// --- Save histograms ---------------------------------------------------------------
void JetTreeAnalyzer::saveHistograms(TFile *outfile, std::string dir){

  std::cout << "Saving histograms ... " << std::endl;
  
  outfile->cd();
  TDirectory *thisdir = outfile->mkdir(dir.c_str());
  thisdir->cd();    // make the "thisdir" directory
  
  hnjets->Write();
  hrawpt->Write();
  hrawpt_pu->Write();
  hrawpt_good->Write();
  hrawpt_response->Write();
  hpt->Write();
  hpt_pu->Write();
  hpt_good->Write();
  hpt_response->Write();

  hrawpt_leadjet->Write();
  hrawpt_pu_leadjet->Write();
  hrawpt_good_leadjet->Write();
  hrawpt_response_leadjet->Write();
  hpt_leadjet->Write();
  hpt_pu_leadjet->Write();
  hpt_good_leadjet->Write();
  hpt_response_leadjet->Write();

  heta->Write();
  heta_pu->Write();
  heta_good->Write();

  heta_leadjet->Write();
  heta_pu_leadjet->Write();
  heta_good_leadjet->Write();

  hmass->Write();
  hmass_response->Write();

  hmass_leadjet->Write();
  hmass_response_leadjet->Write();

  hnparticles->Write();
}
