#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "TClonesArray.h"
#include "Math/ProbFunc.h"

#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include <fastjet/contrib/JetCleanser.hh>

#include "puppiCleanContainer.hh"
#include "NoTrees.hh"
#include "RecoObj.hh"

#include "Bacon/BaconAna/DataFormatsOffline/interface/TGenParticle.hh"
#include "Bacon/BaconAna/DataFormatsOffline/interface/TPFPart.hh"
#include "Bacon/BaconAna/DataFormatsOffline/interface/TJet.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace baconhep;

//-------------------------------------------------------------------
// input Tree variables
//-------------------------------------------------------------------
TClonesArray *fPFPart;
TClonesArray *fGenPart;
TClonesArray *fJet;
TBranch      *fPFPartBr;
TBranch      *fGenPartBr;
TBranch      *fJetBr;

//-------------------------------------------------------------------
// output Tree variables
//-------------------------------------------------------------------
int njets_;
int njets_corr_;

float sumEt_;
float etMissX_;
float etMissY_;

// uncorrected jet variables
std::vector<float> v_jet_m_;
std::vector<float> v_jet_pt_;
std::vector<float> v_jet_eta_;
std::vector<float> v_jet_phi_;

// trimmed info
std::vector<float> v_jet_m_trimmed_;

// 4-vector subtracted
std::vector<float> v_jet_m_4Vcorr_;
std::vector<float> v_jet_pt_4Vcorr_;
std::vector<float> v_jet_eta_4Vcorr_;
std::vector<float> v_jet_phi_4Vcorr_;

// gen level quantities
std::vector<float> v_jet_genm_;
std::vector<float> v_jet_genpt_;
std::vector<float> v_jet_geneta_;
std::vector<float> v_jet_genphi_;

// jet prop.
std::vector<float> v_jet_nParticles_;

// gen matching (index of the gen jet matched to reco)
std::vector<int> v_jet_igenmatch_;

// some randoms
double curRho;
double pfRho;
double pfchsRho;
bool isPFCHS;

//-------------------------------------------------------------------
// Helper functions
//-------------------------------------------------------------------
void initVars();
void addBranches(TTree &tree);
bool FillChain(TChain& chain, const std::string& inputFileList);
void                              setupCMSReadOut(TTree *iTree );
void                              readGenCMSEvent(TTree *iTree, int entry, std::vector< fastjet::PseudoJet > &allParticles );
void                              readCMSEvent   (TTree *iTree, int entry, std::vector< RecoObj >            &allParticles ,bool iUseDeltaZ);
void                              setupCMSSWJetReadOut(TTree *iTree, float R );
void                              readCMSSWJet   (TTree *iTree, int entry, TTree &tree, std::vector<fastjet::PseudoJet> genJets);
std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, double vRparam, bool doGenMatching, std::vector < fastjet::PseudoJet > genJets );
int                               getGenMatchIndex(fastjet::PseudoJet recoJet, std::vector<fastjet::PseudoJet> genJets);
void                              setTDRStyle();

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// begin main
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

int main( int argc, char **argv ) {

    gROOT->ProcessLine("#include <vector>");
    int nEvts = 0;
    
    int useDeltaZ  = true; 

    bool doGenJets     = true;
    bool doPFJets      = true;
    bool doPFCHSJets   = true;
    bool doPuppiJets   = true;
    bool doPFCMSSWJets = true;

    // --- args
    if (argc < 5) {
      cout << "Usage: ./MySimplePuppiAnalyzer <maxEvents> <inputFilesList> <output> <radius>" <<endl;
      exit(0);
    }
    int maxEvents  = atoi(argv[1]);
    char *inputFilesList = argv[2];
    char *olabel   = argv[3];
    float jetR = atof(argv[4]);

    // --- Read list of files to be analyzed and fill TChain
    TChain* fTree = new TChain("Events");
    FillChain(*fTree, inputFilesList);

    if (maxEvents < 0)
      maxEvents = fTree->GetEntries();

    cout << "This analysis will run on "<< maxEvents << " entries" <<endl; 


    // --- Setup input tree readout
    setupCMSReadOut(fTree);
    setupCMSSWJetReadOut(fTree, jetR);
    
    // ---  Output trees
    char tfname[150];
    sprintf( tfname, "outtree_%s.root", olabel );
    TFile fout_(tfname, "RECREATE");

    TTree* tree_gen;
    TTree* tree_pf;
    TTree* tree_pfchs;
    TTree* tree_pf_puppi;
    TTree* tree_pf_cmssw; 
    
    if (doGenJets){
      tree_gen = new TTree("tree_gen", "tree_gen");
      addBranches(*tree_gen);
    }
    
    if (doPFJets){
      tree_pf = new TTree("tree_pf", "tree_pf");
      addBranches(*tree_pf);
    }
    
    if (doPFCHSJets){
      tree_pfchs = new TTree("tree_pfchs", "tree_pfchs");
      addBranches(*tree_pfchs);
    }

    if (doPuppiJets){
      tree_pf_puppi = new TTree("tree_pf_puppi", "tree_pf_puppi");
      addBranches(*tree_pf_puppi);  
    }    

    if (doPFCMSSWJets){
      tree_pf_cmssw = new TTree("tree_pf_cmssw", "tree_pf_cmssw"); // tree for jets from CMSSW objects
      addBranches(*tree_pf_cmssw);
    }

    std::vector < RecoObj > allParticles;
    std::vector < fastjet::PseudoJet > genParticles;
    std::vector < fastjet::PseudoJet > pfParticles;
    std::vector < fastjet::PseudoJet > chsParticles;
    std::vector < fastjet::PseudoJet > puppiParticles;

    allParticles     .clear();
    genParticles     .clear();
    pfParticles      .clear();
    chsParticles     .clear();
    puppiParticles   .clear();


    // --- loop over input tree
    for (int entry = 0; entry < maxEvents ;  entry++){

      if (entry%10 == 0 ) std::cout << "Analysis is " << ((float) entry)*100./maxEvents << "% done" << std::endl;   

      // --- read event
      readCMSEvent(fTree, entry, allParticles, useDeltaZ);
      readGenCMSEvent(fTree, entry, genParticles);
      
      // --- buil particles collections
      puppiCleanContainer curEvent(allParticles);
      puppiParticles      = curEvent.puppiEvent(8,0.5); // puppi metric: 7 = log, 8 = log^2
      pfParticles         = curEvent.pfParticles();
      chsParticles        = curEvent.pfchsParticles();
              
      curRho   = -1;
      pfRho    = -1;
      pfchsRho = -1;
      
      // --- gen jets
      initVars();
      std::vector<fastjet::PseudoJet> dummy;
      std::vector<fastjet::PseudoJet> genJets;
      if (doGenJets) genJets = analyzeEvent( genParticles, *tree_gen, jetR, false, dummy);

      // --- PF jets
      initVars();
      std::vector<fastjet::PseudoJet> pfJets;
      if (doPFJets) pfJets = analyzeEvent( pfParticles, *tree_pf, jetR, true, genJets);
      pfRho = curRho;

      // --- PFCHS jets
      initVars();
      isPFCHS = true;
      std::vector<fastjet::PseudoJet> pfchsJets;
      if (doPFCHSJets) pfchsJets = analyzeEvent( chsParticles, *tree_pfchs, jetR, true, genJets);
      pfchsRho = curRho;
      isPFCHS = false;

      // --- PUPPI jets
      initVars();
      std::vector<fastjet::PseudoJet> puppiJets; 
      if (doPuppiJets) puppiJets = analyzeEvent( puppiParticles, *tree_pf_puppi, jetR, true, genJets);

      // --- read the tree containing standard CMSSW jets and fill directly the output tree   
      initVars();
      if (doPFCMSSWJets) readCMSSWJet(fTree, entry, *tree_pf_cmssw, genJets);
                
      // --- clear all
      genParticles.clear();
      allParticles.clear();
      pfParticles .clear();
      chsParticles.clear();
      puppiParticles.clear();
      
    }
    
    // write output trees
    fout_.cd();
    if (doGenJets) tree_gen->Write();
    if (doPFJets) tree_pf->Write();
    if (doPFCHSJets) tree_pfchs->Write();
    if (doPuppiJets) tree_pf_puppi->Write();
    if (doPFCMSSWJets) tree_pf_cmssw->Write();
    fout_.Close();
    
    return 0;
}
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// end main
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------




// ------------------------------------------------------------------------------------------
void initVars(){
    njets_ = -1.;
    njets_corr_ = -1.;
    
    sumEt_ = -1.;
    etMissX_ = -1.;
    etMissY_ = -1.;
    
    v_jet_m_.clear();
    v_jet_pt_.clear();
    v_jet_eta_.clear();
    v_jet_phi_.clear();
    v_jet_m_trimmed_.clear();
    v_jet_m_4Vcorr_.clear();
    v_jet_pt_4Vcorr_.clear();
    v_jet_eta_4Vcorr_.clear();
    v_jet_phi_4Vcorr_.clear();
    v_jet_nParticles_.clear();
    v_jet_igenmatch_.clear();
    v_jet_genm_.clear();
    v_jet_genpt_.clear();
    v_jet_geneta_.clear();
    v_jet_genphi_.clear();
}
// ------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------
void addBranches( TTree &tree ){
    
    tree.Branch("njetsUncorr",&njets_);
    tree.Branch("njets",&njets_corr_);
    
    tree.Branch("sumEt",&sumEt_);
    tree.Branch("etMissX",&etMissX_);
    tree.Branch("etMissY",&etMissY_);
    
    //tree.Branch("v_jet_m",&v_jet_m_);
    //tree.Branch("v_jet_pt",&v_jet_pt_);
    //tree.Branch("v_jet_eta",&v_jet_eta_);
    //tree.Branch("v_jet_phi",&v_jet_phi_);
    
    //tree.Branch("v_jet_m_trimmed",&v_jet_m_trimmed_);
    
    //tree.Branch("v_jet_m_4Vcorr",&v_jet_m_4Vcorr_);
    //tree.Branch("v_jet_pt_4Vcorr",&v_jet_pt_4Vcorr_);
    //tree.Branch("v_jet_eta_4Vcorr",&v_jet_eta_4Vcorr_);
    //tree.Branch("v_jet_phi_4Vcorr",&v_jet_phi_4Vcorr_);

    //tree.Branch("v_jet_igenmatch",&v_jet_igenmatch_);

    tree.Branch("jetUncorrMass",&v_jet_m_);
    tree.Branch("jetUncorrPt",&v_jet_pt_);
    tree.Branch("jetUncorrEta",&v_jet_eta_);
    tree.Branch("jetUncorrPhi",&v_jet_phi_);
    
    tree.Branch("jetTrimmedMass",&v_jet_m_trimmed_);
    
    tree.Branch("jetMass",&v_jet_m_4Vcorr_);
    tree.Branch("jetPt",&v_jet_pt_4Vcorr_);
    tree.Branch("jetEta",&v_jet_eta_4Vcorr_);
    tree.Branch("jetPhi",&v_jet_phi_4Vcorr_);

    tree.Branch("jetGenMass",&v_jet_genm_);
    tree.Branch("jetGenPt",&v_jet_genpt_);
    tree.Branch("jetGenEta",&v_jet_geneta_);
    tree.Branch("jetGenPhi",&v_jet_genphi_);

    tree.Branch("jetNparticles",&v_jet_nParticles_);

    tree.Branch("jetGenMatchIndex",&v_jet_igenmatch_);
    
}
// ------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------
bool FillChain(TChain& chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open())
    {
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
    }
  
  while(1)
    {
      inFile >> buffer;
      if(!inFile.good()) break;
      chain.Add(buffer.c_str());
    }

  return true;
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void setupCMSReadOut(TTree *iTree ) {
    
    fPFPart  = new TClonesArray("baconhep::TPFPart");
    fGenPart = new TClonesArray("baconhep::TGenParticle");
    iTree->SetBranchAddress("PFPart",       &fPFPart);
    iTree->SetBranchAddress("GenParticle",  &fGenPart);
    fPFPartBr  = iTree->GetBranch("PFPart");
    fGenPartBr = iTree->GetBranch("GenParticle");
    
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void setupCMSSWJetReadOut(TTree *iTree, float R ) {
  
    cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
    fJet  = new TClonesArray("baconhep::TJet");
    iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
    fJetBr  = iTree->GetBranch(Form("Jet0%d",R*10));

}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void readCMSSWJet(TTree *iTree, int entry, TTree &tree, std::vector<fastjet::PseudoJet> genJets) {

  iTree->GetEntry(entry);

  for (int i = 0; i < fJet->GetEntriesFast(); i++){
    TJet *pJet = (TJet*)((*fJet)[i]);
    
    v_jet_pt_.push_back( pJet->ptRaw );

    //use v_jet_*_4Vcorr_ to fill corrected jets even if they have ALL corrections (L1L2L3)
    v_jet_m_4Vcorr_.push_back(pJet->mass);
    v_jet_pt_4Vcorr_.push_back(pJet->pt);
    v_jet_eta_4Vcorr_.push_back(pJet->eta);
    v_jet_phi_4Vcorr_.push_back(pJet->phi);
    v_jet_nParticles_.push_back(pJet->nParticles);

    // gen matching
    int imatch = -1;
    double mindr = 0.25;
    TLorentzVector *recojet = new TLorentzVector();
    recojet->SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    for (int ig = 0; ig < genJets.size(); ig++){
      TLorentzVector *genjet = new TLorentzVector();
      genjet->SetPtEtaPhiM(genJets[ig].pt(), genJets[ig].eta(), genJets[ig].phi(), genJets[ig].m());
      double dr = recojet->DeltaR(*genjet);
      if (dr < mindr){
	mindr = dr;
	imatch = ig;
      }
      delete genjet;
    }

    delete recojet;

    v_jet_igenmatch_.push_back(imatch);
    
    njets_++;
    if (pJet->pt > 25){
      njets_corr_++;
    }
  }

  // event quantities                                                
  tree.Fill();
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++                          
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void readCMSEvent(TTree *iTree, int entry, std::vector< RecoObj > &allParticles,bool iUseDeltaZ) {
    float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
    iTree->GetEntry(entry);
    for (int i = 0; i < fPFPart->GetEntriesFast(); i++){
        TPFPart *pPart = (TPFPart*)((*fPFPart)[i]);
	//Standard PF No PU Association
	//Charged particles are counted if they have a vertex assignment
        //PV particles are all PV chargd particles + all unassociated particles

        //isCh   = (fabs(pPart->q) > 0) && (pPart->vtxId > -1);
	isCh   = (fabs(pPart->q) > 0) && (pPart->vtxId >= -1); // this allows to use vtxId =0 || -1 for ChLV 
        isPV   = (pPart->vtxId <= 0);
        if(iUseDeltaZ) isCh   = (pPart->pfType == 1 || pPart->pfType == 2 || pPart->pfType == 3) && (pPart->vtxId > -1 || fabs(pPart->dz) < 0.2) ;
	if(iUseDeltaZ) isPV   = (pPart->vtxId  == 0  || (fabs(pPart->dz) < 0.2 && isCh));
	
        TLorentzVector pVec; pVec.SetPtEtaPhiM(pPart->pt,pPart->eta,pPart->phi,pPart->m);
        // Id: 0 = NeLV, 1 = NePU, 2 = ChLV, 3 = ChPU
        int lID = -1;
        if (!isCh) lID = 1;
        if (isCh && isPV) lID = 2;
        if (isCh && !isPV) lID = 3;

        RecoObj curPseudoJet;
        curPseudoJet.pt  = pPart->pt;
        curPseudoJet.eta = pPart->eta;
        curPseudoJet.phi = pPart->phi;
        curPseudoJet.m   = pPart->m;
        curPseudoJet.id  = lID;
        curPseudoJet.vtxId = pPart->vtxId;
        curPseudoJet.trkChi2 = pPart->trkChi2;
        curPseudoJet.vtxChi2 = pPart->vtxChi2;
        curPseudoJet.pfType  = pPart->pfType;
        curPseudoJet.depth   = pPart->depth;
        curPseudoJet.time    = pPart->time;
        curPseudoJet.d0        = pPart->d0;
        curPseudoJet.dZ        = pPart->dz;
        allParticles.push_back( curPseudoJet );
    }
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void readGenCMSEvent(TTree *iTree, int entry, std::vector< fastjet::PseudoJet > &allParticles) {
    float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
    iTree->GetEntry(entry);
    for (int i = 0; i < fGenPart->GetEntriesFast(); i++){
        
        //        if(fPartPt == -1) break;
        TGenParticle *pPart = (TGenParticle*)((*fGenPart)[i]);
            
        if(pPart->status != 1) continue;
        if(fabs(pPart->pdgId) == 12 ||
           fabs(pPart->pdgId) == 14 ||
           fabs(pPart->pdgId) == 16) continue;
        
        TLorentzVector pVec;
        pVec.SetPtEtaPhiM(pPart->pt,pPart->eta,pPart->phi,pPart->mass);
        
        fastjet::PseudoJet curPseudoJet (pVec.Px(), pVec.Py(), pVec.Pz(), pVec.E());
        curPseudoJet.set_user_index(-1);
        
        if (curPseudoJet.pt() > 0) allParticles.push_back( curPseudoJet );
    }
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, 
						TTree &tree, double vRparam, bool doGenMatching, 
						std::vector < fastjet::PseudoJet > genJets){
    // recluster on the fly....
    double rParam = vRparam;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;

    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    fastjet::ClusterSequenceArea* thisClustering_ = new fastjet::ClusterSequenceArea(constits, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets_ = sorted_by_pt(thisClustering_->inclusive_jets(25.0));

    // trim jet
    fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05)));

    // 4-vector subtraction
    fastjet::GridMedianBackgroundEstimator lGrid(5.0,0.8);
    lGrid.set_particles(constits);
    
    curRho = lGrid.rho();
    double forwardRho = lGrid.rho();
    double centralRho = lGrid.rho();
    if (isPFCHS){ forwardRho = pfRho; }
    //std::cout << "centralRho = " << centralRho << ", forwardRho = " << forwardRho << std::endl;
    
    njets_corr_ = 0;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    std::vector<fastjet::PseudoJet> corrjets;
    njets_ = (int) out_jets_.size();
    for (unsigned int i = 0; i < out_jets_.size(); i++){
        
        v_jet_m_.push_back( out_jets_[i].m() );
        v_jet_pt_.push_back( out_jets_[i].pt() );
        v_jet_eta_.push_back( out_jets_[i].eta() );
        v_jet_phi_.push_back( out_jets_[i].phi() );
        
        // trim jet
        fastjet::PseudoJet trimmedJet = (trimmer)(out_jets_.at(i));
        v_jet_m_trimmed_.push_back( trimmedJet.m() );
        
        // 4-vector subtraction
        PseudoJet pCorrJet = out_jets_.at(i);
        PseudoJet pArea = out_jets_.at(i).area_4vector();
        
        double tmpRho = -9999;
        if (fabs(pCorrJet.eta()) > 2.5) tmpRho = forwardRho;
        else tmpRho = centralRho;
        
        pCorrJet -= tmpRho * pArea;
        
        v_jet_m_4Vcorr_.push_back(pCorrJet.m());
        v_jet_pt_4Vcorr_.push_back(pCorrJet.pt());
        v_jet_eta_4Vcorr_.push_back(pCorrJet.eta());
        v_jet_phi_4Vcorr_.push_back(pCorrJet.phi());
        
        if (v_jet_pt_4Vcorr_[i] > 25.0){
            njets_corr_++;
            //corrjets.push_back(pCorrJet);
        }

	corrjets.push_back(pCorrJet); //don't cut on pt
	
	// gen matching
	int imatch = -1;
	if (doGenMatching && !(genJets.empty())){
	  imatch = getGenMatchIndex(pCorrJet,genJets);
	}
	v_jet_igenmatch_.push_back(imatch);
	
	if (imatch > -1){
	  v_jet_genm_.push_back( genJets[imatch].m() );
	  v_jet_genpt_.push_back( genJets[imatch].pt() );
	  v_jet_geneta_.push_back( genJets[imatch].eta() );
	  v_jet_genphi_.push_back( genJets[imatch].phi() );
	}
	else {
	  v_jet_genm_.push_back( -999. );
          v_jet_genpt_.push_back( -999. );
          v_jet_geneta_.push_back( -999. );
          v_jet_genphi_.push_back( -999.);
	}

	v_jet_nParticles_.push_back((out_jets_[i].constituents()).size());
	
    }
    
    // MET variables
    float sumEt = 0;
    PseudoJet METvec(0,0,0,0);
    for (unsigned int i = 0; i < constits.size(); i++){
        sumEt += constits[i].Et();
        METvec += constits[i];
    }
    
    sumEt_ = sumEt;
    etMissX_ = fabs(-METvec.px());
    etMissY_ = fabs(-METvec.py());
    
    // event quantities
    tree.Fill();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if (out_jets_.size() > 0) thisClustering_->delete_self_when_unused();
    
    return corrjets;
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
int getGenMatchIndex(fastjet::PseudoJet recoJet, std::vector<fastjet::PseudoJet> genJets){
  
  //  for (int ir = 0; ir < recoJets.size(); ir++){
  int imatch = -1;
  double mindr = 0.25;
  for (int ig = 0; ig < genJets.size(); ig++){
    //double dr = recoJets[ir].delta_R(genJets[ig]);
    double dr = recoJet.delta_R(genJets[ig]);
    if (dr < mindr){
      mindr = dr;
      imatch = ig;
    }
  }
  return imatch;
}
// ------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------
void setTDRStyle() {
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
    
    // For the canvas:
    tdrStyle->SetCanvasBorderMode(0);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(750); //Height of canvas
    tdrStyle->SetCanvasDefW(1050); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);
    
    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);
    
    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);
    
    // For the histo:
    // tdrStyle->SetHistFillColor(1);
    // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
    // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
    // tdrStyle->SetNumberContours(Int_t number = 20);
    
    tdrStyle->SetEndErrorSize(2);
    //  tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(0.);
    
    tdrStyle->SetMarkerStyle(20);
    
    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);
    
    //For the date:
    tdrStyle->SetOptDate(0);
    // tdrStyle->SetDateX(Float_t x = 0.01);
    // tdrStyle->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    tdrStyle->SetOptFile(0);
    tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.010);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.25);
    tdrStyle->SetStatW(0.15);
    // tdrStyle->SetStatStyle(Style_t style = 1001);
    // tdrStyle->SetStatX(Float_t x = 0);
    // tdrStyle->SetStatY(Float_t y = 0);
    
    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.13);
    tdrStyle->SetPadLeftMargin(0.17);
    tdrStyle->SetPadRightMargin(0.17);
    
    // For the Global title:
    
    tdrStyle->SetOptTitle(0);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.005);
    // tdrStyle->SetTitleH(0); // Set the height of the title box
    // tdrStyle->SetTitleW(0); // Set the width of the title box
    // tdrStyle->SetTitleX(0); // Set the position of the title box
    // tdrStyle->SetTitleY(0.985); // Set the position of the title box
    // tdrStyle->SetTitleStyle(Style_t style = 1001);
    // tdrStyle->SetTitleBorderSize(2);
    
    // For the axis titles:
    
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleSize(0.06, "XYZ");
    // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // tdrStyle->SetTitleYSize(Float_t size = 0.02);
    tdrStyle->SetTitleXOffset(0.9);
    tdrStyle->SetTitleYOffset(1.25);
    // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XYZ");
    
    // For the axis:
    
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(505, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);
    
    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);
    
    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
    // tdrStyle->SetLineScalePS(Float_t scale = 3);
    // tdrStyle->SetLineStyleString(Int_t i, const char* text);
    // tdrStyle->SetHeaderPS(const char* header);
    // tdrStyle->SetTitlePS(const char* pstitle);
    
    // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
    // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
    // tdrStyle->SetPaintTextFormat(const char* format = "g");
    // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // tdrStyle->SetTimeOffset(Double_t toffset);
    // tdrStyle->SetHistMinimumZero(kTRUE);
    
    tdrStyle->SetPalette(1);
    
    
    tdrStyle->cd();
    
}
