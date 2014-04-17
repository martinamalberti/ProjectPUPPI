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
//#include "puppiTMVAContainer.hh"
#include "NoTrees.hh"
#include "RecoObj.hh"

#include "Bacon/BaconAna/DataFormatsOffline/interface/TGenParticle.hh"
#include "Bacon/BaconAna/DataFormatsOffline/interface/TPFPart.hh"
#include "Bacon/BaconAna/DataFormatsOffline/interface/TJet.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace baconhep;

ifstream fin;
/////////////////////////////////////////////////////////////////////
// iTree variables
/////////////////////////////////////////////////////////////////////

TClonesArray *fPFPart;
TClonesArray *fGenPart;
TClonesArray *fJet;
TBranch      *fPFPartBr;
TBranch      *fGenPartBr;
TBranch      *fJetBr;

/////////////////////////////////////////////////////////////////////
// oTree variables
/////////////////////////////////////////////////////////////////////
int fCount;
int fCount2;

//ifstream fin;
int njets_;
int njets_corr_;

float sumEt_;
float etMissX_;
float etMissY_;

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

// gen matching (index of the 
std::vector<int> v_jet_igenmatch_;

float p_isPU, p_isCH, p_px, p_py, p_pz, p_e, p_puppiW_pfchs, p_cleansedW, p_puppiW_chLV, p_puppiW_all;
float p_alphas_chLV, p_alphas_all;

// some randoms
double curRho;
double pfRho;
double pfchsRho;
bool isPFCHS;

/////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
void                              setupCMSReadOut(TTree *iTree );
void                              readGenCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles );
void                              readCMSEvent   (TTree *iTree, std::vector< RecoObj >            &allParticles ,bool iUseDeltaZ);
void                              setupCMSSWJetReadOut(TTree *iTree, float R );
void                              readCMSSWJet   (TTree *iTree, TTree &tree, std::vector<fastjet::PseudoJet> genJets);
std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag, double vRparam, bool doGenMatching, std::vector < fastjet::PseudoJet > genJets );
//void getGenMatchIndex(std::vector<fastjet::PseudoJet> recoJets, std::vector<fastjet::PseudoJet> genJets, std::vector<int> &indexes);
int                               getGenMatchIndex(fastjet::PseudoJet recoJet, std::vector<fastjet::PseudoJet> genJets);
void                              plotEvent( std::vector < fastjet::PseudoJet > constits, char* name, std::vector < fastjet::PseudoJet > jets );
void                              setTDRStyle();


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
    v_jet_igenmatch_.clear();
}

void addBranches( TTree &tree ){
    
    tree.Branch("njets",&njets_);
    tree.Branch("njets_corr",&njets_corr_);
    
    tree.Branch("sumEt",&sumEt_);
    tree.Branch("etMissX",&etMissX_);
    tree.Branch("etMissY",&etMissY_);
    
    tree.Branch("v_jet_m",&v_jet_m_);
    tree.Branch("v_jet_pt",&v_jet_pt_);
    tree.Branch("v_jet_eta",&v_jet_eta_);
    tree.Branch("v_jet_phi",&v_jet_phi_);
    
    tree.Branch("v_jet_m_trimmed",&v_jet_m_trimmed_);
    
    tree.Branch("v_jet_m_4Vcorr",&v_jet_m_4Vcorr_);
    tree.Branch("v_jet_pt_4Vcorr",&v_jet_pt_4Vcorr_);
    tree.Branch("v_jet_eta_4Vcorr",&v_jet_eta_4Vcorr_);
    tree.Branch("v_jet_phi_4Vcorr",&v_jet_phi_4Vcorr_);

    tree.Branch("v_jet_igenmatch",&v_jet_igenmatch_);
    
}

void initParticleVars(){
    p_isPU = -1.;
    p_isCH = -1.;
    
    p_px = -1.;
    p_py = -1.;
    p_pz = -1.;
    p_e = -1.;
    
    p_puppiW_pfchs = -1.;
    p_puppiW_chLV = -1.;
    p_puppiW_all = -1.;
    
    p_alphas_chLV = -1.;
    p_alphas_all = -1.;
    
    
    p_cleansedW = -1.;
}

void addParticleBranches( TTree &tree ){

    tree.Branch("p_isPU",&p_isPU);
    tree.Branch("p_isCH",&p_isCH);
    tree.Branch("p_px",&p_px);
    tree.Branch("p_py",&p_py);
    tree.Branch("p_pz",&p_pz);
    tree.Branch("p_e",&p_e);

    tree.Branch("p_puppiW_chLV",&p_puppiW_chLV);
    tree.Branch("p_puppiW_all",&p_puppiW_all);
    tree.Branch("p_alphas_chLV",&p_alphas_chLV);
    tree.Branch("p_alphas_all",&p_alphas_all);

}


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
// ------------------------------------------------------------------------------------------
// begin main
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

int main( int argc, char **argv ) {
    setTDRStyle();
    
    gROOT->ProcessLine("#include <vector>");
    int nEvts = 0;
    
    //int maxEvents  = atoi(argv[1]);
    //char *recoFile = argv[2];;
    //char *olabel   = argv[3];
    ////int useMVA     = atoi(argv[4]);
    //int useDeltaZ  = true; //atoi(argv[5])
    //cout << "==> " << olabel << endl;
    //TFile* lReadout = new TFile(recoFile);
    //TTree* fTree    = (TTree*) lReadout->FindObjectAny("Events");

    int useDeltaZ  = false; 

    // args
    int maxEvents  = atoi(argv[1]);
    char *inputFilesList = argv[2];
    char *olabel   = argv[3];
    float jetR = atof(argv[4]);

    cout << "maxEvents = " << maxEvents <<endl;

    // Read listof files to be analyzed and fill TChain
    TChain* fTree = new TChain("Events");
    FillChain(*fTree, inputFilesList);

    setupCMSReadOut(fTree);
    setupCMSSWJetReadOut(fTree, jetR);

    if (maxEvents < 0)
      maxEvents = fTree->GetEntries();
    
    /////////////////////////
    // otrees
    char tfname[150];
    sprintf( tfname, "outtree_%s.root", olabel );
    TFile fout_(tfname, "RECREATE");
    TTree* tree_gen = new TTree("tree_gen", "tree_gen");
    TTree* tree_pf = new TTree("tree_pf", "tree_pf");
    TTree* tree_pfchs = new TTree("tree_pfchs", "tree_pfchs");
    TTree* tree_pf_puppi = new TTree("tree_pf_puppi", "tree_pf_puppi");
    
    TTree* tree_pf_cmssw = new TTree("tree_pf_cmssw", "tree_pf_cmssw"); // tree for jets from CMSSW objects

    addBranches(*tree_gen);
    addBranches(*tree_pf);
    addBranches(*tree_pfchs);
    addBranches(*tree_pf_puppi);
    addBranches(*tree_pf_cmssw);
    
    TTree* tree_particles = new TTree("tree_particles", "tree_particles");
    addParticleBranches(*tree_particles);
    /////////////////////////
    
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

    while(true){
        
        nEvts++;
	if (nEvts % 10 == 0) std::cout << "file is " << ((float) nEvts)*100./maxEvents << "% done" << std::endl;
        

        if (nEvts == maxEvents){ break; }
        readCMSEvent(fTree, allParticles,useDeltaZ);
        readGenCMSEvent(fTree, genParticles);
	
        puppiCleanContainer curEvent(allParticles);
        puppiParticles      = curEvent.puppiEvent(7,0.5); // puppi metric: 7 = log, 8 = log^2
        pfParticles         = curEvent.pfParticles();
        chsParticles        = curEvent.pfchsParticles();
        
        char canvname[150];
        
        curRho   = -1;
        pfRho    = -1;
        pfchsRho = -1;
        //std::cout << "--------------------" << std::endl;
        initVars();
        //std::cout << "analyze gen" << std::endl;
	std::vector<fastjet::PseudoJet> dummy;
        std::vector<fastjet::PseudoJet> genJets = analyzeEvent( genParticles, *tree_gen, canvname, jetR, false, dummy);
	initVars();
        //std::cout << "analyze pf" << std::endl;
        std::vector<fastjet::PseudoJet> pfJets = analyzeEvent( pfParticles, *tree_pf, canvname, jetR, true, genJets);
	pfRho = curRho;
        initVars();
        //std::cout << "analyze pfchs" << std::endl;
        isPFCHS = true;
	std::vector<fastjet::PseudoJet> pfchsJets = analyzeEvent( chsParticles, *tree_pfchs, canvname, jetR, true, genJets);
	pfchsRho = curRho;
        isPFCHS = false;
        initVars();
        //std::cout << "analyze puppi" << std::endl;
        std::vector<fastjet::PseudoJet> puppiJets = analyzeEvent( puppiParticles, *tree_pf_puppi, canvname, jetR, true, genJets);
	//
	
	// read the tree containing standard CMSSW jets and fill directly the output tree   
	initVars();
	readCMSSWJet(fTree,*tree_pf_cmssw, genJets);

        
        if (nEvts < -100 and nEvts > 0){
            
            std::vector<float> puppiWeights_chLV = curEvent.getPuppiWeights_chLV();
            std::vector<float> puppiWeights_all = curEvent.getPuppiWeights_all();
            std::vector<float> alphas_chLV = curEvent.getPuppiAlphas_chLV();
            std::vector<float> alphas_all = curEvent.getPuppiAlphas_all();

            // fill weights
            for (unsigned int a = 0; a < pfParticles.size(); a++){
                if (pfParticles[a].user_index() == 1 || pfParticles[a].user_index() == 3) p_isPU = 1;
                else p_isPU = 0;
                if (pfParticles[a].user_index() == 2 || pfParticles[a].user_index() == 3) p_isCH = 1;
                else p_isCH = 0;
                p_px = pfParticles[a].px();
                p_py = pfParticles[a].py();
                p_pz = pfParticles[a].pz();
                p_e = pfParticles[a].e();

                p_puppiW_chLV = puppiWeights_chLV[a];
                p_puppiW_all = puppiWeights_all[a];
                p_alphas_chLV = alphas_chLV[a];
                p_alphas_all = alphas_all[a];
                tree_particles->Fill();
            }
            
            sprintf( canvname, "plotting/displays/dis_gen_%i_%s", nEvts, olabel );
            plotEvent( genParticles, canvname, genJets );
            sprintf( canvname, "plotting/displays/dis_pf_%i_%s", nEvts, olabel );
            plotEvent( pfParticles, canvname, pfJets);
            sprintf( canvname, "plotting/displays/dis_pfchs_%i_%s", nEvts, olabel );
            plotEvent( chsParticles, canvname, pfchsJets );
            sprintf( canvname, "plotting/displays/dis_pf_puppi_%i_%s", nEvts, olabel );
            plotEvent( puppiParticles, canvname, puppiJets );
            
        }
        
        genParticles.clear();
        allParticles.clear();
        pfParticles .clear();
        chsParticles.clear();
        puppiParticles.clear();
        
        fCount++;
    }
    
    fout_.cd();
    tree_gen->Write();
    tree_pf->Write();
    tree_pfchs->Write();
    tree_pf_puppi->Write();
    tree_pf_cmssw->Write();
    tree_particles->Write();
    fout_.Close();
    
    return 0;
}

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// end main
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

void setupCMSReadOut(TTree *iTree ) {

    fCount = 0;
    
    fPFPart  = new TClonesArray("baconhep::TPFPart");
    fGenPart = new TClonesArray("baconhep::TGenParticle");
    iTree->SetBranchAddress("PFPart",       &fPFPart);
    iTree->SetBranchAddress("GenParticle",  &fGenPart);
    fPFPartBr  = iTree->GetBranch("PFPart");
    fGenPartBr = iTree->GetBranch("GenParticle");
    
}


void setupCMSSWJetReadOut(TTree *iTree, float R ) {
  
    cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
    fJet  = new TClonesArray("baconhep::TJet");
    iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
    fJetBr  = iTree->GetBranch(Form("Jet0%d",R*10));

    
}


void readCMSSWJet(TTree *iTree, TTree &tree, std::vector<fastjet::PseudoJet> genJets) {

  iTree->GetEntry(fCount);

  for (int i = 0; i < fJet->GetEntriesFast(); i++){
    TJet *pJet = (TJet*)((*fJet)[i]);
    
    v_jet_pt_.push_back( pJet->ptRaw );

    //use v_jet_*_4Vcorr_ to fill corrected jets even if they have ALL corrections (L1L2L3)
    v_jet_m_4Vcorr_.push_back(pJet->mass);
    v_jet_pt_4Vcorr_.push_back(pJet->pt);
    v_jet_eta_4Vcorr_.push_back(pJet->eta);
    v_jet_phi_4Vcorr_.push_back(pJet->phi);

    // gen matching
    int imatch = -1;
    double mindr = 0.4;
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



void readCMSEvent(TTree *iTree, std::vector< RecoObj > &allParticles,bool iUseDeltaZ) {
    float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
    iTree->GetEntry(fCount);
    for (int i = 0; i < fPFPart->GetEntriesFast(); i++){
        TPFPart *pPart = (TPFPart*)((*fPFPart)[i]);
	//Standard PF No PU Association
	//Charged particles are counted if they have a vertex assignment
        // PV particles are all PV chargd particles + all unassociated particles

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
    //fCount++;
}

void readGenCMSEvent(TTree *iTree, std::vector< fastjet::PseudoJet > &allParticles) {
    float px, py, pz, e, pdgid, isCh, isPU,isPV = 0;
    iTree->GetEntry(fCount);
//    std::cout << "fGenPart->GetEntriesFast() = " << fGenPart->GetEntriesFast() << std::endl;
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
    //fCount++;
}

void plotEvent( std::vector < fastjet::PseudoJet > constits, char* name, std::vector < fastjet::PseudoJet > jets ){
    
    double maxCell = 0;
    
    TH2F* h2d = new TH2F( "h2d",";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    TH2F* h2d_PU = new TH2F( "h2d_PU",";#eta;#phi;pT (GeV)",50, -5,5, 50,0,2*TMath::Pi() );
    
    for (unsigned int i = 0; i < constits.size(); i++){
        int curBin = h2d->FindBin(constits[i].eta(),constits[i].phi());
        if (constits[i].pt() > maxCell) maxCell = constits[i].pt();
        h2d->SetBinContent( curBin, h2d->GetBinContent(curBin) + constits[i].pt() );
    }
    
    TCanvas* can = new TCanvas("can","can",1100,800);
    
    h2d->SetMaximum( 1.1*max(h2d->GetMaximum(),h2d_PU->GetMaximum()) );
    h2d->SetMinimum( 1e-1 );
    //    h2d->SetLineWidth( 2 );
    h2d->Draw("COLZ");
    h2d_PU->SetLineWidth( 1 );
    h2d_PU->SetLineColor( 14 );
    h2d_PU->Draw("BOX");
    h2d->Draw("COLZ SAME");
    can->SetLogz();
    
    //draw jets
    for (unsigned j = 0; j < jets.size(); j++){
        TEllipse* cir = new TEllipse(jets[j].eta(),jets[j].phi(),0.7,0.7);
        cir->SetFillStyle(0);
        if (jets[j].pt() > 50 && jets[j].pt() < 200){
            cir->SetLineColor( 7 );
            cir->SetLineWidth( 2 );
        }
        else if (jets[j].pt() > 200){
            cir->SetLineColor( 4 );
            cir->SetLineWidth( 2 );
        }
        else{
            cir->SetLineColor( 6 );
            cir->SetLineWidth( 2 );
        }
        
        cir->Draw("sames");
    }
    
    char oname[96];
    sprintf(oname,"%s.pdf",name);
    can->SaveAs(oname);
    sprintf(oname,"%s.png",name);
    can->SaveAs(oname);
    delete h2d;
    delete h2d_PU;
    delete can;
}

std::vector< fastjet::PseudoJet > analyzeEvent( std::vector < fastjet::PseudoJet > constits, TTree &tree, char* tag, double vRparam, bool doGenMatching, std::vector < fastjet::PseudoJet > genJets){
    
    //std::cout << "constits.size() = " << constits.size() << std::endl;
    
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
    //std::cout << "out_jets_.size() = " << out_jets_.size() << std::endl;
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


//void getGenMatchIndex(std::vector<fastjet::PseudoJet> recoJets, std::vector<fastjet::PseudoJet> genJets, std::vector<int> &indexes){
int getGenMatchIndex(fastjet::PseudoJet recoJet, std::vector<fastjet::PseudoJet> genJets){
  
  //  for (int ir = 0; ir < recoJets.size(); ir++){
  int imatch = -1;
  double mindr = 0.4;
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


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

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
