//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 13 18:03:09 2014 by ROOT version 5.32/00
// from TTree tree_gen/tree_gen
// found on file: outtree_test.root
//////////////////////////////////////////////////////////

#ifndef JetTree_h
#define JetTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTree {
public :
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

   JetTree(TTree *tree=0);
   virtual ~JetTree();
   virtual void  Init(TTree *tree);
   virtual Int_t GetEntry(Long64_t entry);
};
#endif

JetTree::JetTree(TTree *tree) 
{
  Init(tree);
}

JetTree::~JetTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void JetTree::Init(TTree *tree)
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

