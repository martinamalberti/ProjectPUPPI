#include "./JetTreeAnalyzer.h"


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

  TTree *tree_gen   = (TTree *)inputFile->Get("tree_gen");
  TTree *tree_pf    = (TTree *)inputFile->Get("tree_pf");
  TTree *tree_pfchs = (TTree *)inputFile->Get("tree_pfchs");
  TTree *tree_puppi = (TTree *)inputFile->Get("tree_pf_puppi");

  int maxEntries = -1;

  JetTreeAnalyzer *genAnalyzer = new JetTreeAnalyzer(tree_gen);
  genAnalyzer->bookHistograms("gen");
  genAnalyzer->fillHistograms(maxEntries,25.);


  JetTreeAnalyzer *pfAnalyzer = new JetTreeAnalyzer(tree_pf);
  pfAnalyzer->bookHistograms("pf");
  pfAnalyzer->fillHistograms(maxEntries,25.);


  JetTreeAnalyzer *pfchsAnalyzer = new JetTreeAnalyzer(tree_pfchs);
  pfchsAnalyzer->bookHistograms("pfchs");
  pfchsAnalyzer->fillHistograms(maxEntries,25.);


  JetTreeAnalyzer *puppiAnalyzer = new JetTreeAnalyzer(tree_puppi);
  puppiAnalyzer->bookHistograms("puppi");
  puppiAnalyzer->fillHistograms(maxEntries,25.);
  

  // save results in file
  TFile *outfile = new TFile("out.root","RECREATE");
  genAnalyzer->saveHistograms(outfile);
  pfAnalyzer->saveHistograms(outfile);
  pfchsAnalyzer->saveHistograms(outfile);
  puppiAnalyzer->saveHistograms(outfile);

}
