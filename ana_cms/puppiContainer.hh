#include "NoTrees.hh"
#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include "RecoObj.hh"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//......................
class puppiContainer{
public:
    // default ctor
  puppiContainer(std::vector<RecoObj> inParticles);
    ~puppiContainer(); 
    
    std::vector<fastjet::PseudoJet> genParticles(){ return _genParticles; }
    std::vector<fastjet::PseudoJet> pfParticles(){ return _pfParticles; }    
    std::vector<fastjet::PseudoJet> pvParticles(){ return _chargedPV; }        
    std::vector<fastjet::PseudoJet> puParticles(){ return _chargedNoPV; }    
    std::vector<fastjet::PseudoJet> pfchsParticles(){ return _pfchsParticles; }    
    double goodVar(int iId,std::vector<fastjet::PseudoJet> &iParts, int iOpt);    

protected:
    std::vector<PseudoJet> _pfParticles;
    std::vector<PseudoJet> _pfchsParticles;    
    std::vector<PseudoJet> _genParticles;
    std::vector<PseudoJet> _chargedPV;
    std::vector<PseudoJet> _chargedNoPV;
};

//FASTJET_END_NAMESPACE

