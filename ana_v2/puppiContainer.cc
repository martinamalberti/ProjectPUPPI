#include "puppiContainer.hh"

#include "fastjet/internal/base.hh"
#include "Math/ProbFunc.h"

using namespace std;
using namespace fastjet;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


puppiContainer::puppiContainer(std::vector<fastjet::PseudoJet> inParticles, std::vector<int> isPU, std::vector<int> isCh){
    
    //_pfParticles = inParticles;
    _isPU = isPU;
    _isCh = isCh;
    _isPFCHS.resize(0);
    
    _pfParticles.resize(0);
    _genParticles.resize(0);
    _pfchsParticles.resize(0);    
    _neutrals.resize(0);
    _chargedLV.resize(0);    
    _chargedPU.resize(0);    
    
    
    double etaTracker = 5.;
    double userIndex_isPuCh = -1;
    for (unsigned int i = 0; i < inParticles.size(); i++){
        float weightGEN = 0.;
        float weightPFCHS = 0.;
        
        // fill vector of pseudojets
        if (fabs(inParticles[i].eta()) < 5 && inParticles[i].pt() > 0.0){        
            
            PseudoJet curjet(inParticles[i].px(), inParticles[i].py(), inParticles[i].pz(), inParticles[i].e());
            curjet.set_user_index( -99 );
            
            if (isPU[i] == 0){
                weightGEN = 1.;
            }
            
            if ((isPU[i] == 0) || (isPU[i] == 1 && isCh[i] == 0 && fabs(inParticles[i].eta()) < etaTracker) || (isPU[i] == 1 && fabs(inParticles[i].eta()) > etaTracker)){
                weightPFCHS = 1.;
            }
            if (weightGEN > 0){
                PseudoJet curjetGEN( weightGEN*inParticles[i].px(), weightGEN*inParticles[i].py(), weightGEN*inParticles[i].pz(), weightGEN*inParticles[i].e());
                curjetGEN.set_user_index( 0 );
                _genParticles.push_back( curjetGEN );
            }
            if (weightPFCHS > 0){
                PseudoJet curjetPFCHS( weightPFCHS*inParticles[i].px(), weightPFCHS*inParticles[i].py(), weightPFCHS*inParticles[i].pz(), weightPFCHS*inParticles[i].e());            
                _pfchsParticles.push_back( curjetPFCHS );
                
                if (isCh[i] == 0){
                    _neutrals.push_back( curjet );
                    if (isPU[i] == 0) userIndex_isPuCh = -1;
                    else userIndex_isPuCh = 1;
                }
                else{
                    _chargedLV.push_back( curjet );
                    userIndex_isPuCh = 0;                    
                }
            }
            if (weightPFCHS == 0){
                _chargedPU.push_back( curjet );
                userIndex_isPuCh = 2;                
            }
            
            _isPFCHS.push_back( weightPFCHS );
            
            // need this for cleansing
            curjet.set_user_index( userIndex_isPuCh );
            _pfParticles.push_back( curjet );
        }
    }
    
    puppiWeights_pfchs.resize(0);
    puppiWeights_chLV.resize(0);
    puppiWeights_all.resize(0);
    cleansedWeights.resize(0);    
}

puppiContainer::~puppiContainer(){}

std::vector<fastjet::PseudoJet> puppiContainer::trimEvent(double vRjet, double vptcut, double vRsub, double vfcut ){
    
    //std::cout << "trimming event..." << std::endl;
    std::vector<PseudoJet> answer;
    answer.resize(0);
    
    // -- event trimming parameters -- 
    double Rjet = vRjet;
    double ptcut = vptcut;
    double Rsub = vRsub;
    double fcut = vfcut;	
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight0 = 0.0;        
        
        // --- event trimming ---
        double pt_Rjet = pt_within_R(_pfParticles,_pfParticles[i],Rjet);
        double pt_Rsub = 0;
        if (pt_Rjet >= ptcut){
            
            pt_Rsub = pt_within_R(_pfParticles,_pfParticles[i],Rsub);            
            if (pt_Rsub / pt_Rjet > fcut) weight0 = 1.0;
        }
        if (weight0 > 0){
            PseudoJet curjet( weight0*_pfParticles[i].px(), weight0*_pfParticles[i].py(), weight0*_pfParticles[i].pz(), weight0*_pfParticles[i].e());
            answer.push_back( curjet );
        }
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::cleanseEvent( double Rsub ){
    
    std::vector<PseudoJet> answer;
    answer.resize(0);
    cleansedWeights.resize(0);    
    
    // -- event cleansing parameters -- 
    double Rsub_cleansing = Rsub;
    
    std::vector<PseudoJet> charged_PU;
    std::vector<PseudoJet> charged_LV;
    for(unsigned int i = 0; i<_pfParticles.size(); i++){
        if(_isPU[i] == 1 && _isCh[i] == 1) charged_PU.push_back(_pfParticles[i]);
        if(_isPU[i] == 0 && _isCh[i] == 1) charged_LV.push_back(_pfParticles[i]);
    }    
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight1 = 0.0;
        
        // --- event cleansing ---        
        double pt_R_CPU = pt_within_R(charged_PU,_pfParticles[i],Rsub_cleansing);
        double pt_R_CLV = pt_within_R(charged_LV,_pfParticles[i],Rsub_cleansing);        
        
        if (pt_R_CLV + pt_R_CPU > 0) weight1 = pt_R_CLV / (pt_R_CLV + pt_R_CPU );
        if (weight1 > 0){
            PseudoJet curjet( weight1*_pfParticles[i].px(), weight1*_pfParticles[i].py(), weight1*_pfParticles[i].pz(), weight1*_pfParticles[i].e());
            answer.push_back( curjet );
        }
        cleansedWeights.push_back( weight1 );    
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent_V1( double Rsub, double exponent ){
    
    std::vector<PseudoJet> answer;
    answer.resize(0);
    puppiWeights_all.resize(0);
    
    // -- puppi cleansing parameters -- 
    double Rsub_cleansing = Rsub;
    
    std::vector<PseudoJet> charged_PU;
    std::vector<PseudoJet> charged_LV;
    for(unsigned int i = 0; i<_pfParticles.size(); i++){
        if(_isPU[i] == 1 && _isCh[i] == 1) charged_PU.push_back(_pfParticles[i]);
        if(_isPU[i] == 0 && _isCh[i] == 1) charged_LV.push_back(_pfParticles[i]);
    }    
    
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        double weight1 = 0.0;
        
        // --- puppi cleansing ---        
        double pt_R_CPU = ktWeight_within_R(charged_PU,_pfParticles[i],Rsub_cleansing,exponent);
        double pt_R_CLV = ktWeight_within_R(charged_LV,_pfParticles[i],Rsub_cleansing,exponent);        
        
        if (pt_R_CLV + pt_R_CPU > 0) weight1 = pt_R_CLV / (pt_R_CLV + pt_R_CPU );
        if (weight1 > 0){
            PseudoJet curjet( weight1*_pfParticles[i].px(), weight1*_pfParticles[i].py(), weight1*_pfParticles[i].pz(), weight1*_pfParticles[i].e());
            answer.push_back( curjet );
        }
        puppiWeights_all.push_back( weight1 );    
    }
    
    return answer;
}

std::vector<fastjet::PseudoJet> puppiContainer::puppiEvent(int iPU){
    
    double Rsub = 0.3;///Using 0.3 radius
    std::vector<PseudoJet> answer;
    answer.resize(0);
    std::vector<double> lVals;
    std::vector<double> lValsPV;
    std::vector<double> lValsAll;
    
    std::vector<double> lValsAPV;
    std::vector<double> lValsAPU;
    std::vector<double> lValsAAll;    
    double lAvgPV = 0; 
    double lAvgPU = 0; 
    double lAvgAll = 0; 
    for(unsigned int i=0; i<_pfParticles.size(); i++){
        
        double puppi_Rsub = puppi_within_R(_pfchsParticles,_pfParticles[i],Rsub,false);      
        double puppi_RsubPV = puppi_within_R(_chargedLV,_pfParticles[i],Rsub,false);   
        double puppi_RsubAll = puppi_within_R(_pfParticles,_pfParticles[i],Rsub,false);   
        
//        // downweighting for soft particles near hard particles
//        double pt_Rsub = pt_within_R(_pfParticles,_pfParticles[i],Rsub);
//        double weightI = _pfParticles[i].pt()/pt_Rsub;
//        
//        puppi_Rsub*=weightI;
//        puppi_RsubPV*=weightI;        
        
        lVals.push_back(puppi_Rsub);
        lValsPV.push_back(puppi_RsubPV);
        lValsAll.push_back(puppi_RsubAll);
                
        //Was applying a pt cut
        if(_pfParticles[i].pt() > 0.2 && _isCh[i]) { 
            if(_isPU[i]) lValsAPV.push_back(puppi_RsubPV);
            if(_isPU[i]) lValsAPU.push_back(puppi_Rsub);
            if(_isPU[i]) lValsAAll.push_back(puppi_RsubAll);
        }
        
    }
    
    //    //shift to make positive-definite
    //    for(int i0 = 0; i0 < lVals.size(); i0++) { 
    //        lVals[i0] -= min_lVals;
    //        lValsPV[i0] -= min_lValsPV;
    //    }
    
    //Median
    std::sort (lValsAPV.begin(),lValsAPV.end());  
    std::sort (lValsAPU.begin(),lValsAPU.end());  
    std::sort (lValsAAll.begin(),lValsAAll.end());  
    lAvgPV = lValsAPV[int(lValsAPV.size()/2.+0.5)];
    lAvgPU = lValsAPU[int(lValsAPU.size()/2.+0.5)];
    lAvgAll = lValsAAll[int(lValsAAll.size()/2.+0.5)];
    //RMS
    double lRMSPU = 0; 
    for(int i0 = 0 ;i0 < lValsAPU.size(); i0++) {lRMSPU += (lValsAPU[i0]-lAvgPU)*(lValsAPU[i0]-lAvgPU);}
    double lRMSPV = 0; 
    for(int i0 = 0 ;i0 < lValsAPV.size(); i0++) {lRMSPV += (lValsAPV[i0]-lAvgPV)*(lValsAPV[i0]-lAvgPV);}
    double lRMSAll = 0; 
    for(int i0 = 0 ;i0 < lValsAAll.size(); i0++) {lRMSAll += (lValsAAll[i0]-lAvgAll)*(lValsAAll[i0]-lAvgAll);}
    
    if(lValsAPV.size() > 0) lRMSPV/=lValsAPV.size(); 
    if(lValsAPU.size() > 0) lRMSPU/=lValsAPU.size(); 
    if(lValsAAll.size() > 0) lRMSAll/=lValsAAll.size(); 
    
    //    std::cout << "min_lVals = " << min_lVals << std::endl;
    //    std::cout << "min_lValsPV = " << min_lValsPV << std::endl;
    
    //Particles
    for(int i0 = 0; i0 < lVals.size(); i0++) { 
        
        double pWeight = lVals[i0];
        double pWeight1 = lValsPV[i0];
        double pWeight2 = lValsAll[i0];
        
        double chiWeight = (lVals[i0]-lAvgPU)/sqrt(lRMSPU)*((lVals[i0]-lAvgPU)/sqrt(lRMSPU)); 
        double chiWeight1 = (lValsPV[i0]-lAvgPV)/sqrt(lRMSPV)*(lValsPV[i0]-lAvgPV)/sqrt(lRMSPV); 
        double chiWeight2 = (lValsAll[i0]-lAvgAll)/sqrt(lRMSAll)*((lValsAll[i0]-lAvgAll)/sqrt(lRMSAll)); 
        pWeight = ROOT::Math::chisquared_cdf(chiWeight,1.);
        pWeight1 = ROOT::Math::chisquared_cdf(chiWeight1,1.);
        pWeight2 = ROOT::Math::chisquared_cdf(chiWeight2,1.);
        
        //std::cout << "pWeight2 = " << pWeight2 << std::endl;
        
        puppiWeights_pfchs.push_back(pWeight);
        puppiWeights_chLV.push_back(pWeight1);    
        puppiWeights_all.push_back(pWeight2);    
        
        int isChPU = 0;
        if (_isPU[i0]==1 && _isCh[i0]==1) isChPU = 1;
        
        if (pWeight1 > 0.80 and isChPU == 0){
            //pWeight2 = 1.;
            PseudoJet curjet(pWeight*_pfParticles[i0].px(), 
                             pWeight*_pfParticles[i0].py(), 
                             pWeight*_pfParticles[i0].pz(), 
                             pWeight*_pfParticles[i0].e());
            curjet.set_user_index(_isPU[i0]);            
            answer.push_back( curjet );
        }
    }
    
    return answer;
}


double puppiContainer::pt_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R){
    
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    double answer = 0.0;
    //std::cout << "near particles (pt) = " << near_particles.size() << std::endl;
    
    for(unsigned int i=0; i<near_particles.size(); i++){
        answer += near_particles[i].pt();
    }
    return(answer);
}

double puppiContainer::ktWeight_within_R(const vector<PseudoJet> & particles, const PseudoJet& centre, double R, double exponent){
    
    fastjet::Selector sel = fastjet::SelectorCircle(R);
    sel.set_reference(centre);
    vector<PseudoJet> near_particles = sel(particles);
    
    double sumWeight = 0.0;
    //std::cout << "near particles (kt) = " << near_particles.size() << std::endl;
    for(unsigned int i=0; i<near_particles.size(); i++){
        double deta = (centre.eta() - near_particles[i].eta());
        double dphi = (centre.phi() - near_particles[i].phi());
        double drij2 = deta*deta + dphi*dphi;
        double drij = sqrt(drij2);
        if (drij2 > 0){
            //float ktFact = min( pow(near_particles[i].pt(),2), pow(centre.pt(),2) )*drij2;
            float weight = min( pow(near_particles[i].pt(),exponent),pow(centre.pt(),exponent) ) / pow(drij,exponent);
            //std::cout << "ktFact = " << ktFact << std::endl;
            sumWeight += weight;
        }
    }
    return(sumWeight);
}

double puppiContainer::getAverage(const vector<double> & particles){
    //std::cout << "get mean" << std::endl;
    double ltotal = 0.;
    for(int i0 = 0 ;i0 < particles.size(); i0++) {ltotal += particles[i0];}
    return ltotal/particles.size();
}

double puppiContainer::getRMS(const vector<double> & particles){
    
    double average = getAverage( particles );
    double lRMS = 0;
    for(int i0 = 0 ;i0 < particles.size(); i0++) {lRMS += (particles[i0]-average)*(particles[i0]-average);}
    return lRMS/particles.size(); 
}


//FASTJET_END_NAMESPACE
