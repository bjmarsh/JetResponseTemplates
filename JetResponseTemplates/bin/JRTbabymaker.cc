#include <vector>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"

using namespace std;

// define what muon you are using; this is necessary as FWLite is not 
// capable of reading edm::Views
using reco::GenJet;
using reco::GenParticle;
using reco::PFJet;
using reco::Vertex;

vector<int> getCandidates(vector<const GenParticle*> constituents, vector<const GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi);
int getFlavourBennett(vector<int> match_cands, vector<const GenParticle*> genps);
int getFlavourCMSSW(vector<const GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi);

int main(int argc, char* argv[]) 
{

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    FWLiteEnabler::enable();

    // parse arguments
    if ( argc < 2 ) {
        cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
        return 0;
    }

    if( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ){
        cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << endl; exit(0);
    }
    // get the python configuration
    const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
    fwlite::InputSource inputHandler_(process); fwlite::OutputFiles outputHandler_(process);

    edm::InputTag genJetsTag ("ak4GenJets");
    edm::InputTag recoJetsTag ("ak4PFJetsCHS");
    edm::InputTag genParticlesTag ("genParticles");
    edm::InputTag fixedGridRhoTag ("fixedGridRhoFastjetAll");
    edm::InputTag vertexTag ("offlinePrimaryVertices");

    TFile *fout = new TFile(outputHandler_.file().c_str(),"RECREATE");

    JRTTree t;
    t.Init();

    // loop the events
    int ievt=0;  
    int maxEvents_( inputHandler_.maxEvents() );
    for(unsigned int iFile=0; iFile<inputHandler_.files().size(); ++iFile){

        TFile* inFile = TFile::Open(inputHandler_.files()[iFile].c_str());
        if( inFile ){

            fwlite::Event ev(inFile);
            for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
                edm::EventBase const & event = ev;
                // break loop if maximal number of events is reached 
                if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
                // simple event counter
                if(inputHandler_.reportAfter()!=0 ? (ievt>0 && ievt%inputHandler_.reportAfter()==0) : false) 
                    cout << "  processing event: " << ievt << endl;
	
                ////////////////////////
                //// ANALYSIS CODE /////
                ////////////////////////

                t.Reset();
                
                t.evt_run   = (unsigned int) ev.getRun().run();
                t.evt_lumi  = (unsigned int) ev.getLuminosityBlock().luminosityBlock();
                t.evt_event = (ULong64_t)    ev.id().event();

                // cout << "\n EVENT: " << ievt << endl;

                // Handle to the genjet collection
                edm::Handle<vector<GenJet> > genJets;
                event.getByLabel(genJetsTag, genJets);
                // Handle to the recojet collection
                edm::Handle<vector<PFJet> > recoJets;
                event.getByLabel(recoJetsTag, recoJets);
                // Handle to the genparticle collection
                edm::Handle<vector<GenParticle> > genParticles;
                event.getByLabel(genParticlesTag, genParticles);

                // Handle to the fixedGrigRho
                edm::Handle<double> fixedGridRhoHandle;
                event.getByLabel(fixedGridRhoTag, fixedGridRhoHandle);
                double fixedGridRho = *fixedGridRhoHandle.product();
                t.evt_fixgridfastjet_all_rho = fixedGridRho;

                // Handle to the vertex collection
                edm::Handle< vector<Vertex> > vertexHandle;
                event.getByLabel(vertexTag, vertexHandle);
                t.evt_nvertices = vertexHandle.product()->size();

                // fill a vector of genps for convenience
                vector<const GenParticle*> genps;
                for(vector<GenParticle>::const_iterator p = genParticles->begin(); p!=genParticles -> end(); ++p){
                    genps.push_back(& * p);
                }

                // loop genjet collection
                for(vector<GenJet>::const_iterator gj=genJets->begin(); gj!=genJets->end(); ++gj){
                    
                    float jet_eta = gj->eta();
                    float jet_phi = gj->phi();
                    float jet_pt  = gj->pt();
                    if(jet_pt < 10) continue;

                    // cout << "gen jet pt: " << jet_pt << ", eta: " << jet_eta << ", phi: " << jet_phi << endl;

                    // fill the tree variables for later filling;
                    t.genjet_pt->push_back(jet_pt);
                    t.genjet_eta->push_back(jet_eta);
                    t.genjet_phi->push_back(jet_phi);
                    
                    // do the flavour matching
                    vector<const GenParticle*> constituents = gj->getGenConstituents();
                    vector<int> match_cands = getCandidates(constituents, genps, jet_pt, jet_eta, jet_phi);
                    int flavourBennett = getFlavourBennett(match_cands, genps);
                    int flavourCMSSW = getFlavourCMSSW(genps, jet_pt, jet_eta, jet_phi);

                    t.genjet_flavour_bennett->push_back(flavourBennett);
                    t.genjet_flavour_cmssw->push_back(flavourCMSSW);
                }

                // set up JECs
                string jecDir = string(getenv("CMSSW_BASE")) + "/src/JetResponseTemplates/JetResponseTemplates/jecs/";                    
                JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(jecDir+"Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt"); 
                JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(jecDir+"Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt");
                JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(jecDir+"Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt");
                JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(jecDir+"Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt");
                vector<JetCorrectorParameters> vPar;
                vPar.push_back(*L1JetPar);
                vPar.push_back(*L2JetPar);
                vPar.push_back(*L3JetPar);
                vPar.push_back(*ResJetPar);
                
                FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);

                // loop recojet collection
                for(vector<PFJet>::const_iterator rj=recoJets->begin(); rj!=recoJets->end(); ++rj){

                    JetCorrector->setRho(fixedGridRho);
                    JetCorrector->setJetPt(rj->pt());
                    JetCorrector->setJetEta(rj->eta());
                    JetCorrector->setJetA(rj->jetArea());

                    double correction = JetCorrector->getCorrection();

                    // fill the tree variables for later filling;
                    t.recojet_pt->push_back(rj->pt() * correction);
                    t.recojet_pt_uncor->push_back(rj->pt());
                    t.recojet_eta->push_back(rj->eta());
                    t.recojet_phi->push_back(rj->phi());
                    t.recojet_area->push_back(rj->jetArea());
                    
                }

                t.n_genjet = t.genjet_pt->size();
                t.n_recojet = t.recojet_pt->size();
                
                t.Fill();
            }  
            // close input file
            inFile->Close();
        }

        // break loop if maximal number of events is reached:
        // this has to be done twice to stop the file loop as well
        if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
    }

    t.Write(fout);
    fout->Close();
    
    return 0;
}

vector<int> getCandidates(vector<const GenParticle*> constituents, vector<const GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi){

    vector<int> cands;

    float dr_cut = 0.1;
    float pt_ratio_cut = 0.3;

    for(unsigned int i=0; i<constituents.size(); i++){
        const GenParticle *p = constituents.at(i);
        float dr = sqrt((jet_eta-p->eta())*(jet_eta-p->eta()) + (jet_phi-p->phi())*(jet_phi-p->phi()));

        //ignore tops
        if(abs(p->pdgId())==6) continue;

        if(dr<dr_cut && p->pt()/jet_pt>pt_ratio_cut){
            int idx = -1;
            vector<const GenParticle*>::const_iterator found = find(genps.begin(), genps.end(), p);
            if(found!= genps.end()) idx = found - genps.begin();
            cands.push_back(idx);
        }

        for(unsigned int im=0; im < p->numberOfMothers(); im++){
            if(p->mother(im)->pdgId()==2212 && p->mother(im)->status()==4)
                continue;
            vector<const GenParticle*> mother (1,(GenParticle*)p->mother(im));
            vector<int> mcands = getCandidates(mother, genps, jet_pt, jet_eta, jet_phi);
            for(unsigned int j=0; j<mcands.size(); j++) cands.push_back(mcands.at(j));
        }            
    }

    vector<int> unique_cands;
    for(unsigned int i=0; i<cands.size(); i++){
        if(find(unique_cands.begin(),unique_cands.end(),cands.at(i)) == unique_cands.end())
            unique_cands.push_back(cands.at(i));
    }
    
    return unique_cands;

}

bool isBottomFlavoured(int id){
    id = abs(id);
    if(id == 5 || (id/100)%10 == 5 || (id/1000)%10 == 5)
        return true;
    return false;
}
bool isCharmFlavoured(int id){
    id = abs(id);
    if(id == 4 || (id/100)%10 == 4 || (id/1000)%10 == 4)
        return true;
    return false;
}

// my own method based on "jet constituents". The above method gets every gen particle
// that is part of the jet, and finds the ones that most likely drive the jet's kinematics
// (within dr=0.1, minimum pt threshold). This uses the one with the highest pt
// to determine the flavour. KEY: 5=bottom, 4=charm, 1=light-flavor, 0=leptonic
int getFlavourBennett(vector<int> match_cands, vector<const GenParticle*> genps){
    int highpt_idx = -1;
    // cout << "CANDIDATES: ";
    bool foundB = false;
    bool foundC = false;
    for(unsigned int i=0; i<match_cands.size(); i++){
        // cout << genps.at(match_cands.at(i))->pdgId() << ", ";
        if(isBottomFlavoured(genps.at(match_cands.at(i))->pdgId()))
            foundB = true;
        if(isCharmFlavoured(genps.at(match_cands.at(i))->pdgId()))
            foundC = true;        
        if(highpt_idx==-1 || genps.at(highpt_idx)->pt() < genps.at(match_cands.at(i))->pt())
            highpt_idx = match_cands.at(i);
    }
    // if(highpt_idx != -1)
    //     cout << "HIGH: " << genps.at(highpt_idx)->pdgId();
    // cout << endl;
                    
    int highpt_id;
    if(highpt_idx != -1){
        highpt_id = abs(genps.at(highpt_idx)->pdgId());
    }else{
        // just call it a gluon jet
        highpt_id = 21;
    }

    if(foundB) return 5;
    if(foundC) return 4;

    int flavour;
    if(highpt_id == 5 || (highpt_id/100)%10 == 5 || (highpt_id/1000)%10 == 5){
        // call it a b if a b quark or a b meson (xx5xx) or b baryon (5xxx)
        flavour = 5;
    }else if(highpt_id == 4 || (highpt_id/100)%10 == 4 || (highpt_id/1000)%10 == 4){
        // call it a c if a c quark or a c meson (xx4xx) or c baryon (4xxx)
        flavour = 4;
    }else if(highpt_id>=11 && highpt_id<=16){
        //leptonic jet
        flavour = 0;
    }else{
        //call it a light-flavour
        flavour = 1;
    }

    return flavour;

}

// implementation of one of the official CMSSW jet-tagging algorithms.
// Uses status=2 partons and status=3 leptons within a dr=0.3 jet cone
// If we find any b's, call it a b-jet. Else if we find any c's, call it 
// a c-jet. Else use the hightest pT parton/lepton to determine the flavour.
// Return numbers are the same as the above method.
int getFlavourCMSSW(vector<const GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi){

    float dr_cut = 0.3;

    vector<const GenParticle*> partons;
    for(unsigned int ip=0; ip<genps.size(); ip++){
        const GenParticle *p = genps.at(ip);
        int id = abs(p->pdgId());
        bool isAParton=false, isALepton=false;
        if(id==1 || id==2 || id==3 || id==4 || id==5 || id==21)
            isAParton = true;
        if(id>=11 && id<=16)
            isALepton = true;

        // // add all status 3 partons
        // if(p->status()==3 && isAParton)
        //     partons.push_back(p);

        // add non-status-3 partons with no parton daughters
        if(p->status()!=3 && p->numberOfDaughters() > 0 && isAParton){
            int npartondaughters = 0;
            for(unsigned int j=0; j<p->numberOfDaughters(); j++){
                int dID = abs(p->daughter(j)->pdgId());
                if((dID>=1 && dID<=6) || dID==21)
                    npartondaughters++;
            }
            if(npartondaughters == 0){
                partons.push_back(p);
            }
        }

        // add status 1 leptons
        if(p->status()==1 && isALepton)
            partons.push_back(p);
            ;

    }

    int match = -1;
    int match_highpt = -1;
    float max_pt = 0;
    for(unsigned int ip=0; ip<partons.size(); ip++){
        const GenParticle *p = partons.at(ip);
        int id = abs(p->pdgId());
        float dr = sqrt((jet_eta-p->eta())*(jet_eta-p->eta()) + (jet_phi-p->phi())*(jet_phi-p->phi()));
        if(dr < dr_cut){
            if(match==-1 && id==4) match = 4;
            if(             id==5) match = 5;
            if(p->pt() > max_pt){
                match_highpt = id;
                max_pt = p->pt();
            }
        }
    }

    if(match==-1) match = match_highpt;

    // cout << match << endl;

    if(match==5) return 5;  // b
    if(match==4) return 4;  // c
    if((match>=1 && match <=3) || match==21) return 1; //light-flavor
    if(match>=11 && match<=16) return 0; // leptonic
    return -1; //shouldn't ever get here.
   
}
