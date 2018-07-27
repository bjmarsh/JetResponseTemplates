/****************************************************************

Debugger that prints contents of reco, gen jets and filter info

****************************************************************/

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"


using namespace std;

int main(int argc, char* argv[]) 
{
    // define what muon you are using; this is necessary as FWLite is not 
    // capable of reading edm::Views
    using reco::GenJet;
    using reco::GenParticle;

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

    edm::InputTag recoJetsTag ("ak4PFJetsCHS");
    edm::InputTag pfCandsTag ("particleFlow");
    edm::InputTag genJetsTag ("ak4GenJets");
    edm::InputTag genParticlesTag ("genParticles");
    edm::InputTag muonsTag ("muons");

    // edm::InputTag genJetsTag ("slimmedGenJets");
    // edm::InputTag genParticlesTag ("prunedGenParticles");


    // loop the events
    int startEvent = process.getParameter<int>("startEvent");
    int recoJetIdx = process.getParameter<int>("recoJetIdx");
    int genJetIdx = process.getParameter<int>("genJetIdx");

    TFile* inFile = TFile::Open(inputHandler_.files()[0].c_str());
    if( inFile ){

        fwlite::Event ev(inFile);
        ev.to(startEvent);
        edm::EventBase const & event = ev;
        // break loop if maximal number of events is reached 
                
        ////////////////////////
        //// ANALYSIS CODE /////
        ////////////////////////
                
        // Handle to the genjet collection
        edm::Handle<vector<reco::PFJet> > recoJets;
        event.getByLabel(recoJetsTag, recoJets);
                
        // Handle to the genjet collection
        edm::Handle<vector<GenJet> > genJets;
        event.getByLabel(genJetsTag, genJets);

        // Handle to the genparticle collection
        edm::Handle<vector<GenParticle> > genParticles;
        event.getByLabel(genParticlesTag, genParticles);

        // Handle to the genparticle collection
        edm::Handle<vector<reco::PFCandidate> > pfcands;
        event.getByLabel(pfCandsTag, pfcands);

        // Handle to the genparticle collection
        edm::Handle<vector<reco::Muon> > muons;
        event.getByLabel(muonsTag, muons);
                
        // fill a vector of genps for convenience
        vector<const GenParticle*> genps;
        for(vector<GenParticle>::const_iterator p = genParticles->begin(); p!=genParticles -> end(); ++p){
            genps.push_back(& * p);
        }

        vector<GenJet>::const_iterator gj = genJets->begin() + genJetIdx;
                    
        float gj_eta = gj->eta();
        float gj_phi = gj->phi();
        float gj_pt  = gj->pt();
        cout << "gen jet pt: " << gj_pt << ", eta: " << gj_eta << ", phi: " << gj_phi << endl;
        
        vector<const GenParticle*> constituents = gj->getGenConstituents();
        cout << "\tNumber constituents: " << constituents.size() << endl;
        
        for(unsigned int ipart=0; ipart<constituents.size(); ipart++){
            const GenParticle* p = constituents.at(ipart);
            
            // get the index for identification
            int idx = -1;
            vector<const GenParticle*>::const_iterator found = find(genps.begin(), genps.end(), p);
            if(found!= genps.end()) idx = found - genps.begin();

            cout << "\t";            
            cout << p->pdgId() << " (" << idx << ": " << p->pt() << ", " << p->eta() << ", " << p->phi() << ", " << p->status() << ") \n";
        }

        vector<reco::PFJet>::const_iterator rj = recoJets->begin() + recoJetIdx;
        float rj_eta = rj->eta();
        float rj_phi = rj->phi();
        float rj_pt = rj->pt();
        cout << "\nreco jet pt: " << rj_pt << ", eta: " << rj_eta << ", phi: " << rj_phi << endl;
        
        uint nPFCands = rj->getPFConstituents().size();
        for(uint ipf=0; ipf < nPFCands; ipf++){
            const reco::PFCandidate* p = rj->getPFConstituents().at(ipf).get();
            cout << "\t" << p->pdgId() << " (" << ipf << ": " << p->pt() << ", " << p->eta() << ", " << p->phi() << ") \n";
        
}
        cout << "\n gen particles within dR<0.4 of gen jet" << endl;
        for(vector<GenParticle>::const_iterator p = genParticles->begin(); p!= genParticles->end(); p++){
            if(p->status() != 1)
                continue;
            if(p->pt() < 20)
                continue;
            float dR = deltaR(p->eta(), p->phi(), gj_eta, gj_phi) > 0.4;
            if(dR > 0.4)
                continue;
            cout << "\t" << p->pdgId() << " (" << (int)(p-genParticles->begin()) << ": " << p->pt() << ", " << p->eta() << ", " << p->phi() << ", " << p->status() << ") \n";
            
        }

        cout << "\n reco particles within dR<0.4 of gen jet" << endl;
        for(vector<reco::PFCandidate>::const_iterator p = pfcands->begin(); p!= pfcands->end(); p++){
            if(p->pt() < 20)
                continue;
            float dR = deltaR(p->eta(), p->phi(), gj_eta, gj_phi) > 0.4;
            if(dR > 0.4)
                continue;
            cout << "\t" << p->pdgId() << " (" << (int)(p-pfcands->begin()) << ": " << p->pt() << ", " << p->eta() << ", " << p->phi() << ") \n";
            
        }


        cout << "\n muons" << endl;
        for(vector<reco::Muon>::const_iterator m = muons->begin(); m!= muons->end(); m++){
            reco::TrackRef t = m->innerTrack();
            reco::TrackRef bt = m->muonBestTrack();
            float tk_pt = t.isNonnull() ? t->pt() : -9999;
            float tk_ptErr = t.isNonnull() ? t->ptError() : -9999;
            int tk_algo = t.isNonnull() ? t->algo() : -9999;
            int tk_algoOrig = t.isNonnull() ? t->originalAlgo() : -9999;
            float m_segmComp = muon::segmentCompatibility(*m);
            float bt_pt = bt.isNonnull() ? bt->pt() : -9999;
            float bt_ptErr = bt.isNonnull() ? bt->ptError() : -9999;
            cout << "\t" << m->pdgId() << ", " << m->pt() << ", " << tk_pt << ", " << tk_ptErr << ", " << bt_pt << ", " << bt_ptErr << ", " << tk_algo << ", " << tk_algoOrig << ", " << m->type() << ", " << m_segmComp << ", " << m->isPFMuon() << endl;
            
            bool fails_seg_compatibility = m_segmComp < 0.3;
            bool fails_best_track_ptrelerr = bt_ptErr / bt_pt > 2.0;
            bool fails_inner_track_ptrelerr = tk_ptErr / tk_pt > 1.0;
            cout << "\t\t" << fails_seg_compatibility << " " << fails_best_track_ptrelerr << " " << fails_inner_track_ptrelerr << endl;

            if(fails_seg_compatibility || fails_best_track_ptrelerr || fails_inner_track_ptrelerr){
                for(vector<reco::PFCandidate>::const_iterator p = pfcands->begin(); p!= pfcands->end(); p++){
                    if(abs(p->pdgId()) != 13  || p->pt() < 100.0)
                        continue;
                    float dr = deltaR(p->eta(), p->phi(), m->eta(), m->phi());
                    if(dr < 0.001)
                        cout << "\t\t" << "Fails badMuonFilter!!" << endl;
                }
            }

        }

    }  
    // close input file
    inFile->Close();

    return 0;
}
