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

    edm::InputTag genJetsTag ("ak4GenJets");
    edm::InputTag genParticlesTag ("genParticles");

    // edm::InputTag genJetsTag ("slimmedGenJets");
    // edm::InputTag genParticlesTag ("prunedGenParticles");


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

                cout << "\n EVENT: " << ievt << endl;

                // Handle to the genjet collection
                edm::Handle<vector<GenJet> > genJets;
                event.getByLabel(genJetsTag, genJets);
                // Handle to the genparticle collection
                edm::Handle<vector<GenParticle> > genParticles;
                event.getByLabel(genParticlesTag, genParticles);
                
                // fill a vector of genps for convenience
                vector<const GenParticle*> genps;
                for(vector<GenParticle>::const_iterator p = genParticles->begin(); p!=genParticles -> end(); ++p){
                    genps.push_back(& * p);
                }

                // loop genjet collection and fill histograms
                for(vector<GenJet>::const_iterator gj=genJets->begin(); gj!=genJets->end(); ++gj){
                    
                    float jet_eta = gj->eta();
                    float jet_phi = gj->phi();
                    float jet_pt  = gj->pt();
                    cout << "gen jet pt: " << jet_pt << ", eta: " << jet_eta << ", phi: " << jet_phi << endl;

                    vector<const GenParticle*> constituents = gj->getGenConstituents();
                    cout << "\tNumber constituents: " << constituents.size() << endl;

                    vector<int> flavorMatchCands;

                    for(unsigned int ipart=0; ipart<constituents.size(); ipart++){

                        const GenParticle* p = constituents.at(ipart);

                        int level = 0;
                        GenParticle *m = (GenParticle*) p;

                        while(!(m->pdgId()==2212 && m->status()==4)){

                            // get the index for identification
                            int idx = -1;
                            vector<const GenParticle*>::const_iterator found = find(genps.begin(), genps.end(), m);
                            if(found!= genps.end()) idx = found - genps.begin();

                            float dr = sqrt((jet_eta-m->eta())*(jet_eta-m->eta()) + (jet_phi-m->phi())*(jet_phi-m->phi()));
                            float pt_ratio = m->pt() / jet_pt;
                            
                            if((dr<0.1 && pt_ratio>0.3 && pt_ratio<2) && find(flavorMatchCands.begin(),flavorMatchCands.end(),idx)==flavorMatchCands.end())
                                flavorMatchCands.push_back(idx);
                            
                            cout << "\t";
                            for(int i=0; i<level; i++) cout << "   ";

                            cout << m->pdgId() << " (" << idx << ": " << m->pt() << ", " << m->eta() << ", " << m->phi() << ") \n";
                            while(m->pdgId() == m->mother()->pdgId())
                                m = (GenParticle*) m->mother();
                            m = (GenParticle*) m->mother();
                            level++;
                        }
                        cout << endl;                        
                    } 

                    cout << "\t***MATCH CANDIDATES***" << endl;
                    for(unsigned int ic=0; ic<flavorMatchCands.size(); ic++){
                        const GenParticle *c = genps.at(flavorMatchCands.at(ic));
                        cout << "\t" << c->pdgId() << " (" << flavorMatchCands.at(ic) << ": " << c->pt() << ", " << c->eta() << ", " << c->phi() << ")" << endl;
                    }
                }
                
            }  
            // close input file
            inFile->Close();
        }

        // break loop if maximal number of events is reached:
        // this has to be done twice to stop the file loop as well
        if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
    }
    
    return 0;
}
