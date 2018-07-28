

// system include files
#include <memory>
#include <vector>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TBenchmark.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"

using namespace std;

vector<int> getCandidates(vector<const reco::GenParticle*> constituents, vector<const reco::GenParticle*> genps, float jet_pt, 
                          float jet_eta, float jet_phi, vector<const reco::GenParticle*>* found);
int getFlavourBennett(vector<int> match_cands, vector<const reco::GenParticle*> genps);
int getFlavourCMSSW(vector<const reco::GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi);

//
// class declaration
//

class JRTbabymaker : public edm::one::EDAnalyzer<>  {
   public:
      explicit JRTbabymaker(const edm::ParameterSet&);
      ~JRTbabymaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

    const edm::EDGetTokenT<reco::GenJetCollection > gj_;
    const edm::EDGetTokenT<reco::PFJetCollection > rj_;
    const edm::EDGetTokenT<reco::GenParticleCollection > gp_;
    const edm::EDGetTokenT<reco::PFCandidateCollection > pfc_;
    const edm::EDGetTokenT<reco::VertexCollection > vx_;
    const edm::EDGetTokenT<reco::PFMETCollection > pfmet_;
    const edm::EDGetTokenT<reco::GenMETCollection > genmet_;
    const edm::EDGetTokenT<reco::MuonCollection> mus_;
    const edm::EDGetTokenT<double > fixedGridRho_;
    const edm::EDGetTokenT<bool> ecalTP_;
    const edm::EDGetTokenT<bool> hbheNoise_;
    const edm::EDGetTokenT<bool> hbheNoiseIso_;
    const edm::EDGetTokenT<bool> eeBadSC_;
    const edm::EDGetTokenT<bool> badPFMuon_;
    const edm::EDGetTokenT<bool> badChargedCandidate_;
    const edm::EDGetTokenT<bool> globalTightHalo2016_;
    const edm::EDGetTokenT<bool> ecalBadCalib_;

    JRTTree t;
    TFile *fout;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JRTbabymaker::JRTbabymaker(const edm::ParameterSet& iConfig) : 
    gj_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjets"))),
    rj_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfjets"))),
    gp_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
    pfc_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcands"))),
    vx_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    pfmet_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("pfmet"))),
    genmet_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genmet"))),
    mus_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    fixedGridRho_(consumes<double>(iConfig.getParameter<edm::InputTag>("fixedGridRho"))),
    ecalTP_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalTP"))),
    hbheNoise_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoise"))),
    hbheNoiseIso_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseIso"))),
    eeBadSC_(consumes<bool>(iConfig.getParameter<edm::InputTag>("eeBadSC"))),
    badPFMuon_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuon"))),
    badChargedCandidate_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChargedCandidate"))),
    globalTightHalo2016_(consumes<bool>(iConfig.getParameter<edm::InputTag>("globalTightHalo2016"))),
    ecalBadCalib_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalBadCalib")))
{

    fout = new TFile(iConfig.getParameter<string>("outFile").c_str(), "RECREATE");
    fout->cd();
    t.Init();

}


JRTbabymaker::~JRTbabymaker()
{
    t.Write(fout);
    fout->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
JRTbabymaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    t.Reset();
    
    t.evt_run   = (unsigned int) iEvent.getRun().run();
    t.evt_lumi  = (unsigned int) iEvent.getLuminosityBlock().luminosityBlock();
    t.evt_event = (ULong64_t)    iEvent.id().event();

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(gj_, genJets);

    edm::Handle<reco::PFJetCollection> recoJets;
    iEvent.getByToken(rj_, recoJets);

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(gp_, genParticles);

    edm::Handle<reco::PFCandidateCollection> pfcands;
    iEvent.getByToken(pfc_, pfcands);

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vx_, vertexHandle);
    t.evt_nvertices = vertexHandle.product()->size();

    edm::Handle<reco::PFMETCollection> pfmetHandle;
    iEvent.getByToken(pfmet_, pfmetHandle);
    reco::PFMET pfmet = pfmetHandle.product()->at(0);
    double pfmet_pt = pfmet.pt();
    double pfmet_phi = pfmet.phi();

    edm::Handle<reco::GenMETCollection> genmetHandle;
    iEvent.getByToken(genmet_, genmetHandle);
    reco::GenMET genmet = genmetHandle.product()->at(0);
    double genmet_pt = genmet.pt();
    double genmet_phi = genmet.phi();

    edm::Handle<reco::MuonCollection> muons;
    iEvent.getByToken(mus_, muons);

    edm::Handle<double> fixedGridRhoHandle;
    iEvent.getByToken(fixedGridRho_, fixedGridRhoHandle);
    double fixedGridRho = *fixedGridRhoHandle.product();
    t.evt_fixgridfastjet_all_rho = fixedGridRho;                

    edm::Handle<bool> ecalTPHandle;
    iEvent.getByToken(ecalTP_, ecalTPHandle);
    t.Flag_ecalDeadCellTriggerPrimitiveFilter = *ecalTPHandle.product();

    edm::Handle<bool> hbheNoiseHandle;
    iEvent.getByToken(hbheNoise_, hbheNoiseHandle);
    t.Flag_hbheNoiseFilter = *hbheNoiseHandle.product();

    edm::Handle<bool> hbheNoiseIsoHandle;
    iEvent.getByToken(hbheNoiseIso_, hbheNoiseIsoHandle);
    t.Flag_hbheNoiseIsoFilter = *hbheNoiseIsoHandle.product();

    edm::Handle<bool> eeBadSCHandle;
    iEvent.getByToken(eeBadSC_, eeBadSCHandle);
    t.Flag_eeBadScFilter = *eeBadSCHandle.product();

    edm::Handle<bool> badPFMuonHandle;
    iEvent.getByToken(badPFMuon_, badPFMuonHandle);
    t.Flag_badPFMuonFilter = *badPFMuonHandle.product();

    edm::Handle<bool> badChargedCandidateHandle;
    iEvent.getByToken(badChargedCandidate_, badChargedCandidateHandle);
    t.Flag_badChargedCandidateFilter = *badChargedCandidateHandle.product();

    edm::Handle<bool> globalTightHalo2016Handle;
    iEvent.getByToken(globalTightHalo2016_, globalTightHalo2016Handle);
    t.Flag_globalTightHalo2016Filter = *globalTightHalo2016Handle.product();

    edm::Handle<bool> ecalBadCalibHandle;
    iEvent.getByToken(ecalBadCalib_, ecalBadCalibHandle);
    t.Flag_ecalBadCalibFilter = *ecalBadCalibHandle.product();

    // fill a vector of genps for convenience
    vector<const reco::GenParticle*> genps;
    for(vector<reco::GenParticle>::const_iterator p = genParticles->begin(); p!=genParticles -> end(); ++p){
        genps.push_back(& * p);
    }

    // loop genjet collection
    for(reco::GenJetCollection::const_iterator gj=genJets->begin(); gj!=genJets->end(); ++gj){
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
        vector<const reco::GenParticle*> constituents = gj->getGenConstituents();
        vector<const reco::GenParticle*>* found = new vector<const reco::GenParticle*>;
        vector<int> match_cands = getCandidates(constituents, genps, jet_pt, jet_eta, jet_phi, found);
        int flavourBennett = getFlavourBennett(match_cands, genps);
        // int flavourBennett = 0;
        int flavourCMSSW = getFlavourCMSSW(genps, jet_pt, jet_eta, jet_phi);
        delete found;

        t.genjet_flavour_bennett->push_back(flavourBennett);
        t.genjet_flavour_cmssw->push_back(flavourCMSSW);

        float gj_muFrac = 0.0;
        for(uint ic=0; ic<constituents.size(); ic++){
            if(constituents.at(ic)->status()==1 && abs(constituents.at(ic)->pdgId())==13)
                gj_muFrac += constituents.at(ic)->energy();
        }
        gj_muFrac /= gj->energy();
        t.genjet_muFrac->push_back((char)round(100*gj_muFrac));
    }

    // set up JECs
    string jecDir = string(getenv("CMSSW_BASE")) + "/src/JetResponseTemplates/JetResponseTemplates/jecs/";                    
    JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(jecDir+"Fall17_17Nov2017_V4_MC_L2L3Residual_AK4PFchs.txt"); 
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(jecDir+"Fall17_17Nov2017_V4_MC_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(jecDir+"Fall17_17Nov2017_V4_MC_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(jecDir+"Fall17_17Nov2017_V4_MC_L1FastJet_AK4PFchs.txt");
    vector<JetCorrectorParameters> vPar;
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    vPar.push_back(*ResJetPar);
                
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);

    // get the x and y components of pfmet, so we can correct when we do JECs
    double pfmet_x = pfmet_pt * cos(pfmet_phi);
    double pfmet_y = pfmet_pt * sin(pfmet_phi);

    // loop recojet collection
    for(reco::PFJetCollection::const_iterator rj=recoJets->begin(); rj!=recoJets->end(); ++rj){

        JetCorrector->setRho(fixedGridRho);
        JetCorrector->setJetPt(rj->pt());
        JetCorrector->setJetEta(rj->eta());
        JetCorrector->setJetA(rj->jetArea());

        vector<float> corrs = JetCorrector->getSubCorrections();
        float correction = corrs.at(corrs.size()-1);
        // float corr_l1 = corrs.at(0);
                    
        float pt = rj->pt();
        float pt_cor = rj->pt() * correction;
        float phi = rj->phi();
        float eta = rj->eta();

        // add on the original pt and subtract back off the corrected to get corrected pfmet;
        pfmet_x += pt * cos(phi);
        pfmet_y += pt * sin(phi);
        pfmet_x -= pt_cor * cos(phi);
        pfmet_y -= pt_cor * sin(phi);

        float energy = rj->p4().energy();
        float chf = rj->chargedHadronEnergy()/energy;
        float nhf = rj->neutralHadronEnergy()/energy;
        float cef = rj->chargedEmEnergy()/energy;
        float nef = rj->neutralEmEnergy()/energy;
        float muf = rj->muonEnergy()/energy;
        float elf = rj->electronEnergy()/energy;
        float cm = rj->chargedMultiplicity();
        float nm = rj->neutralMultiplicity();

        // do jet ID
        bool isLoosePFJet = true;
        bool isTightPFJet = true;
        if(eta < 3.0){
            if(nhf >= 0.99) isLoosePFJet = false;
            if(nef >= 0.99) isLoosePFJet = false;
            if(cm + nm < 2) isLoosePFJet = false;
            if(nef >= 0.9) isTightPFJet = false;
            if(nhf >= 0.9) isTightPFJet = false;
            if(eta < 2.4){
                if(!(cm > 0))  isLoosePFJet = false;
                if(!(chf > 0)) isLoosePFJet = false;
                if(!(cef < 0.99)) isLoosePFJet = false;
            }
        }else{
            if(!(nef < 0.9)) isLoosePFJet = false;
            if(!(nm > 10)) isLoosePFJet = false;
        }
        isTightPFJet = isTightPFJet && isLoosePFJet;

        int leadingPFCandId;
        if(rj->getPFConstituents().size() > 0)
            leadingPFCandId = rj->getPFConstituents().at(0).get()->pdgId();
        else
            leadingPFCandId = 0;

        // fill the tree variables for later filling;
        t.recojet_pt->push_back(pt_cor);
        t.recojet_pt_uncor->push_back(pt);
        t.recojet_eta->push_back(eta);
        t.recojet_phi->push_back(phi);
        t.recojet_area->push_back(rj->jetArea());
        t.recojet_isLoosePFJet->push_back(isLoosePFJet);
        t.recojet_isTightPFJet->push_back(isTightPFJet);                    
        t.recojet_cemFrac->push_back((char)round(100*cef));
        t.recojet_nemFrac->push_back((char)round(100*nef));
        t.recojet_chFrac->push_back((char)round(100*chf));
        t.recojet_nhFrac->push_back((char)round(100*nhf));
        t.recojet_muFrac->push_back((char)round(100*muf));
        t.recojet_elFrac->push_back((char)round(100*elf));
        t.recojet_leadingPFCandId->push_back(leadingPFCandId);

    } // recojet loop

    //clean up
    delete ResJetPar;
    delete L3JetPar;
    delete L2JetPar;
    delete L1JetPar;
    delete JetCorrector;
    
    t.n_genjet = t.genjet_pt->size();
    t.n_recojet = t.recojet_pt->size();
    
    t.pfmet_pt = sqrt(pfmet_x*pfmet_x + pfmet_y*pfmet_y);
    t.pfmet_phi = atan2(pfmet_y, pfmet_x);
    t.pfmet_pt_uncor = pfmet_pt;
    t.pfmet_phi_uncor = pfmet_phi;
    t.genmet_pt = genmet_pt;
    t.genmet_phi = genmet_phi;
    

    // filters for badly reconstructed things
    t.Flag_badMuonFilter2016 = true;
    t.Flag_badMuonFilter2016_loose = true;
    t.Flag_badChargedCandidateFilter2016 = true;
    for(vector<reco::Muon>::const_iterator m = muons->begin(); m!= muons->end(); m++){
        reco::TrackRef tk = m->innerTrack();
        reco::TrackRef bt = m->muonBestTrack();
        float tk_pt = tk.isNonnull() ? tk->pt() : -9999;
        float tk_ptErr = tk.isNonnull() ? tk->ptError() : -9999;
        int tk_algo = tk.isNonnull() ? tk->algo() : -9999;
        int tk_algoOrig = tk.isNonnull() ? tk->originalAlgo() : -9999;
        float m_segmComp = muon::segmentCompatibility(*m);
        float bt_pt = bt.isNonnull() ? bt->pt() : -9999;
        float bt_ptErr = bt.isNonnull() ? bt->ptError() : -9999;
            
        bool fails_seg_compatibility = m_segmComp < 0.3;
        bool fails_best_track_ptrelerr = bt_ptErr / bt_pt > 2.0;
        bool fails_inner_track_ptrelerr = tk_ptErr / tk_pt > 1.0;                    

        if(fails_seg_compatibility || fails_best_track_ptrelerr || fails_inner_track_ptrelerr){
            // badMuonFilter first
            if((m->pt() > 100.0 || tk_pt > 100.0) &&   // require one of pT > 100.0
               (m->type() & (1<<1)) ){               // require global muon

                for(vector<reco::PFCandidate>::const_iterator p = pfcands->begin(); p!= pfcands->end(); p++){
                    if(abs(p->pdgId()) != 13  || p->pt() < 100.0)
                        continue;
                    float dr = deltaR(p->eta(), p->phi(), m->eta(), m->phi());
                    if(dr < 0.001){  //fails filter!
                        t.Flag_badMuonFilter2016_loose = false;
                        if (tk_algo==14 && tk_algoOrig==14)    //suspicious algorithm
                            t.Flag_badMuonFilter2016 = false;
                    }
                }
            }
            // now badChargedCandidateFilter
            if((m->pt() > 100.0 || tk_pt > 100.0) &&   // require one of pT > 100.0
               (!m->isPFMuon()) &&                   // muon is not PF muon
               (m->type() & (1<<1)) ){               // require global muon

                for(vector<reco::PFCandidate>::const_iterator p = pfcands->begin(); p!= pfcands->end(); p++){
                    if(abs(p->pdgId()) != 211)
                        continue;

                    float dr = deltaR(p->eta(), p->phi(), m->eta(), m->phi());
                    if(dr >= 0.00001)
                        continue;
                                
                    float diffPt = p->pt() - tk_pt;
                    float avgPt = 0.5*(p->pt() + tk_pt);
                    if(fabs(diffPt) / avgPt >= 0.00001)
                        continue;

                    //fails filter!
                    t.Flag_badChargedCandidateFilter2016 = false;
                }
            }
        }
    }


    t.Fill();
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
JRTbabymaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JRTbabymaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JRTbabymaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


vector<int> getCandidates(vector<const reco::GenParticle*> constituents, vector<const reco::GenParticle*> genps, float jet_pt, 
                          float jet_eta, float jet_phi, vector<const reco::GenParticle*>* found){

    vector<int> cands;

    float dr_cut = 0.1;
    float pt_ratio_cut = 0.0;

    for(unsigned int i=0; i<constituents.size(); i++){
        const reco::GenParticle *p = constituents.at(i);

        if(find(found->begin(), found->end(), p) == found->end())
            found->push_back(p);
        else
            continue;

        float dr = deltaR(jet_eta, jet_phi, p->eta(), p->phi());

        int idx = -1;
        vector<const reco::GenParticle*>::const_iterator loc = find(genps.begin(), genps.end(), p);
        if(loc!= genps.end()) idx = loc - genps.begin();


        //ignore tops
        if(abs(p->pdgId())==6) continue;

        if(dr<dr_cut && p->pt()/jet_pt>pt_ratio_cut){
            cands.push_back(idx);
        }

        for(unsigned int im=0; im < p->numberOfMothers(); im++){
            if(p->mother(im)->pdgId()==2212 && p->mother(im)->status()==4)
                continue;
            vector<const reco::GenParticle*> mother (1,(reco::GenParticle*)p->mother(im));
            vector<int> mcands = getCandidates(mother, genps, jet_pt, jet_eta, jet_phi, found);
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
int getFlavourBennett(vector<int> match_cands, vector<const reco::GenParticle*> genps){
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
// Uses status!=3 partons and status=1 leptons within a dr=0.3 jet cone
// If we find any b's, call it a b-jet. Else if we find any c's, call it 
// a c-jet. Else use the hightest pT parton/lepton to determine the flavour.
// Return numbers are the same as the above method.
int getFlavourCMSSW(vector<const reco::GenParticle*> genps, float jet_pt, float jet_eta, float jet_phi){

    float dr_cut = 0.3;

    vector<const reco::GenParticle*> partons;
    for(unsigned int ip=0; ip<genps.size(); ip++){
        const reco::GenParticle *p = genps.at(ip);
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
        const reco::GenParticle *p = partons.at(ip);
        int id = abs(p->pdgId());
        float dr = deltaR(jet_eta, jet_phi, p->eta(), p->phi());
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

//define this as a plug-in
DEFINE_FWK_MODULE(JRTbabymaker);
