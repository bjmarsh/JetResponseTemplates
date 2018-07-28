
#include "../interface/JRTTree.h"
#include <vector>
#include "TTree.h"
#include "TDirectory.h"

JRTTree::JRTTree(TTree *tree){
    if(tree != NULL)
        Init(tree);
}

void JRTTree::Init(TTree *tree){
    if(tree == NULL){
        // if we don't pass a tree, open in "write" mode. pointers have to be initialized,
        // and then create the branches.
        this->tree = new TTree("Events","");

        this->genjet_pt   = new vector<float>;
        this->genjet_eta  = new vector<float>;
        this->genjet_phi  = new vector<float>;
        this->genjet_muFrac  = new vector<char>;
        this->recojet_pt_uncor  = new vector<float>;
        this->recojet_pt  = new vector<float>;
        this->recojet_eta = new vector<float>;
        this->recojet_phi = new vector<float>;
        this->recojet_area = new vector<float>;
        this->recojet_isLoosePFJet = new vector<bool>;
        this->recojet_isTightPFJet = new vector<bool>;
        this->recojet_cemFrac = new vector<char>;
        this->recojet_nemFrac = new vector<char>;
        this->recojet_chFrac = new vector<char>;
        this->recojet_nhFrac = new vector<char>;
        this->recojet_muFrac = new vector<char>;
        this->recojet_elFrac = new vector<char>;
        this->recojet_leadingPFCandId = new vector<int>;
        this->genjet_flavour_bennett = new vector<int>;
        this->genjet_flavour_cmssw = new vector<int>;

        b_evt_run                     = this->tree->Branch("evt_run", &evt_run);
        b_evt_lumi                    = this->tree->Branch("evt_lumi", &evt_lumi);
        b_evt_event                   = this->tree->Branch("evt_event", &evt_event, "evt_event/l");
        b_Flag_badMuonFilter2016                  = this->tree->Branch("Flag_badMuonFilter2016", &Flag_badMuonFilter2016);
        b_Flag_badMuonFilter2016_loose            = this->tree->Branch("Flag_badMuonFilter2016_loose", &Flag_badMuonFilter2016_loose);
        b_Flag_badChargedCandidateFilter2016      = this->tree->Branch("Flag_badChargedCandidateFilter2016", &Flag_badChargedCandidateFilter2016);
        b_Flag_ecalDeadCellTriggerPrimitiveFilter = this->tree->Branch("Flag_ecalDeadCellTriggerPrimitiveFilter", &Flag_ecalDeadCellTriggerPrimitiveFilter);
        b_Flag_hbheNoiseFilter                    = this->tree->Branch("Flag_hbheNoiseFilter", &Flag_hbheNoiseFilter);
        b_Flag_hbheNoiseIsoFilter                 = this->tree->Branch("Flag_hbheNoiseIsoFilter", &Flag_hbheNoiseIsoFilter);
        b_Flag_eeBadScFilter                      = this->tree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter);
        b_Flag_badPFMuonFilter                    = this->tree->Branch("Flag_badPFMuonFilter", &Flag_badPFMuonFilter);
        b_Flag_badChargedCandidateFilter          = this->tree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter);
        b_Flag_globalTightHalo2016Filter          = this->tree->Branch("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
        b_Flag_ecalBadCalibFilter          = this->tree->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
        b_evt_fixgridfastjet_all_rho  = this->tree->Branch("evt_fixgridfastjet_all_rho", &evt_fixgridfastjet_all_rho);
        b_evt_nvertices               = this->tree->Branch("evt_nvertices", &evt_nvertices);
        b_pfmet_pt                    = this->tree->Branch("pfmet_pt", &pfmet_pt);
        b_pfmet_phi                   = this->tree->Branch("pfmet_phi", &pfmet_phi);
        b_pfmet_pt_uncor              = this->tree->Branch("pfmet_pt_uncor", &pfmet_pt_uncor);
        b_pfmet_phi_uncor             = this->tree->Branch("pfmet_phi_uncor", &pfmet_phi_uncor);
        b_genmet_pt                   = this->tree->Branch("genmet_pt", &genmet_pt);
        b_genmet_phi                  = this->tree->Branch("genmet_phi", &genmet_phi);
        b_genjet_pt                   = this->tree->Branch("genjet_pt", genjet_pt);
        b_genjet_eta                  = this->tree->Branch("genjet_eta", genjet_eta);
        b_genjet_phi                  = this->tree->Branch("genjet_phi", genjet_phi);
        b_genjet_muFrac               = this->tree->Branch("genjet_muFrac", genjet_muFrac);
        b_n_genjet                    = this->tree->Branch("n_genjet", &n_genjet);
        b_recojet_pt_uncor            = this->tree->Branch("recojet_pt_uncor", recojet_pt_uncor);
        b_recojet_pt                  = this->tree->Branch("recojet_pt", recojet_pt);
        b_recojet_eta                 = this->tree->Branch("recojet_eta", recojet_eta);
        b_recojet_phi                 = this->tree->Branch("recojet_phi", recojet_phi);
        b_recojet_area                = this->tree->Branch("recojet_area", recojet_area);
        b_recojet_isLoosePFJet        = this->tree->Branch("recojet_isLoosePFJet", recojet_isLoosePFJet);
        b_recojet_isTightPFJet        = this->tree->Branch("recojet_isTightPFJet", recojet_isTightPFJet);
        b_recojet_cemFrac             = this->tree->Branch("recojet_cemFrac", recojet_cemFrac);
        b_recojet_nemFrac             = this->tree->Branch("recojet_nemFrac", recojet_nemFrac);
        b_recojet_chFrac              = this->tree->Branch("recojet_chFrac", recojet_chFrac);
        b_recojet_nhFrac              = this->tree->Branch("recojet_nhFrac", recojet_nhFrac);
        b_recojet_muFrac              = this->tree->Branch("recojet_muFrac", recojet_muFrac);
        b_recojet_elFrac              = this->tree->Branch("recojet_elFrac", recojet_elFrac);
        b_recojet_leadingPFCandId     = this->tree->Branch("recojet_leadingPFCandId", recojet_leadingPFCandId);
        b_n_recojet                   = this->tree->Branch("n_recojet", &n_recojet);
        b_genjet_flavour_bennett      = this->tree->Branch("genjet_flavour_bennett", genjet_flavour_bennett);
        b_genjet_flavour_cmssw        = this->tree->Branch("genjet_flavour_cmssw", genjet_flavour_cmssw);
        
        Reset();
    }else{
        // if we do pass a tree, open in "read" mode
        this->tree = tree;
        this->tree->SetMakeClass(1);
        this->tree->SetBranchAddress("evt_run", &evt_run, &b_evt_run);
        this->tree->SetBranchAddress("evt_lumi", &evt_lumi, &b_evt_lumi);
        this->tree->SetBranchAddress("evt_event", &evt_event, &b_evt_event);
        this->tree->SetBranchAddress("Flag_badMuonFilter2016", &Flag_badMuonFilter2016, &b_Flag_badMuonFilter2016);
        this->tree->SetBranchAddress("Flag_badMuonFilter2016_loose", &Flag_badMuonFilter2016_loose, &b_Flag_badMuonFilter2016_loose);
        this->tree->SetBranchAddress("Flag_badChargedCandidateFilter2016", &Flag_badChargedCandidateFilter2016, &b_Flag_badChargedCandidateFilter2016);
        this->tree->SetBranchAddress("Flag_ecalDeadCellTriggerPrimitiveFilter", &Flag_ecalDeadCellTriggerPrimitiveFilter, &b_Flag_ecalDeadCellTriggerPrimitiveFilter);
        this->tree->SetBranchAddress("Flag_hbheNoiseFilter", &Flag_hbheNoiseFilter, &b_Flag_hbheNoiseFilter);
        this->tree->SetBranchAddress("Flag_hbheNoiseIsoFilter", &Flag_hbheNoiseIsoFilter, &b_Flag_hbheNoiseIsoFilter);
        this->tree->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
        this->tree->SetBranchAddress("Flag_badPFMuonFilter", &Flag_badPFMuonFilter, &b_Flag_badPFMuonFilter);
        this->tree->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, &b_Flag_badChargedCandidateFilter);
        this->tree->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
        this->tree->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
        this->tree->SetBranchAddress("evt_fixgridfastjet_all_rho", &evt_fixgridfastjet_all_rho, &b_evt_fixgridfastjet_all_rho);
        this->tree->SetBranchAddress("evt_nvertices", &evt_nvertices, &b_evt_nvertices);
        this->tree->SetBranchAddress("pfmet_pt", &pfmet_pt, &b_pfmet_pt);
        this->tree->SetBranchAddress("pfmet_phi", &pfmet_phi, &b_pfmet_phi);
        this->tree->SetBranchAddress("pfmet_pt_uncor", &pfmet_pt_uncor, &b_pfmet_pt_uncor);
        this->tree->SetBranchAddress("pfmet_phi_uncor", &pfmet_phi_uncor, &b_pfmet_phi_uncor);
        this->tree->SetBranchAddress("genmet_pt", &genmet_pt, &b_genmet_pt);
        this->tree->SetBranchAddress("genmet_phi", &genmet_phi, &b_genmet_phi);
        this->tree->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
        this->tree->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
        this->tree->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
        this->tree->SetBranchAddress("genjet_muFrac", &genjet_muFrac, &b_genjet_muFrac);
        this->tree->SetBranchAddress("n_genjet", &n_genjet, &b_n_genjet);
        this->tree->SetBranchAddress("recojet_pt_uncor", &recojet_pt_uncor, &b_recojet_pt_uncor);
        this->tree->SetBranchAddress("recojet_pt", &recojet_pt, &b_recojet_pt);
        this->tree->SetBranchAddress("recojet_eta", &recojet_eta, &b_recojet_eta);
        this->tree->SetBranchAddress("recojet_phi", &recojet_phi, &b_recojet_phi);
        this->tree->SetBranchAddress("recojet_area", &recojet_area, &b_recojet_area);
        this->tree->SetBranchAddress("recojet_isLoosePFJet", &recojet_isLoosePFJet, &b_recojet_isLoosePFJet);
        this->tree->SetBranchAddress("recojet_isTightPFJet", &recojet_isTightPFJet, &b_recojet_isTightPFJet);
        this->tree->SetBranchAddress("recojet_cemFrac", &recojet_cemFrac, &b_recojet_cemFrac);
        this->tree->SetBranchAddress("recojet_nemFrac", &recojet_nemFrac, &b_recojet_nemFrac);
        this->tree->SetBranchAddress("recojet_chFrac", &recojet_chFrac, &b_recojet_chFrac);
        this->tree->SetBranchAddress("recojet_nhFrac", &recojet_nhFrac, &b_recojet_nhFrac);
        this->tree->SetBranchAddress("recojet_muFrac", &recojet_muFrac, &b_recojet_muFrac);
        this->tree->SetBranchAddress("recojet_elFrac", &recojet_elFrac, &b_recojet_elFrac);
        this->tree->SetBranchAddress("recojet_leadingPFCandId", &recojet_leadingPFCandId, &b_recojet_leadingPFCandId);
        this->tree->SetBranchAddress("n_recojet", &n_recojet, &b_n_recojet);
        this->tree->SetBranchAddress("genjet_flavour_bennett", &genjet_flavour_bennett, &b_genjet_flavour_bennett);
        this->tree->SetBranchAddress("genjet_flavour_cmssw", &genjet_flavour_cmssw, &b_genjet_flavour_cmssw);
    }
}

void JRTTree::Fill(){
    tree->Fill();
}

void JRTTree::Reset(){
    evt_run = -1;
    evt_lumi = -1;
    evt_event = -1;
    Flag_badMuonFilter2016 = 1;
    Flag_badMuonFilter2016_loose = 1;
    Flag_badChargedCandidateFilter2016 = 1;
    Flag_ecalDeadCellTriggerPrimitiveFilter = 1;
    Flag_hbheNoiseFilter = 1;
    Flag_hbheNoiseIsoFilter = 1;
    Flag_eeBadScFilter = 1;
    Flag_badPFMuonFilter = 1;
    Flag_badChargedCandidateFilter = 1;
    Flag_globalTightHalo2016Filter = 1;
    Flag_ecalBadCalibFilter = 1;
    evt_fixgridfastjet_all_rho = -1;
    evt_nvertices = -1;
    pfmet_pt = -1;
    pfmet_phi = -1;
    pfmet_pt_uncor = -1;
    pfmet_phi_uncor = -1;
    genmet_pt = -1;
    genmet_phi = -1;
    genjet_pt->clear();
    genjet_eta->clear();
    genjet_phi->clear();
    genjet_muFrac->clear();
    n_genjet = 0;
    recojet_pt_uncor->clear();
    recojet_pt->clear();
    recojet_eta->clear();
    recojet_phi->clear();
    recojet_area->clear();
    recojet_isLoosePFJet->clear();
    recojet_isTightPFJet->clear();
    recojet_cemFrac->clear();
    recojet_nemFrac->clear();
    recojet_chFrac->clear();
    recojet_nhFrac->clear();
    recojet_muFrac->clear();
    recojet_elFrac->clear();
    recojet_leadingPFCandId->clear();
    n_recojet = 0;
    genjet_flavour_bennett->clear();
    genjet_flavour_cmssw->clear();
}

void JRTTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}

void JRTTree::GetEntry(int i){
    this->tree->GetEntry(i);
}
