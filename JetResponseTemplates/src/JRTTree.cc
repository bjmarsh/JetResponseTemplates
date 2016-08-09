
#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"
#include <vector>
#include "TTree.h"
#include "TDirectory.h"

JRTTree::JRTTree(TTree *tree){
    if(tree != NULL)
        Init(tree);
}

void JRTTree::Init(TTree *tree){
    if(tree == NULL){
        this->tree = new TTree("Events","");
        b_evt_run       = this->tree->Branch("evt_run", &evt_run);
        b_evt_lumi      = this->tree->Branch("evt_lumi", &evt_lumi);
        b_evt_event     = this->tree->Branch("evt_event", &evt_event, "evt_event/l");
        b_genjet_pt     = this->tree->Branch("genjet_pt", &genjet_pt);
        b_genjet_eta    = this->tree->Branch("genjet_eta", &genjet_eta);
        b_genjet_phi    = this->tree->Branch("genjet_phi", &genjet_phi);
        b_n_genjet      = this->tree->Branch("n_genjet", &n_genjet);
        b_recojet_pt    = this->tree->Branch("recojet_pt", &recojet_pt);
        b_recojet_eta   = this->tree->Branch("recojet_eta", &recojet_eta);
        b_recojet_phi   = this->tree->Branch("recojet_phi", &recojet_phi);
        b_n_recojet     = this->tree->Branch("n_recojet", &n_recojet);
        b_genjet_flavour_bennett     = this->tree->Branch("genjet_flavour_bennett", &genjet_flavour_bennett);
        b_genjet_flavour_cmssw       = this->tree->Branch("genjet_flavour_cmssw", &genjet_flavour_cmssw);
        
        Reset();
    }else{
        this->tree = tree;
        this->tree->SetBranchAddress("evt_run", &evt_run, &b_evt_run);
        this->tree->SetBranchAddress("evt_lumi", &evt_lumi, &b_evt_lumi);
        this->tree->SetBranchAddress("evt_event", &evt_event, &b_evt_event);
        this->tree->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
        this->tree->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
        this->tree->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
        this->tree->SetBranchAddress("n_genjet", &n_genjet, &b_n_genjet);
        this->tree->SetBranchAddress("recojet_pt", &recojet_pt, &b_recojet_pt);
        this->tree->SetBranchAddress("recojet_eta", &recojet_eta, &b_recojet_eta);
        this->tree->SetBranchAddress("recojet_phi", &recojet_phi, &b_recojet_phi);
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
    genjet_pt.clear();
    genjet_eta.clear();
    genjet_phi.clear();
    n_genjet = 0;
    recojet_pt.clear();
    recojet_eta.clear();
    recojet_phi.clear();
    n_recojet = 0;
    genjet_flavour_bennett.clear();
    genjet_flavour_cmssw.clear();
}

void JRTTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}
