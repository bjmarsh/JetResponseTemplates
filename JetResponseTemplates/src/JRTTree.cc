
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
        // if we don't pass a tree, open in "write" mode. pointers have to be initialized,
        // and then create the branches.
        this->tree = new TTree("Events","");

        this->genjet_pt   = new vector<float>;
        this->genjet_eta  = new vector<float>;
        this->genjet_phi  = new vector<float>;
        this->recojet_pt_uncor  = new vector<float>;
        this->recojet_pt  = new vector<float>;
        this->recojet_eta = new vector<float>;
        this->recojet_phi = new vector<float>;
        this->recojet_area = new vector<float>;
        this->genjet_flavour_bennett = new vector<int>;
        this->genjet_flavour_cmssw = new vector<int>;

        b_evt_run       = this->tree->Branch("evt_run", &evt_run);
        b_evt_lumi      = this->tree->Branch("evt_lumi", &evt_lumi);
        b_evt_event     = this->tree->Branch("evt_event", &evt_event, "evt_event/l");
        b_evt_fixgridfastjet_all_rho = this->tree->Branch("evt_fixgridfastjet_all_rho", &evt_fixgridfastjet_all_rho);
        b_evt_nvertices = this->tree->Branch("evt_nvertices", &evt_nvertices);
        b_genjet_pt     = this->tree->Branch("genjet_pt", genjet_pt);
        b_genjet_eta    = this->tree->Branch("genjet_eta", genjet_eta);
        b_genjet_phi    = this->tree->Branch("genjet_phi", genjet_phi);
        b_n_genjet      = this->tree->Branch("n_genjet", &n_genjet);
        b_recojet_pt_uncor    = this->tree->Branch("recojet_pt_uncor", recojet_pt_uncor);
        b_recojet_pt    = this->tree->Branch("recojet_pt", recojet_pt);
        b_recojet_eta   = this->tree->Branch("recojet_eta", recojet_eta);
        b_recojet_phi   = this->tree->Branch("recojet_phi", recojet_phi);
        b_recojet_area   = this->tree->Branch("recojet_area", recojet_area);
        b_n_recojet     = this->tree->Branch("n_recojet", &n_recojet);
        b_genjet_flavour_bennett     = this->tree->Branch("genjet_flavour_bennett", genjet_flavour_bennett);
        b_genjet_flavour_cmssw       = this->tree->Branch("genjet_flavour_cmssw", genjet_flavour_cmssw);
        
        Reset();
    }else{
        // if we do pass a tree, open in "read" mode
        this->tree = tree;
        this->tree->SetMakeClass(1);
        this->tree->SetBranchAddress("evt_run", &evt_run, &b_evt_run);
        this->tree->SetBranchAddress("evt_lumi", &evt_lumi, &b_evt_lumi);
        this->tree->SetBranchAddress("evt_event", &evt_event, &b_evt_event);
        this->tree->SetBranchAddress("evt_fixgridfastjet_all_rho", &evt_fixgridfastjet_all_rho, &b_evt_fixgridfastjet_all_rho);
        this->tree->SetBranchAddress("evt_nvertices", &evt_nvertices, &b_evt_nvertices);
        this->tree->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
        this->tree->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
        this->tree->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
        this->tree->SetBranchAddress("n_genjet", &n_genjet, &b_n_genjet);
        this->tree->SetBranchAddress("recojet_pt_uncor", &recojet_pt_uncor, &b_recojet_pt_uncor);
        this->tree->SetBranchAddress("recojet_pt", &recojet_pt, &b_recojet_pt);
        this->tree->SetBranchAddress("recojet_eta", &recojet_eta, &b_recojet_eta);
        this->tree->SetBranchAddress("recojet_phi", &recojet_phi, &b_recojet_phi);
        this->tree->SetBranchAddress("recojet_area", &recojet_area, &b_recojet_area);
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
    evt_fixgridfastjet_all_rho = -1;
    evt_nvertices = -1;
    genjet_pt->clear();
    genjet_eta->clear();
    genjet_phi->clear();
    n_genjet = 0;
    recojet_pt_uncor->clear();
    recojet_pt->clear();
    recojet_eta->clear();
    recojet_phi->clear();
    recojet_area->clear();
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
