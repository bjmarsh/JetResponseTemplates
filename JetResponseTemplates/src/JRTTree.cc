
#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"
#include <vector>
#include "TTree.h"
#include "TDirectory.h"

JRTTree::JRTTree(TTree *tree){
    this->tree = tree;
    if(tree==0){
        this->tree = new TTree("Events","");
        Init();
    }
}

void JRTTree::Init(){
    b_evt_run = tree->Branch("evt_run", &evt_run);
    b_evt_lumi = tree->Branch("evt_lumi", &evt_lumi);
    b_evt_event = tree->Branch("evt_event", &evt_event, "evt_event/l");
    b_genjet_pt = tree->Branch("genjet_pt", &genjet_pt);
    b_genjet_eta = tree->Branch("genjet_eta", &genjet_eta);
    b_genjet_phi = tree->Branch("genjet_phi", &genjet_phi);
    b_n_genjet = tree->Branch("n_genjet", &n_genjet);
    b_recojet_pt = tree->Branch("recojet_pt", &recojet_pt);
    b_recojet_eta = tree->Branch("recojet_eta", &recojet_eta);
    b_recojet_phi = tree->Branch("recojet_phi", &recojet_phi);
    b_n_recojet = tree->Branch("n_recojet", &n_recojet);
}

void JRTTree::Fill(){
    tree->Fill();
}

void JRTTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}
