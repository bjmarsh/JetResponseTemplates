
#include "../interface/JRTrebalTree.h"
#include <vector>
#include "TTree.h"
#include "TDirectory.h"

JRTrebalTree::JRTrebalTree(TTree *tree){
    if(tree != NULL)
        Init(tree);
}

void JRTrebalTree::Init(TTree *tree){
    if(tree == NULL){
        // if we don't pass a tree, open in "write" mode. pointers have to be initialized,
        // and then create the branches.
        this->tree = new TTree("rebalance","");

        this->rebalanceFactors   = new vector<float>;
        this->useJet  = new vector<int>;

        b_status            = this->tree->Branch("status", &status);
        b_prescale          = this->tree->Branch("prescale", &prescale);
        b_new_met           = this->tree->Branch("new_met", &new_met, "new_met/l");
        b_rebalanceFactors  = this->tree->Branch("rebalanceFactors", &rebalanceFactors);
        b_useJet            = this->tree->Branch("useJet", &useJet);
        
        Reset();
    }else{
        // if we do pass a tree, open in "read" mode
        this->tree = tree;
        this->tree->SetMakeClass(1);
        this->tree->SetBranchAddress("status", &status, &b_status);
        this->tree->SetBranchAddress("prescale", &prescale, &b_prescale);
        this->tree->SetBranchAddress("new_met", &new_met, &b_new_met);
        this->tree->SetBranchAddress("rebalanceFactors", &rebalanceFactors, &b_rebalanceFactors);
        this->tree->SetBranchAddress("useJet", &useJet, &b_useJet);
    }
}

void JRTrebalTree::Fill(){
    tree->Fill();
}

void JRTrebalTree::Reset(){
    status = -999;
    prescale = -1;
    new_met = -1;
    rebalanceFactors->clear();
    useJet->clear();
}

void JRTrebalTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}

void JRTrebalTree::GetEntry(int i){
    this->tree->GetEntry(i);
}
