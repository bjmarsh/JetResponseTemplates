
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
    b_genjet_pt = tree->Branch("genjet_pt", &genjet_pt);
}

void JRTTree::Fill(){
    tree->Fill();
}

void JRTTree::Write(TDirectory *d){
    d->cd();
    tree->Write();
}
