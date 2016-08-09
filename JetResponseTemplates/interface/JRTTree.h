#ifndef JRTTREE
#define JRTTREE

#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class JRTTree {
 public:

    int              evt_run;
    int              evt_lumi;
    int              evt_event;
    vector<float>    genjet_pt;
    vector<float>    genjet_eta;
    vector<float>    genjet_phi;
    int              n_genjet;
    vector<float>    recojet_pt;
    vector<float>    recojet_eta;
    vector<float>    recojet_phi;
    int              n_recojet;
    vector<int>      genjet_flavour_bennett;
    vector<int>      genjet_flavour_cmssw;

    JRTTree(TTree *tree=0);
    void Init(TTree *tree=0);
    void Fill();
    void Reset();
    void Write(TDirectory *d);

 private:
    TTree *tree;
    
    TBranch *b_evt_run;
    TBranch *b_evt_lumi;
    TBranch *b_evt_event;
    TBranch *b_genjet_pt;
    TBranch *b_genjet_eta;
    TBranch *b_genjet_phi;
    TBranch *b_n_genjet;
    TBranch *b_recojet_pt;
    TBranch *b_recojet_eta;
    TBranch *b_recojet_phi;
    TBranch *b_n_recojet;
    TBranch *b_genjet_flavour_bennett;
    TBranch *b_genjet_flavour_cmssw;
};

#endif
