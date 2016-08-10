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

    int               evt_run;
    int               evt_lumi;
    ULong64_t         evt_event;
    double            evt_fixgridfastjet_all_rho;
    double            evt_nvertices;
    vector<float>*    genjet_pt = 0;
    vector<float>*    genjet_eta = 0;
    vector<float>*    genjet_phi = 0;
    int               n_genjet;
    vector<float>*    recojet_pt_uncor = 0;
    vector<float>*    recojet_pt = 0;
    vector<float>*    recojet_eta = 0;
    vector<float>*    recojet_phi = 0;
    vector<float>*    recojet_area = 0;
    int               n_recojet;
    vector<int>*      genjet_flavour_bennett = 0;
    vector<int>*      genjet_flavour_cmssw = 0;

    JRTTree(TTree *tree=0);
    void Init(TTree *tree=0);
    void Fill();
    void Reset();
    void Write(TDirectory *d);
    void GetEntry(int entry);

 private:
    TTree *tree;
    
    TBranch *b_evt_run = 0;
    TBranch *b_evt_lumi = 0;
    TBranch *b_evt_fixgridfastjet_all_rho = 0;
    TBranch *b_evt_nvertices = 0;
    TBranch *b_evt_event = 0;
    TBranch *b_genjet_pt = 0;
    TBranch *b_genjet_eta = 0;
    TBranch *b_genjet_phi = 0;
    TBranch *b_n_genjet = 0;
    TBranch *b_recojet_pt_uncor = 0;
    TBranch *b_recojet_pt = 0;
    TBranch *b_recojet_eta = 0;
    TBranch *b_recojet_phi = 0;
    TBranch *b_recojet_area = 0;
    TBranch *b_n_recojet = 0;
    TBranch *b_genjet_flavour_bennett = 0;
    TBranch *b_genjet_flavour_cmssw = 0;
};

#endif
