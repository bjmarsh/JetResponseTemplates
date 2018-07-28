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
    bool              Flag_badMuonFilter2016;
    bool              Flag_badMuonFilter2016_loose;
    bool              Flag_badChargedCandidateFilter2016;
    bool              Flag_ecalDeadCellTriggerPrimitiveFilter;
    bool              Flag_hbheNoiseFilter;
    bool              Flag_hbheNoiseIsoFilter;
    bool              Flag_eeBadScFilter;
    bool              Flag_badPFMuonFilter;
    bool              Flag_badChargedCandidateFilter;
    bool              Flag_globalTightHalo2016Filter;
    bool              Flag_ecalBadCalibFilter;
    double            evt_fixgridfastjet_all_rho;
    double            evt_nvertices;
    double            pfmet_pt;
    double            pfmet_phi;
    double            pfmet_pt_uncor;
    double            pfmet_phi_uncor;
    double            genmet_pt;
    double            genmet_phi;
    vector<float>*    genjet_pt = 0;
    vector<float>*    genjet_eta = 0;
    vector<float>*    genjet_phi = 0;
    vector<char>*     genjet_muFrac = 0;
    int               n_genjet;
    vector<float>*    recojet_pt_uncor = 0;
    vector<float>*    recojet_pt = 0;
    vector<float>*    recojet_eta = 0;
    vector<float>*    recojet_phi = 0;
    vector<float>*    recojet_area = 0;
    vector<bool>*     recojet_isLoosePFJet = 0;
    vector<bool>*     recojet_isTightPFJet = 0;
    vector<char>*     recojet_cemFrac = 0;
    vector<char>*     recojet_nemFrac = 0;
    vector<char>*     recojet_chFrac = 0;
    vector<char>*     recojet_nhFrac = 0;
    vector<char>*     recojet_muFrac = 0;
    vector<char>*     recojet_elFrac = 0;
    vector<int>*      recojet_leadingPFCandId = 0;
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
    TBranch *b_evt_event = 0;
    TBranch *b_Flag_badMuonFilter2016 = 0;
    TBranch *b_Flag_badMuonFilter2016_loose = 0;
    TBranch *b_Flag_badChargedCandidateFilter2016 = 0;
    TBranch *b_Flag_ecalDeadCellTriggerPrimitiveFilter = 0;
    TBranch *b_Flag_hbheNoiseFilter = 0;
    TBranch *b_Flag_hbheNoiseIsoFilter = 0;
    TBranch *b_Flag_eeBadScFilter = 0;
    TBranch *b_Flag_badPFMuonFilter = 0;
    TBranch *b_Flag_badChargedCandidateFilter = 0;
    TBranch *b_Flag_globalTightHalo2016Filter = 0;
    TBranch *b_Flag_ecalBadCalibFilter = 0;
    TBranch *b_evt_fixgridfastjet_all_rho = 0;
    TBranch *b_evt_nvertices = 0;
    TBranch *b_pfmet_pt = 0;
    TBranch *b_pfmet_phi = 0;
    TBranch *b_pfmet_pt_uncor = 0;
    TBranch *b_pfmet_phi_uncor = 0;
    TBranch *b_genmet_pt = 0;
    TBranch *b_genmet_phi = 0;
    TBranch *b_genjet_pt = 0;
    TBranch *b_genjet_eta = 0;
    TBranch *b_genjet_phi = 0;
    TBranch *b_genjet_muFrac = 0;
    TBranch *b_n_genjet = 0;
    TBranch *b_recojet_pt_uncor = 0;
    TBranch *b_recojet_pt = 0;
    TBranch *b_recojet_eta = 0;
    TBranch *b_recojet_phi = 0;
    TBranch *b_recojet_area = 0;
    TBranch *b_recojet_isLoosePFJet = 0;
    TBranch *b_recojet_isTightPFJet = 0;
    TBranch *b_recojet_cemFrac = 0;
    TBranch *b_recojet_nemFrac = 0;
    TBranch *b_recojet_chFrac = 0;
    TBranch *b_recojet_nhFrac = 0;
    TBranch *b_recojet_muFrac = 0;
    TBranch *b_recojet_elFrac = 0;
    TBranch *b_recojet_leadingPFCandId = 0;
    TBranch *b_n_recojet = 0;
    TBranch *b_genjet_flavour_bennett = 0;
    TBranch *b_genjet_flavour_cmssw = 0;
};

#endif
