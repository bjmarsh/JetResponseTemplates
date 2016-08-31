#ifndef JRTREBALTREE
#define JRTREBALTREE

#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class JRTrebalTree {
 public:

    int               status;
    int               prescale;
    float             new_met;
    vector<float>*    rebalanceFactors = 0;
    vector<int>*      useJet = 0;

    JRTrebalTree(TTree *tree=0);
    void Init(TTree *tree=0);
    void Fill();
    void Reset();
    void Write(TDirectory *d);
    void GetEntry(int entry);

 private:
    TTree *tree;
    
    TBranch *b_status = 0;
    TBranch *b_prescale = 0;
    TBranch *b_new_met = 0;
    TBranch *b_rebalanceFactors = 0;
    TBranch *b_useJet = 0;
};

#endif
