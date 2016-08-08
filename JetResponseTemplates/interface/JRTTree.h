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

    vector<float> genjet_pt;

    JRTTree(TTree *tree=0);
    void Init();
    void Fill();
    void Write(TDirectory *d);

 private:
    TTree *tree;

    TBranch *b_genjet_pt;

};

#endif
