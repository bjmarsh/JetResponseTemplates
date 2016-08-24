#ifndef JRTREADER
#define JRTREADER

#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TH1D.h"

#define N_PT_BINS 23
#define N_ETA_BINS 12

using namespace std;

class JRTreader {
 public:
    
    JRTreader(char *fname=0);
    ~JRTreader();
    int Init(char *fname);
    float GetResponse(float pt, float eta, bool isBjet);
    static int GetPtBin(float pt);
    static int GetEtaBin(float eta);
    void UseRawHistograms(bool use);

 private:
    vector< vector<TH1D*>* > *fits_b;
    vector< vector<TH1D*>* > *fits_nonb;
    bool useFits = true;
};

#endif
