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
    float GetRandomResponse(float pt, float eta, bool isBjet);
    float GetValue(float pt, float eta, bool isBjet, float smearfact);
    static int GetPtBin(float pt);
    static int GetEtaBin(float eta);
    static void GetModifiedBins(int ptbin, int etabin, bool isBjet, int *new_ptbin, int *new_etabin);
    void UseRawHistograms(bool use=true);


 private:
    vector< vector<TH1D*>* > *fits_b;
    vector< vector<TH1D*>* > *fits_nonb;
    bool useFits = true;
};

#endif
