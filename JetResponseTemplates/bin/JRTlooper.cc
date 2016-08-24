#include <vector>
#include <iostream>

#include <TH1D.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTreeCache.h>
#include <TString.h>
#include <TDirectory.h>

#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"

using namespace std;

int getBin(float val, float *bins, int n_bins);
int getSampleID(const char* filename);

int main(int argc, char* argv[]) 
{

    if(argc==1){
        cout << "USAGE: JRTlooper <files to run over>" << endl;
        return 1;
    }

    // pt and eta binning
    float pt_bins[] = {0,20,30,50,80,120,170,230,300,380,470,570,680,800,1000,1300,1700,2200,2800,3500,4300,5200,6500,-1};
    int n_pt_bins = 23;
    float eta_bins[] = {0,0.3,0.5,0.8,1.1,1.4,1.7,2.3,2.8,3.2,4.1,5.0,-1};
    int n_eta_bins = 12;


    // setup histograms
    TH1D* h_n_matches = new TH1D("h_n_matches","N reco jets matched to each gen jet", 5, 0, 5);
    TH1D* h_pileup = new TH1D("h_pileup",";n vertices",40,0,40);

    vector< vector<TH1D*>* > *JRTs_b = new vector< vector<TH1D*>* >;
    vector< vector<TH1D*>* > *JRTs_l = new vector< vector<TH1D*>* >;
    vector< vector<TH1D*>* > *JRTs_a = new vector< vector<TH1D*>* >;
    for(int ipt = 0; ipt < n_pt_bins; ipt++){
        JRTs_b->push_back(new vector<TH1D*>);
        JRTs_l->push_back(new vector<TH1D*>);
        JRTs_a->push_back(new vector<TH1D*>);
        for(int ieta = 0; ieta < n_eta_bins; ieta++){
            TString hname_b = Form("JRT_pt%d_eta%d_bjet", ipt, ieta);
            TString hname_l = Form("JRT_pt%d_eta%d_nonbjet", ipt, ieta);
            TString hname_a = Form("JRT_pt%d_eta%d_alljet", ipt, ieta);
            JRTs_b->at(ipt)->push_back(new TH1D(hname_b.Data(), ";p_{T}^{reco}/p_{T}^{gen}",150,0,3));
            JRTs_l->at(ipt)->push_back(new TH1D(hname_l.Data(), ";p_{T}^{reco}/p_{T}^{gen}",150,0,3));
            JRTs_a->at(ipt)->push_back(new TH1D(hname_a.Data(), ";p_{T}^{reco}/p_{T}^{gen}",150,0,3));
            JRTs_b->at(ipt)->at(ieta)->Sumw2();
            JRTs_l->at(ipt)->at(ieta)->Sumw2();
            JRTs_a->at(ipt)->at(ieta)->Sumw2();
        }
    }


    TChain *chain = new TChain("Events");

    for(int i=1; i<argc; i++)
        chain->Add(argv[i]);

    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    cout << "Processing " << nEventsChain << " events." << endl;

    JRTTree t;

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        
        cout << currentFile->GetTitle() << endl;
        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128*1024*1024);
        t.Init(tree);

        // branch for scale1fb is not in JRTTree
        float scale;
        tree->SetBranchAddress("evt_scale1fb", &scale);

        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {
    
            // Get Event Content
            Long64_t tentry = tree->LoadTree(event);
            t.GetEntry(tentry);
            ++nEventsTotal;

            //ANALYSIS CODE HERE
            for(unsigned int igj=0; igj<t.genjet_pt->size(); igj++){
                // float gj_pt = t.genjet_pt->at(igj);
                float gj_eta = t.genjet_eta->at(igj);
                float gj_phi = t.genjet_phi->at(igj);
                
                int nmatches = 0;
                int match_idx = -1;
                for(unsigned int irj=0; irj<t.recojet_pt->size(); irj++){
                    float rj_pt = t.recojet_pt->at(irj);
                    float rj_eta = t.recojet_eta->at(irj);
                    float rj_phi = t.recojet_phi->at(irj);
                    
                    // reject recojets with pT < 10 GeV
                    if(rj_pt < 10) continue;
                    
                    float dr = sqrt((rj_eta-gj_eta)*(rj_eta-gj_eta) + (rj_phi-gj_phi)*(rj_phi-gj_phi));
                    if(dr < 0.3){
                        // found a match!
                        nmatches++;
                        match_idx = irj;
                    }
                }//recojet loop

                h_n_matches->Fill(nmatches);

                if(nmatches == 1){
                    float ratio = t.recojet_pt->at(match_idx) / t.genjet_pt->at(igj);
                    int pt_bin = getBin(t.genjet_pt->at(igj), pt_bins, n_pt_bins);
                    int eta_bin = getBin(fabs(t.genjet_eta->at(igj)), eta_bins, n_eta_bins);
                    if(pt_bin==-1 || eta_bin==-1){
                        cout << "WARNING: bad genjet at event # " << event << endl;
                        continue;
                    }
                    if(t.genjet_flavour_cmssw->at(igj) == 5){
                        JRTs_b->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                    }else if(t.genjet_flavour_cmssw->at(igj) >= 1){
                        JRTs_l->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                    }
                    
                    if(t.genjet_flavour_cmssw->at(igj) > 0)
                        JRTs_a->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);

                }//recojet loop

            }//genjet loop

            h_pileup->Fill(t.evt_nvertices, scale);

        }//event loop
    }//file loop

    TFile *fout = new TFile("JetResponseTemplates.root","RECREATE");
    TDirectory *d_aux = fout->mkdir("auxiliary");
    d_aux->cd();
    h_n_matches->Write();
    h_pileup->Write();

    vector<TDirectory*>*            pt_dirs  = new vector<TDirectory*>;
    vector< vector<TDirectory*>* >* eta_dirs = new vector< vector<TDirectory*>* >;
    for(int ipt = 0; ipt < n_pt_bins; ipt++){
        TString pt_name = Form("pt%d",ipt);
        pt_dirs->push_back(fout->mkdir(pt_name.Data()));
        eta_dirs->push_back(new vector<TDirectory*>);
        for(int ieta = 0; ieta < n_eta_bins; ieta++){
            TString eta_name = Form("pt%d_eta%d", ipt, ieta);
            eta_dirs->at(ipt)->push_back(pt_dirs->at(ipt)->mkdir(eta_name.Data()));
            eta_dirs->at(ipt)->at(ieta)->cd();

            // normalize
            float bw = JRTs_b->at(ipt)->at(ieta)->GetBinWidth(1);
            if(JRTs_b->at(ipt)->at(ieta)->Integral(0,-1) > 0)
                JRTs_b->at(ipt)->at(ieta)->Scale(1.0 / JRTs_b->at(ipt)->at(ieta)->Integral(0,-1) / bw);
            if(JRTs_l->at(ipt)->at(ieta)->Integral(0,-1) > 0)
                JRTs_l->at(ipt)->at(ieta)->Scale(1.0 / JRTs_l->at(ipt)->at(ieta)->Integral(0,-1) / bw);
            if(JRTs_a->at(ipt)->at(ieta)->Integral(0,-1) > 0)
                JRTs_a->at(ipt)->at(ieta)->Scale(1.0 / JRTs_a->at(ipt)->at(ieta)->Integral(0,-1) / bw);

            //write to file
            JRTs_b->at(ipt)->at(ieta)->Write();
            JRTs_l->at(ipt)->at(ieta)->Write();
            JRTs_a->at(ipt)->at(ieta)->Write();
        }
    }

    fout->Close();

}

// return the index of the pt or eta bin
int getBin(float val, float *bins, int n_bins){
    for(int ibin=0; ibin<n_bins; ibin++){
        if(val >= bins[ibin] && val < bins[ibin+1])
            return ibin;
    }
    if(val >= bins[n_bins-1])
        return n_bins-1;

    return -1;
}

int getSampleID(const char* filename){
    
    if(strstr(filename, "qcd_pt15to30") != NULL)     return 1;
    if(strstr(filename, "qcd_pt30to50") != NULL)     return 2;
    if(strstr(filename, "qcd_pt50to80") != NULL)     return 3;
    if(strstr(filename, "qcd_pt80to120") != NULL)    return 4;
    if(strstr(filename, "qcd_pt120to170") != NULL)   return 5;
    if(strstr(filename, "qcd_pt170to300") != NULL)   return 6;
    if(strstr(filename, "qcd_pt300to470") != NULL)   return 7;
    if(strstr(filename, "qcd_pt470to600") != NULL)   return 8;
    if(strstr(filename, "qcd_pt600to800") != NULL)   return 9;
    if(strstr(filename, "qcd_pt800to1000") != NULL)  return 10;
    if(strstr(filename, "qcd_pt1000to1400") != NULL) return 11;
    if(strstr(filename, "qcd_pt1400to1800") != NULL) return 12;
    if(strstr(filename, "qcd_pt1800to2400") != NULL) return 13;
    if(strstr(filename, "qcd_pt2400to3200") != NULL) return 14;
    if(strstr(filename, "qcd_pt3200toInf") != NULL)  return 15;

    return -1;
}
