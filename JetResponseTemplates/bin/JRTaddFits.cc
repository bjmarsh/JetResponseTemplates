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
#include <TF1.h>

using namespace std;

void SetInitialParameters(TH1D *h, TF1 *f);

int main(int argc, char* argv[]) 
{

    if(argc==1){
        cout << "USAGE: JRTaddFits <template file>" << endl;
        return 1;
    }

    // pt and eta binning
    int n_pt_bins = 23;
    int n_eta_bins = 12;

    vector< vector<TH1D*>* > *JRTs_b = new vector< vector<TH1D*>* >;
    vector< vector<TH1D*>* > *JRTs_l = new vector< vector<TH1D*>* >;
    vector< vector<TH1D*>* > *JRTs_a = new vector< vector<TH1D*>* >;
    for(int ipt = 0; ipt < n_pt_bins; ipt++){
        JRTs_b->push_back(new vector<TH1D*>);
        JRTs_l->push_back(new vector<TH1D*>);
        JRTs_a->push_back(new vector<TH1D*>);
    }

    TFile fin (argv[1]);
    
    for(int ipt = 0; ipt < n_pt_bins; ipt++){
        for(int ieta = 0; ieta < n_eta_bins; ieta++){
            
            TString histname;
      
            histname = Form("pt%d/pt%d_eta%d/JRT_pt%d_eta%d_bjet",ipt,ipt,ieta,ipt,ieta);
            JRTs_b->at(ipt)->push_back((TH1D*)fin.Get(histname));
            JRTs_b->at(ipt)->at(ieta)->SetDirectory(0);
            
            histname = Form("pt%d/pt%d_eta%d/JRT_pt%d_eta%d_nonbjet",ipt,ipt,ieta,ipt,ieta);
            JRTs_l->at(ipt)->push_back((TH1D*)fin.Get(histname));
            JRTs_l->at(ipt)->at(ieta)->SetDirectory(0);

            histname = Form("pt%d/pt%d_eta%d/JRT_pt%d_eta%d_alljet",ipt,ipt,ieta,ipt,ieta);
            JRTs_a->at(ipt)->push_back((TH1D*)fin.Get(histname));
            JRTs_a->at(ipt)->at(ieta)->SetDirectory(0);
        }
    }

    fin.Close();

    TFile fout (argv[1], "UPDATE");

    for(int ipt = 0; ipt < n_pt_bins; ipt++){
        for(int ieta = 0; ieta < n_eta_bins; ieta++){

            TString dirname;
            dirname = Form("pt%d/pt%d_eta%d",ipt,ipt,ieta);
            
            fout.cd(dirname);

            TF1 *f;
            TH1D *h;
            TString fname;
            // TString func = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]*exp(-0.5*((x-[4])/[5])*((x-[4])/[5]))";
            TString func = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]*exp(-0.5*((x-[4])/[5])*((x-[4])/[5])) + [6]*exp(-0.5*((x-[7])/[8])*((x-[7])/[8]))";
            
            string jetregs[] = {"bjet","nonbjet","alljet"};
            vector< vector<TH1D*>* >* histvecs[] = {JRTs_b, JRTs_l, JRTs_a};

            for(int ij=0; ij<3; ij++){
                h = histvecs[ij]->at(ipt)->at(ieta);
                if(h->GetEntries()==0)
                    continue;
                
                for(int ibin=1; ibin <= h->GetNbinsX(); ibin++){
                    if(h->GetBinContent(ibin) == 0){
                        h->SetBinContent(ibin,0.0000001);
                        h->SetBinError(ibin,0.001);
                    }
                }

                fname = Form("fit_pt%d_eta%d_%s",ipt,ieta, jetregs[ij].c_str());
                gDirectory->Delete(fname+";*");
                f = new TF1(fname, func, 0.0, 3.0);
                SetInitialParameters(h, f);
                    h->Fit(f,"QNR","goff");
                    f->Write();
                delete f;
            }

        }
    }    

    fout.Close();

}


void SetInitialParameters(TH1D *h, TF1 *f){
    
    TH1D *h2 = (TH1D*)h->Clone();
    for(int i=5; i<=h2->GetNbinsX()-4; i++){
        h2->SetBinContent(i,h->Integral(i-4,i+4));
    }

    int maxbin = h2->GetMaximumBin();
    float max = h->Integral(maxbin-1, maxbin+1) / 3.0;
    float center = h->GetBinCenter(maxbin);
    float rms = h->GetRMS();
    f->SetParameters(max, center, rms, max/3, center+0.2, 0.3, max/3, center-0.2, 0.3);

    f->SetParLimits(0,max-0.003,max+0.003);
    f->SetParLimits(1, center-0.02, center+0.02);

    f->SetParLimits(3,0,max/2);
    f->SetParLimits(4,0,3);
    f->SetParLimits(5,0.1,5);

    f->SetParLimits(6,0,max/2);
    f->SetParLimits(7,0,3);
    f->SetParLimits(8,0.1,5);

    delete h2;
}
