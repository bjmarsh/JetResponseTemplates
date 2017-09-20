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
#include <TLinearFitter.h>

using namespace std;

TH1D* FitMethod1(TH1D *h);
TH1D* FitMethod2(TH1D *h, int degree=1, int window=8);

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

            TH1D *h;
            TString fname;
            
            string jetregs[] = {"bjet","nonbjet","alljet"};
            vector< vector<TH1D*>* >* histvecs[] = {JRTs_b, JRTs_l, JRTs_a};
            
            // loop over bjet, nonbjet, alljet
            for(int ij=0; ij<3; ij++){
                h = histvecs[ij]->at(ipt)->at(ieta);

                fname = Form("fit_pt%d_eta%d_%s",ipt,ieta, jetregs[ij].c_str());
                gDirectory->Delete(fname+";*");
                gDirectory->Delete(fname+"_core;*");
                gDirectory->Delete(fname+"_tail;*");
                
                // dont try to do a fit if the hist is empty
                if(h->GetEntries()==0)
                    continue;

                float zeroerror = h->GetMaximum() / 100;
                
                for(int ibin=1; ibin <= h->GetNbinsX(); ibin++){
                    if(h->GetBinContent(ibin) == 0){
                        h->SetBinContent(ibin,0.0000001);
                        h->SetBinError(ibin,zeroerror);
                    }
                }
                
                TH1D* h_fit;

                // h_fit = FitMethod1(h);

                int deg = 2;
                int window = 8;
                // bjets
                if(ij==0){
                    deg = 2;
                    window = 10;
                    if(ipt>=0 && ieta>=8){ deg=1; window=8; }
                    if(ipt>=5 && ieta>=8){ deg=2; window=8; }
                    if(ipt>=7)           { deg=2; window=6; }
                }
                // non bjets
                if(ij==1 || ij==2){
                    deg = 2;
                    window = 8;
                    if(ipt>=1 && ieta==11){ deg=1;            }
                    if(ipt>=3 && ieta==11){ deg=1; window=15; }
                    if(ipt>=4)            {      ; window=5;  }
                    if(ipt>=5 && ieta==10){ deg=2; window=6;  }
                }
                h_fit = FitMethod2(h, deg, window);

                h_fit->SetName(fname);
                h_fit->Write();

                //get core and tails
                TF1 *core = new TF1("fcore","[0]*exp(-0.5*(x-[1])/[2]*(x-[1])/[2])",0,3);
                float max = h_fit->GetBinContent(h_fit->GetMaximumBin());
                float mean = h_fit->GetBinCenter(h_fit->GetMaximumBin());
                float rms = h_fit->GetRMS();
                if(rms > 0.1)
                    rms *= 0.5;
                core->SetRange(mean-rms, mean+rms);
                core->SetParameters(max, mean, 0.1);
                core->SetParLimits(0, max-0.0005, max+0.0005);
                core->SetParLimits(1, mean-0.02, mean+0.02);
                h_fit->Fit(core, "QNR", "goff");
                core->SetRange(0,3);
                core->SetNpx(300);

                TH1D* h_core = (TH1D*)core->GetHistogram()->Clone(fname+"_core");
                TH1D* h_tail = (TH1D*)h_fit->Clone(fname+"_tail");
                h_tail->Add(h_core, -1);

                h_core->Write();
                h_tail->Write();
                   
                delete h_fit;
                delete core;
                delete h_core;
                delete h_tail;
            }



        }
    }    

    fout.Close();

}


// fit to a triple gaussian
TH1D* FitMethod1(TH1D *h){
    
    TString func = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]*exp(-0.5*((x-[4])/[5])*((x-[4])/[5])) + [6]*exp(-0.5*((x-[7])/[8])*((x-[7])/[8]))";
    // TString func = "0.3989*([0]/[2]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + [3]/[5]*exp(-0.5*((x-[4])/[5])*((x-[4])/[5])) + (1-[0]-[3])/[7]*exp(-0.5*((x-[6])/[7])*((x-[6])/[7])))";
    TF1* f = new TF1("ftmp", func, 0.0, 3.0);

    TH1D *h2 = (TH1D*)h->Clone();
    for(int i=5; i<=h2->GetNbinsX()-4; i++){
        h2->SetBinContent(i,h->Integral(i-4,i+4));
    }

    int maxbin = h2->GetMaximumBin();
    float max = h->Integral(maxbin-1, maxbin+1) / 3.0;
    float center = h->GetBinCenter(maxbin);
    float rms = h->GetRMS();

    f->SetParameters(max, center, rms, max/3, center+0.2, 0.3, max/3, center-0.2, 0.3);
    // max = max*rms/0.3989;
    // f->SetParameters(max, center, rms, max/3, center+0.2, 0.3, center-0.2, 0.3);

    f->SetParLimits(0,max-0.003,max+0.003);
    f->SetParLimits(1, center-0.02, center+0.02);

    f->SetParLimits(3,0,max/2);
    f->SetParLimits(4,0,3);
    f->SetParLimits(5,0.1,5);

    f->SetParLimits(6,0,max/2);
    f->SetParLimits(7,0,3);
    f->SetParLimits(8,0.1,5);

    // f->SetParLimits(6,0,3);
    // f->SetParLimits(7,0.1,5);

    h->Fit(f,"QNR","goff");
    f->SetNpx(300);
    TH1D* hfit = (TH1D*)f->GetHistogram()->Clone();

    delete f;
    delete h2;
    
    return hfit;
}

TH1D* FitMethod2(TH1D* h, int degree, int window){
    TH1D *h_pred = (TH1D*)h->Clone("h_pred");
    TLinearFitter *fitter = new TLinearFitter();
    if(degree==1)
        fitter->SetFormula("1 ++ x");
    else if(degree==2)
        fitter->SetFormula("1 ++ x ++ x*x");

    double minerror = h->GetMaximum()/100;
        
    for(int i=1; i<=h->GetNbinsX(); i++){
        fitter->ClearPoints();
        for(int j=i-window+1; j<=i+window-1; j++){
            if(j<1 || j>h->GetNbinsX())
                continue;
            double x = h->GetBinCenter(j);
            double y = h->GetBinContent(j);
            double dist = fabs((float)(j-i))/window;
            double w = pow(1-dist*dist*dist, 3);
            double e = max(h->GetBinError(j), minerror) / sqrt(w);
            fitter->AddPoint(&x, y, e);
        }
        fitter->Eval();

        double x = h->GetBinCenter(i);
        double pred = fitter->GetParameter(0) + fitter->GetParameter(1) * x;
        if(degree==2)
            pred += fitter->GetParameter(2)*x*x;

        h_pred->SetBinContent(i,pred);
        if(pred<0)
            h_pred->SetBinContent(i,0);
        h_pred->SetBinError(i,0);
    }

    h_pred->SetBinContent(0,0); // shouldn't need to do this, but just in case
    h_pred->SetBinContent(h_pred->GetNbinsX()+1, 0);
    TH1D *h_pred_fine = new TH1D("h_pred_fine","",300, 0, 3);
    for(int i=1; i<=h_pred_fine->GetNbinsX(); i++){
        if(i%2 == 1){
            int j = (i-1)/2;
            h_pred_fine->SetBinContent(i, 0.25*h_pred->GetBinContent(j) + 0.75*h_pred->GetBinContent(j+1));
        }
        if(i%2 == 0){
            int j = (i+1)/2;
            h_pred_fine->SetBinContent(i, 0.75*h_pred->GetBinContent(j) + 0.25*h_pred->GetBinContent(j+1));
        }
    }

    if(h_pred_fine->Integral() > 0)
        h_pred_fine->Scale(1.0/h_pred_fine->Integral()/h_pred_fine->GetBinWidth(1));

    delete fitter;
    delete h_pred;
    
    return h_pred_fine;

}
