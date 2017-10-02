#include <vector>
#include <iostream>

#include <TH1D.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTreeCache.h>
#include <TString.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TLinearFitter.h>

using namespace std;

TH1D* FitMethod1(TH1D *h, TF1* core=NULL);
TH1D* FitMethod2(TH1D *h, int degree=1, int window=8);
TH1D* FitMethod3(TH1D *h, TF1* core=NULL);

double DoubleCrystalBall(double *xp, double *par);

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
            
            string jetregs[] = {"nonbjet","bjet","alljet"};
            vector< vector<TH1D*>* >* histvecs[] = {JRTs_l, JRTs_b, JRTs_a};
            
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
                
                // Uncomment to add pseudo-content to bins with no entries
                // prevents spurious "bumps" in the fit
                float zeroerror = h->GetMaximum() / 10;                             
                for(int ibin=1; ibin <= h->GetNbinsX(); ibin++){
                    if(h->GetBinContent(ibin) == 0){
                        h->SetBinContent(ibin,0.0000001);
                        h->SetBinError(ibin,zeroerror);
                    }
                }

                // set a minimum error so as not to weight small-content bins too highly
                float minerror = h->GetMaximum() * 2.0e-3;                           
                for(int ibin=1; ibin <= h->GetNbinsX(); ibin++){
                    if(h->GetBinError(ibin) < minerror){
                        h->SetBinError(ibin,minerror);
                    }
                }
                

                TF1 *core = new TF1("fcore","[0]*exp(-0.5*(x-[1])/[2]*(x-[1])/[2])",0,3);
                TH1D* h_fit;                

                ////// METHOD 1 (triple gaussian) ///////////////////////////////////

                // h_fit = FitMethod1(h, core);

                ////// METHOD 2 (smoothing) /////////////////////////////////////////

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
                // will fit the core if the FitMethod doesn't do this already (FitMethod2)
                float max = h_fit->GetBinContent(h_fit->GetMaximumBin());
                float mean = h_fit->GetBinCenter(h_fit->GetMaximumBin());
                float rms = h_fit->GetRMS();
                // if(rms > 0.1)
                //     rms *= 0.5;
                core->SetRange(mean-rms, mean+rms);
                core->SetParameters(max, mean, 0.1);
                core->SetParLimits(0, max-0.003, max+0.003);
                core->SetParLimits(1, mean-0.02, mean+0.02);
                h_fit->Fit(core, "QNR", "goff");

                ////// METHOD 3 (double crystal ball) /////////////////////////////////

                // h_fit = FitMethod3(h, core);

                //////////////////////////////////////////////////////////////////////

                h_fit->SetName(fname);
                h_fit->Write();

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

                // return 0;
            }



        }
    }    

    fout.Close();

}


// fit to a triple gaussian
TH1D* FitMethod1(TH1D *h, TF1* core){
    
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

    if(core != NULL){
        core->SetParameter(0, f->GetParameter(0));
        core->SetParameter(1, f->GetParameter(1));
        core->SetParameter(2, f->GetParameter(2));
    }

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
        h_pred_fine->SetBinError(i, 0.01);
    }

    if(h_pred_fine->Integral() > 0)
        h_pred_fine->Scale(1.0/h_pred_fine->Integral()/h_pred_fine->GetBinWidth(1));

    delete fitter;
    delete h_pred;
    
    return h_pred_fine;

}

// fit to a DoubleCrystalBall
TH1D* FitMethod3(TH1D *h, TF1* core){
    
    TF1* f = new TF1("fdcb", DoubleCrystalBall, 0.0, 3.0, 7);

    TH1D *h2 = (TH1D*)h->Clone();
    for(int i=5; i<=h2->GetNbinsX()-4; i++){
        h2->SetBinContent(i,h->Integral(i-4,i+4));
    }

    int maxbin = h2->GetMaximumBin();
    float max = h->Integral(maxbin-1, maxbin+1) / 3.0;
    float center = h->GetBinCenter(maxbin);
    float rms = h->GetRMS();

    // fit to gaussian first
    TF1 *fgaus = new TF1("fgaus", "gaus", center-1.0*rms, center+1.0*rms);
    fgaus->SetParameters(max, center, rms);
    fgaus->SetParLimits(0, 0.97*max, 1.03*max);
    fgaus->SetParLimits(1, center-0.03, center+0.03);
    // fgaus->SetParLimits(2, 0., 1.);
    // cout << max << " " << fgaus->GetParameter(0) << " " << 0.97*max << " " << 1.03*max << endl;
    // cout << center << " " << fgaus->GetParameter(1) << " " << center-0.03 << " " << center+0.03 <<  endl;
    // double up, down;
    // fgaus->GetParLimits(0, up, down);
    // cout << rms << " " << fgaus->GetParameter(2) << " " << up << " " << down << endl;
    h->Fit(fgaus, "QNRM", "goff");
    
    f->FixParameter(0, fgaus->GetParameter(0));
    f->FixParameter(1, fgaus->GetParameter(1));
    f->FixParameter(2, fgaus->GetParameter(2));

    f->SetParameter(3, 2.);
    f->FixParameter(4, 2.);
    f->SetParameter(5, 10.);
    f->FixParameter(6, 10.);

    f->SetParLimits(3,0.5,5.);
    f->SetParLimits(5,0.5,25);

    h->Fit(f, "QNRM", "goff", 0.0, 1.0);
    f->FixParameter(3, f->GetParameter(3));
    f->FixParameter(5, f->GetParameter(5));
    f->ReleaseParameter(4);
    f->ReleaseParameter(6);
    f->SetParLimits(4,0.5,5.);
    f->SetParLimits(6,0.5,25);
    h->Fit(f, "QNRM", "goff", 1.0, 3.0);

    f->SetNpx(300);
    TH1D* hfit = (TH1D*)f->GetHistogram()->Clone();

    // cout << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << " " << f->GetParameter(5) << " " << f->GetParameter(6) << endl;

    // TCanvas *c = new TCanvas();
    // f->Draw();
    // h->Draw("PE SAME");
    // TLine line;
    // line.SetLineColor(kGreen);
    // line.SetLineWidth(2);
    // float a1 = f->GetParameter(1) - f->GetParameter(2)*f->GetParameter(3);
    // float a2 = f->GetParameter(1) + f->GetParameter(2)*f->GetParameter(4);
    // line.DrawLine(a1, 0, a1, 1);
    // line.DrawLine(a2, 0, a2, 1);
    // c->SaveAs("~/public_html/test.pdf");

    if(core != NULL){
        core->SetParameter(0, f->GetParameter(0));
        core->SetParameter(1, f->GetParameter(1));
        core->SetParameter(2, f->GetParameter(2));
    }

    delete f;
    delete h2;
    
    return hfit;
}

double DoubleCrystalBall(double *xp, double *par) {
    // par = {scale, mean, width, alpha1, alpha2, n1, n2}

    double scale = par[0];
    double mean = par[1];
    double width = par[2];
    double alpha1 = par[3];
    double alpha2 = par[4];
    double n1 = par[5];
    double n2 = par[6];

    // cout << mean << " " << width << " " << alpha1 << " " << alpha2 << endl;

    double x = xp[0];
    double t = (x-mean)/width;
    if(t >= -alpha1 && t <= alpha2){
        return scale*exp(-0.5*t*t);
    }else if(t < -alpha1){
        double A1 = pow(n1/fabs(alpha1), n1)*exp(-alpha1*alpha1/2);
        double B1 = n1/fabs(alpha1) - fabs(alpha1);
        return scale*A1*pow(B1-t, -n1);
    }else if(t > alpha2){
        double A2 = pow(n2/fabs(alpha2), n2)*exp(-alpha2*alpha2/2);
        double B2 = n2/fabs(alpha2) - fabs(alpha2);
        return scale*A2*pow(B2+t, -n2);
    }else{
        cout << "ERROR: shouldn't get here !!!!!! " << mean << " " << width << " " << alpha1 << " " << alpha2 << endl;
        return 9999;
    }
        
}
