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

TH1D* FitMethod1(TH1D *h, TF1* core=NULL, bool verbose=false);
TH1D* FitMethod2(TH1D *h, int degree=1, int window=8);
TH1D* FitMethod3(TH1D *h, TF1* core=NULL);
TH1D* FitMethod4(TH1D *h, TF1* core=NULL);
TH1D* StraightTemplate(TH1D *h, TF1* core=NULL);

double DoubleCrystalBall(double *xp, double *par);
double BennettFunc(double *xp, double *par);

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

    // for(int ipt = 0; ipt <= 0; ipt++){
    //     for(int ieta = 1; ieta <= 1; ieta++){

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
                float zeroerror = h->GetMaximum() / 1000;   
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
                float coresigma = 999;
                
                ////// METHOD 1 (triple gaussian) ///////////////////////////////////

                // bool verbose = false;
                // if(ij==0 && ipt==8 && ieta==6)
                //     verbose = true;
                // h_fit = FitMethod1(h, core, verbose);

                ////// METHOD 2 (smoothing) /////////////////////////////////////////

                // int deg = 2;
                // int window = 8;
                // // bjets
                // if(ij==0){
                //     deg = 2;
                //     window = 10;
                //     if(ipt>=0 && ieta>=8){ deg=1; window=8; }
                //     if(ipt>=5 && ieta>=8){ deg=2; window=8; }
                //     if(ipt>=7)           { deg=2; window=6; }
                // }
                // // non bjets
                // if(ij==1 || ij==2){
                //     deg = 2;
                //     window = 8;
                //     if(ipt>=1 && ieta==11){ deg=1;            }
                //     if(ipt>=3 && ieta==11){ deg=1; window=15; }
                //     if(ipt>=4)            {      ; window=5;  }
                //     if(ipt>=5 && ieta==10){ deg=2; window=6;  }
                // }
                // h_fit = FitMethod2(h, deg, window);
                // // will fit the core if the FitMethod doesn't do this already (FitMethod2)
                // float max = h_fit->GetBinContent(h_fit->GetMaximumBin());
                // float mean = h_fit->GetBinCenter(h_fit->GetMaximumBin());
                // float rms = h_fit->GetRMS();
                // // if(rms > 0.1)
                // //     rms *= 0.5;
                // core->SetRange(mean-rms, mean+rms);
                // core->SetParameters(max, mean, 0.1);
                // core->SetParLimits(0, max-0.003, max+0.003);
                // core->SetParLimits(1, mean-0.02, mean+0.02);
                // h_fit->Fit(core, "QNR", "goff");

                ////// METHOD 3 (double crystal ball) /////////////////////////////////

                // h_fit = FitMethod3(h, core);

                //////////////////////////////////////////////////////////////////////

                ////// METHOD 4 (bennett func) /////////////////////////////////

                // h_fit = FitMethod4(h, core);

                //////////////////////////////////////////////////////////////////////

                ////// METHOD 5 (straight template w/ core fit) /////////////////////////////////

                h_fit = StraightTemplate(h, core);
                coresigma = 2.0;

                //////////////////////////////////////////////////////////////////////

                ////// METHOD 6 (bennett on left tail, straight template on right) /////////////////////////////////

                // TH1D* h_fit2 = StraightTemplate(h, NULL);
                // h_fit = FitMethod4(h, core);

                // int lowidx = 101;
                // if(ij==1) lowidx = 1;  // for b-jets, use straight templates everywhere.
                // for(int i=lowidx; i<=300; i++){
                //     h_fit->SetBinContent(i, h_fit2->GetBinContent(i));
                //     h_fit->SetBinError(i, h_fit2->GetBinError(i));
                // }

                // delete h_fit2;
                
                //////////////////////////////////////////////////////////////////////

                h_fit->SetName(fname);

                core->SetRange(0,3);
                core->SetNpx(300);

                TH1D* h_core = (TH1D*)core->GetHistogram()->Clone(fname+"_core");
                float mean=h_core->GetMean(), rms=h_core->GetRMS();
                for(int i=1; i<h_core->GetNbinsX()+1; i++){
                    float sigma = fabs(h_core->GetBinCenter(i) - mean)/rms;
                    if (coresigma-1 <= sigma  && sigma < coresigma)
                        h_core->SetBinContent(i, (coresigma-sigma) * h_core->GetBinContent(i));
                    if(coresigma <= sigma)                            
                        h_core->SetBinContent(i, 0);
                    h_fit->SetBinError(i, 0);
                    h_core->SetBinError(i, 0);
                }

                TH1D* h_tail = (TH1D*)h_fit->Clone(fname+"_tail");
                h_tail->SetLineWidth(2);
                h_tail->SetLineColor(kGreen+2);
                h_tail->Add(h_core, -1);

                h_fit->Add(h_tail, -1);
                for(int i=1; i<h_tail->GetNbinsX()+1; i++){
                    float sigma = fabs(h_core->GetBinCenter(i) - mean)/rms;
                    if (sigma < coresigma-1)
                        h_tail->SetBinContent(i, 0);
                    h_tail->SetBinError(i, 0);                        
                }
                h_fit->Add(h_tail);

                h_fit->Write();
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
TH1D* FitMethod1(TH1D *h, TF1* core, bool verbose){
    
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


    // fit to gaussian first
    TF1 *fgaus = new TF1("fgaus", "gaus", center-1.0*rms, center+1.0*rms);
    fgaus->SetParameters(max, center, rms);
    fgaus->SetParLimits(0, 0.97*max, 1.03*max);
    fgaus->SetParLimits(1, center-0.03, center+0.03);
    h->Fit(fgaus, "QNRM", "goff");

    max = fgaus->GetParameter(0);
    center = fgaus->GetParameter(1);
    rms = fgaus->GetParameter(2);

    // f->SetParameters(max, center, rms, max/3, center-rms, rms/2, max/3, center+rms, rms/2);
    f->SetParameters(max, center, rms, max/3, center-rms, rms/2, max/3, center+rms, rms);
    // max = max*rms/0.3989;
    // f->SetParameters(max, center, rms, max/3, center+0.2, 0.3, center-0.2, 0.3);

    f->SetParLimits(0,max-0.003,max+0.003);
    f->SetParLimits(1, center-0.02, center+0.02);

    f->SetParLimits(3,0,max/2);
    f->SetParLimits(4,0.4,1-0.3*rms);
    f->SetParLimits(5,0.1,min(2.0*rms,0.5));

    f->SetParLimits(6,0,max/2);
    f->SetParLimits(7,1+0.5*rms,2);
    f->SetParLimits(8,0.1,2.0*rms);

    // f->SetParLimits(6,0,3);
    // f->SetParLimits(7,0.1,5);

    h->Fit(f,"QNR","goff");
    f->SetNpx(300);
    TH1D* hfit = (TH1D*)f->GetHistogram()->Clone();

    if(verbose){
        cout << "-------------------------------------------------------------" << endl;
        cout << max << " " << center << " " << rms << endl;
        cout << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << " " << 
            f->GetParameter(5) << " " << f->GetParameter(6) << " " << f->GetParameter(7) << " " << f->GetParameter(8) << endl;
    }

    if(core != NULL){
        core->SetParameter(0, f->GetParameter(0));
        core->SetParameter(1, f->GetParameter(1));
        core->SetParameter(2, f->GetParameter(2));
    }

    // core->SetParameters(max, center, rms);
    // core->SetParLimits(0, max-0.0005, max+0.0005);
    // core->SetParLimits(1, center-0.02, center+0.02);
    // if(rms > 0.1) rms *= 0.5;
    // hfit->Fit(core, "QN", "goff", center-rms, center+rms);

    delete f;
    delete fgaus;
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

    if(core != NULL){
        core->SetParameter(0, f->GetParameter(0));
        core->SetParameter(1, f->GetParameter(1));
        core->SetParameter(2, f->GetParameter(2));
    }

    delete f;
    delete h2;
    
    return hfit;
}
// fit to a BennettFunc
TH1D* FitMethod4(TH1D *h, TF1* core){
    
    TF1* f = new TF1("fbennett", BennettFunc, 0.0, 3.0, 7);

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
    h->Fit(fgaus, "QNRM", "goff");
    
    f->FixParameter(0, fgaus->GetParameter(0));
    f->FixParameter(1, fgaus->GetParameter(1));
    f->FixParameter(2, fgaus->GetParameter(2));

    // fit left tail
    f->SetParameter(3, 2.);
    f->FixParameter(4, 2.);
    f->SetParameter(5, rms / 2.);
    f->FixParameter(6, rms / 2.);
    f->SetParLimits(3,0.5,5.);
    f->SetParLimits(5,0.,3.);
    h->Fit(f, "QNRM", "goff", 0.0, 1.0);

    // fit right tail
    f->FixParameter(3, f->GetParameter(3));
    f->ReleaseParameter(4);
    f->FixParameter(5, f->GetParameter(5));
    f->ReleaseParameter(6);
    f->SetParLimits(4,0.5,5.);
    f->SetParLimits(6,0.,3.);
    h->Fit(f, "QNRM", "goff", 1.0, 3.0);

    f->SetNpx(300);
    TH1D* hfit = (TH1D*)f->GetHistogram()->Clone();

    // cout << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    // double xp[] = {0.5};
    // double par[] = {1.83, 1.01, 0.19, 2, 2};
    // cout << BennettFunc(xp, par) << endl;
    // cout << f->Eval(0.5) << endl;

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

TH1D* StraightTemplate(TH1D* h, TF1* core){

    // remove spikes from low stats
    TH1D *h2 = (TH1D*)h->Clone();
    for(int i=2; i<=h2->GetNbinsX()-1; i++){
        if(h->GetBinContent(i-1)<1e-5 && h->GetBinContent(i+1)<1e-5)
            h2->SetBinContent(i,0);
    }

    TH1D* hfit = new TH1D("hfit", "", 300,0,3);
    hfit->SetBinContent(1, h2->GetBinContent(1) / 2.0);
    hfit->SetBinError(1, h2->GetBinError(1) / 2.0);
    hfit->SetBinContent(300, h2->GetBinContent(150) / 2.0);
    hfit->SetBinError(300, h2->GetBinError(150) / 2.0);
    for(int i=2; i<300; i++){
        if(i%2==0){
            hfit->SetBinContent(i, 3./4*h2->GetBinContent(i/2) + 1./4*h2->GetBinContent(i/2+1));
            hfit->SetBinError(i, 3./4*h2->GetBinError(i/2) + 1./4*h2->GetBinError(i/2+1));
        }else{
            hfit->SetBinContent(i, 1./4*h2->GetBinContent((i-1)/2) + 3./4*h2->GetBinContent((i-1)/2+1));
            hfit->SetBinError(i, 1./4*h2->GetBinError((i-1)/2) + 3./4*h2->GetBinError((i-1)/2+1));
        }        
    }

    hfit->Scale(1.0 / hfit->Integral() / hfit->GetBinWidth(1));

    if(core != NULL){
        float max = hfit->GetBinContent(hfit->GetMaximumBin());
        float mean = hfit->GetBinCenter(hfit->GetMaximumBin());
        float rms = hfit->GetRMS();
        // if(rms > 0.1)
        //     rms *= 0.5;
        core->SetRange(mean-rms, mean+rms);
        core->SetParameters(max, mean, min(rms,0.1f));
        core->SetParLimits(0, max-0.003, max+0.003);
        core->SetParLimits(1, mean-0.02, mean+0.02);
        hfit->Fit(core, "QNR", "goff");
    }

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

double BennettFunc(double *xp, double *par) {
    // par = {scale, mean, width, alpha1, alpha2, beta1, beta2}

    double scale = par[0];
    double mean = par[1];
    double width = par[2];
    double alpha1 = par[3];
    double alpha2 = par[4];
    double beta1 = par[5];
    double beta2 = par[6];

    double x = xp[0];
    double t = (x-mean)/width;
    if(t >= -alpha1 && t <= alpha2){
        return scale*exp(-0.5*t*t);
    }else if(t < -alpha1){
        double A1 = scale * exp(-alpha1*alpha1/2.0);
        // double B1 = width / alpha1;
        double B1 = beta1;
        return A1*exp((x + alpha1*width - mean) / B1);
    }else if(t > alpha2){
        double A2 = scale * exp(-alpha2*alpha2/2.0);
        // double B2 = width / alpha2;
        double B2 = beta2;
        return A2*exp(-(x - alpha2*width - mean) / B2);
    }else{
        cout << "ERROR: shouldn't get here !!!!!! " << mean << " " << width << " " << alpha1 << " " << alpha2 << endl;
        return 9999;
    }
        
}

