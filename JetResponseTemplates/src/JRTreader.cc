
// #include "JetResponseTemplates/JetResponseTemplates/interface/JRTreader.h"
#include "../interface/JRTreader.h"
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"

JRTreader::JRTreader(char *fname){
    if(fname != NULL){
        Init(fname);
    }
}

JRTreader::~JRTreader(){
    if(fits_b != NULL){
        for(int i=0; i<N_PT_BINS; i++){
            if(fits_b->at(i) != NULL){
                for(int j=0; j<N_ETA_BINS; j++){
                    if(fits_b->at(i)->at(j) != NULL)
                        delete fits_b->at(i)->at(j);
                }
                delete fits_b->at(i);
            }
        }
        delete fits_b;
    }
    if(fits_nonb != NULL){
        for(int i=0; i<N_PT_BINS; i++){
            if(fits_nonb->at(i) != NULL){
                for(int j=0; j<N_ETA_BINS; j++){
                    if(fits_nonb->at(i)->at(j) != NULL)
                        delete fits_nonb->at(i)->at(j);
                }
                delete fits_nonb->at(i);
            }
        }
        delete fits_nonb;
    }
}

int JRTreader::Init(char *fname){
    TFile *f = TFile::Open(fname);
    if(f->IsZombie())
        return 0;

    fits_b = new vector< vector<TH1D*>* >;
    fits_nonb = new vector< vector<TH1D*>* >;

    TString prefix = "fit";
    if(!useFits)
        prefix = "JRT";

    for(int ipt=0; ipt<N_PT_BINS; ipt++){
        fits_b->push_back(new vector<TH1D*>);
        fits_nonb->push_back(new vector<TH1D*>);
        for(int ieta=0; ieta<N_ETA_BINS; ieta++){
            TString fitname;
            fitname = Form("pt%d/pt%d_eta%d/%s_pt%d_eta%d_bjet",ipt,ipt,ieta,prefix.Data(),ipt,ieta);
            fits_b->at(ipt)->push_back((TH1D*)f->Get(fitname));
            if(fits_b->at(ipt)->at(ieta) != NULL) 
               fits_b->at(ipt)->at(ieta)->SetDirectory(0);

            fitname = Form("pt%d/pt%d_eta%d/%s_pt%d_eta%d_nonbjet",ipt,ipt,ieta,prefix.Data(),ipt,ieta);
            fits_nonb->at(ipt)->push_back((TH1D*)f->Get(fitname));
            if(fits_nonb->at(ipt)->at(ieta) != NULL)
                fits_nonb->at(ipt)->at(ieta)->SetDirectory(0);
        }
    }

    f->Close();
    delete f;
        
    return 1;
}

float JRTreader::GetRandomResponse(float pt, float eta, bool isBjet){
    

    int ptbin_orig = JRTreader::GetPtBin(pt);
    int etabin_orig = JRTreader::GetEtaBin(fabs(eta));

    int ptbin, etabin;
    JRTreader::GetModifiedBins(ptbin_orig, etabin_orig, isBjet, &ptbin, &etabin);

    TH1D *fit;
    if(isBjet)
        fit = fits_b->at(ptbin)->at(etabin);
    else
        fit = fits_nonb->at(ptbin)->at(etabin);

    if(fit==NULL){
        cout << "WARNING: could not find fit histogram for bins pt=" << ptbin << " eta=" << etabin << endl;
        return 1.;
    }

    float response = fit->GetRandom();
    return response;
}

float JRTreader::GetValue(float pt, float eta, bool isBjet, float smearfact){
    
    int ptbin_orig = JRTreader::GetPtBin(pt);
    int etabin_orig = JRTreader::GetEtaBin(fabs(eta));

    int ptbin, etabin;
    JRTreader::GetModifiedBins(ptbin_orig, etabin_orig, isBjet, &ptbin, &etabin);

    TH1D *fit;
    if(isBjet)
        fit = fits_b->at(ptbin)->at(etabin);
    else
        fit = fits_nonb->at(ptbin)->at(etabin);

    if(fit==NULL){
        cout << "WARNING: could not find fit histogram for bins pt=" << ptbin << " eta=" << etabin << endl;
        return 1.;
    }

    float response = fit->GetBinContent(fit->FindBin(smearfact));
    return response;
}

void JRTreader::UseRawHistograms(bool use){
    useFits = !use;
}

void JRTreader::GetModifiedBins(int ptbin, int etabin, bool isBjet, int *new_ptbin, int *new_etabin){

    if(isBjet){
        if(ptbin > 19) ptbin = 19;

        if(ptbin>=0 && ptbin<=1)
            ;
        if(ptbin>=2 && ptbin<=4)
            if(etabin >= 11) etabin = 10;
        if(ptbin>=5 && ptbin<=7)
            if(etabin >= 10) etabin = 9;
        if(ptbin>=8 && ptbin<=9)
            if(etabin >=9) etabin = 8;
        if(ptbin>=10 && ptbin<=12)
            if(etabin >= 8) etabin = 7;
        if(ptbin>=13 && ptbin<=15)
            if(etabin >= 7) etabin = 6;
        if(ptbin>=16 && ptbin<=17)
            if(etabin >= 5) etabin = 4;
        if(ptbin>=18 && ptbin<=19)
            if(etabin >=4) etabin = 3;

    }
    if(!isBjet){
        if(ptbin > 19) ptbin = 19;

        if(ptbin>=0 && ptbin<=3)
            ;
        if(ptbin>=4 && ptbin<=5)
            if(etabin >= 11) etabin = 10;
        if(ptbin>=6 && ptbin<=8)
            if(etabin >= 10) etabin = 9;
        if(ptbin>=9 && ptbin<=10)
            if(etabin >=9) etabin = 8;
        if(ptbin>=11 && ptbin<=13)
            if(etabin >= 8) etabin = 7;
        if(ptbin>=14 && ptbin<=15)
            if(etabin >= 7) etabin = 6;
        if(ptbin>=16 && ptbin<=17)
            if(etabin >= 6) etabin = 5;
        if(ptbin>=18 && ptbin<=18)
            if(etabin >=5) etabin = 4;
        if(ptbin>=19 && ptbin<=19)
            if(etabin >=4) etabin = 3;

    }

    *new_ptbin = ptbin;
    *new_etabin = etabin;
    
}

int JRTreader::GetPtBin(float pt){
    float pt_bins[] = {0,20,30,50,80,120,170,230,300,380,470,570,680,800,1000,1300,1700,2200,2800,3500,4300,5200,6500};
    int n_pt_bins = 23;

    for(int i=0; i<n_pt_bins-1; i++){
        if(pt >= pt_bins[i] && pt < pt_bins[i+1])
            return i;
    }

    return n_pt_bins-1;
}

int JRTreader::GetEtaBin(float eta){
    float eta_bins[] = {0,0.3,0.5,0.8,1.1,1.4,1.7,2.3,2.8,3.2,4.1,5.0,-1};
    int n_eta_bins = 12;

    for(int i=0; i<n_eta_bins-1; i++){
        if(eta >= eta_bins[i] && eta < eta_bins[i+1])
            return i;
    }

    return n_eta_bins-1;
}

