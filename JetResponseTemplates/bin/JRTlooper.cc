#include <vector>
#include <utility>
#include <iostream>
#include <map>

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTreeCache.h>
#include <TString.h>
#include <TDirectory.h>

#include "DataFormats/Math/interface/deltaR.h"

#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"
#include "JetResponseTemplates/JetResponseTemplates/interface/BTagCalibrationStandalone.h"

using namespace std;

int getBin(float val, float *bins, int n_bins);
int getSampleID(const char* filename);
bool inDeadEcalCell(float eta, float phi);
bool passTightPFJetID_2017(float eta, char chf, char nhf, char cef, char nef, int cm, int nm, int npfcands);
double getBtagEffFromHist(float pt, float eta, BTagEntry::JetFlavor flavor, TH2D* effs);

std::map<string,float> scale1fbs;

// manual veto windows for dead ecal cells
const uint NDEADCELLS = 17;
float dead_cell_boundaries[NDEADCELLS][4] = {
    {-2.4,-1.6,-2.15,-1.3},
    {-2.15,-2.6,-1.95,-2.35},
    {-1.65,0.5,-1.5,0.7},
    {-1.35,-2.75,-1.25,-2.6},
    {-1.3,0.5,-1.2,0.6},
    {-1.15,0.5,-1.05,0.6},
    {-0.4,0.1,-0.2,0.3},
    {-1.0,2.9,-0.85,3.1},
    {0.85,2.8,1.0,2.9},
    {0.25,-0.8,0.35,-0.7},
    {0.6,-0.8,0.7,-0.7},
    {0.95,-0.8,1.05,-0.7},
    {1.35,-2.65,1.45,-2.5},
    {1.65,0.7,1.85,0.95},
    {1.65,-0.7,1.8,-0.5},
    {1.5,-1.6,1.65,-1.4},
    {1.7,-2.2,1.9,-1.95}
};


int main(int argc, char* argv[]) 
{
    //80x_v1
    scale1fbs["qcd_pt15to30"] = 165146.970476;
    scale1fbs["qcd_pt30to50"] = 14121.3721374;
    scale1fbs["qcd_pt50to80"] = 1952.64830653;
    scale1fbs["qcd_pt80to120"] = 395.396136109;
    scale1fbs["qcd_pt170to300"] = 16.853128483;
    scale1fbs["qcd_pt300to470"] = 0.348160902229;
    scale1fbs["qcd_pt470to600"] = 0.16368744738;
    scale1fbs["qcd_pt600to800"] = 0.0479672067533;
    scale1fbs["qcd_pt1000to1400"] = 0.00314040790659;
    scale1fbs["qcd_pt1400to1800"] = 0.00212570854849;
    scale1fbs["qcd_pt1800to2400"] = 5.79923291077e-05;
    scale1fbs["qcd_pt2400to3200"] = 1.71076282607e-05;
    scale1fbs["qcd_pt3200toInf"] = 4.22339081267e-07;

    if(argc<3){
        cout << "USAGE: JRTlooper <tag> <files to run over>" << endl;
        return 1;
    }

    // pt and eta binning
    float pt_bins[] = {0,20,30,50,80,120,170,230,300,380,470,570,680,800,1000,1300,1700,2200,2800,3500,4300,5200,6500,-1};
    int n_pt_bins = 23;
    float eta_bins[] = {0,0.3,0.5,0.8,1.1,1.4,1.7,2.3,2.8,3.2,4.1,5.0,-1};
    int n_eta_bins = 12;


    // setup histograms
    TH1D* h_pileup = new TH1D("h_pileup",";n vertices",100,0,100);
    TH1D* h_dr = new TH1D("h_dr",";dR(reco,gen)",100,0,0.5);
    
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


    BTagCalibration* calib = new BTagCalibration("csvv2", "btagsf/CSVv2_Moriond17_B_H.csv");
    BTagCalibrationReader* reader_btagsf = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
    reader_btagsf->load(*calib, BTagEntry::JetFlavor::FLAV_B, "comb");
    reader_btagsf->load(*calib, BTagEntry::JetFlavor::FLAV_C, "comb");
    reader_btagsf->load(*calib, BTagEntry::JetFlavor::FLAV_UDSG, "incl");
    TFile* f_btag_eff = new TFile("btagsf/btageff__ttbar_powheg_pythia8_25ns_Moriond17.root");
    TH2D* h_btag_eff_b = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_b")->Clone("h_btag_eff_b");
    TH2D* h_btag_eff_c = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_c")->Clone("h_btag_eff_c");
    TH2D* h_btag_eff_udsg = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_udsg")->Clone("h_btag_eff_udsg");
    h_btag_eff_b->SetDirectory(0);
    h_btag_eff_c->SetDirectory(0);
    h_btag_eff_udsg->SetDirectory(0);
    f_btag_eff->Close();

    TChain *chain = new TChain("Events");

    for(int i=2; i<argc; i++)
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

        float scale;

        // branch for scale1fb is not in JRTTree
        // tree->SetBranchAddress("evt_scale1fb", &scale);

        TString title (currentFile->GetTitle());
        int i1 = title.Index("qcd_pt");
        int i2 = title.Index("_",i1+5);
        if(i2==-1)
            i2 = title.Index(".",i1);
        TString samp(title(i1,i2-i1));
        // scale = scale1fbs[samp.Data()];
        scale = 1.0;
        cout << "    " << scale << endl;

        unsigned int nEventsTree = tree->GetEntriesFast();
        // nEventsTree = 10000.;
        for( unsigned int event = 0; event < nEventsTree; ++event) {
    
            // Get Event Content
            Long64_t tentry = tree->LoadTree(event);
            t.GetEntry(tentry);
            ++nEventsTotal;
            
            if(nEventsTotal%100000==0)
            // if(nEventsTotal%1==0)
                cout << "Processed " << nEventsTotal << " / " << nEventsChain << " events\n";

            if(!t.Flag_badMuonFilter2016 || !t.Flag_badChargedCandidateFilter2016)
                continue;
            if(!t.Flag_badPFMuonFilter || !t.Flag_badChargedCandidateFilter)
                continue;
            if(!t.Flag_ecalDeadCellTriggerPrimitiveFilter)
                continue;
            if(!t.Flag_hbheNoiseFilter)
                continue;
            if(!t.Flag_hbheNoiseIsoFilter)
                continue;
            if(!t.Flag_eeBadScFilter)
                continue;
            // if(!t.Flag_ecalBadCalibFilter)
            //     continue;
            if(!t.Flag_globalTightHalo2016Filter)
                continue;

            vector<pair<uint,uint> > matches;
            vector<uint> counts_gj (t.n_genjet, 0);
            vector<uint> counts_rj (t.n_recojet, 0);

            //ANALYSIS CODE HERE
            for(unsigned int igj=0; igj<t.genjet_pt->size(); igj++){
                // float gj_pt = t.genjet_pt->at(igj);
                float gj_eta = t.genjet_eta->at(igj);
                float gj_phi = t.genjet_phi->at(igj);
                
                for(unsigned int irj=0; irj<t.recojet_pt->size(); irj++){
                    float rj_pt = t.recojet_pt->at(irj);
                    float rj_eta = t.recojet_eta->at(irj);
                    float rj_phi = t.recojet_phi->at(irj);
                    
                    // reject recojets with pT < 10 GeV
                    if(rj_pt < 10) continue;
                    
                    float dr = deltaR(rj_eta, rj_phi, gj_eta, gj_phi);
                    if(dr < 0.3){
                        matches.push_back(pair<uint,uint>(igj,irj));
                    }
                    if(dr < 0.5){
                        counts_gj[igj]++;
                        counts_rj[irj]++;
                    }
                }//recojet loop

            } // genjet loop

            uint nmatches = matches.size();
            for(uint imatch=0; imatch<nmatches; imatch++){

                uint igj = matches[imatch].first;
                uint irj = matches[imatch].second;
                if(counts_gj[igj] != 1 || counts_rj[irj] != 1)
                    continue;
                
                float ratio = t.recojet_pt->at(irj) / t.genjet_pt->at(igj);
                int pt_bin = getBin(t.genjet_pt->at(igj), pt_bins, n_pt_bins);
                int eta_bin = getBin(fabs(t.genjet_eta->at(igj)), eta_bins, n_eta_bins);
                if(pt_bin==-1 || eta_bin==-1){
                    cout << "WARNING: bad genjet at event # " << event << endl;
                    continue;
                }

                // if(inDeadEcalCell(t.genjet_eta->at(igj), t.genjet_phi->at(igj)))
                //     continue;

                // apply jet ID
                // if(!isTightPFJet_2017(t.recojet_eta->at(irj), t.recojet_chFrac->at(irj), t.recojet_nhFrac->at(irj),
                //                       t.recojet_cemFrac->at(irj), t.recojet_nemFrac->at(irj),
                //                       t.recojet_chargedMult->at(irj), t.recojet_neutralMult->at(irj),
                //                       t.recojet_npfcands->at(irj)))
                //     continue;
                if(!t.recojet_isTightPFJet->at(irj))
                    continue;

                // apply pileup jet ID
                if(!(t.recojet_puId_ID->at(irj) & (1<<0))) //loose
                // if(!(t.recojet_puId_ID->at(irj) & (1<<1))) //medium
                // if(!(t.recojet_puId_ID->at(irj) & (1<<2))) //tight
                    continue;

                h_dr->Fill(deltaR(t.genjet_eta->at(igj), t.genjet_phi->at(igj), t.recojet_eta->at(irj), t.recojet_phi->at(irj)), scale);

                // if(abs(t.recojet_leadingPFCandId->at(irj)) == 13 && t.recojet_muFrac->at(irj) > 50 && t.genjet_muFrac->at(igj) == 0 && t.recojet_pt->at(irj)*0.01*t.recojet_muFrac->at(irj) > 100.0){
                //         cout << "[JRTlooper] matched muonic rj/non-muonic gj: " << t.evt_run << ":" << t.evt_lumi << ":" << t.evt_event << " " << ratio << " " << t.recojet_leadingPFCandId->at(irj) << " " << 
                //             t.Flag_badMuonFilter2016 << " " << t.Flag_badMuonFilter2016_loose << " " << t.Flag_badChargedCandidateFilter2016 << " " << 
                //             (int)t.recojet_muFrac->at(irj) << " " << t.genjet_flavour_cmssw->at(igj) << " " << t.genjet_flavour_bennett->at(igj) << endl;
                // }

                // if(pt_bin>=6 && eta_bin<=3 && ratio >= 2.0){
                //     cout << "[JRTlooper] high smear event: " << t.evt_run << ":" << t.evt_lumi << ":" << t.evt_event << " " << 
                //         t.Flag_badMuonFilter2016 << " " << t.Flag_badMuonFilter2016_loose << " " << t.Flag_badChargedCandidateFilter2016 << " " << ratio << " " << 
                //         t.recojet_leadingPFCandId->at(irj) << " " << (int)t.recojet_muFrac->at(irj) << " " << (int)t.genjet_muFrac->at(igj) << endl;
                // }

                if(pt_bin>=8 && eta_bin<=11 && ratio <= 0.5){
                    cout << "[JRTlooper]  low smear event: " << t.evt_run << ":" << t.evt_lumi << ":" << t.evt_event << " " << 
                        t.Flag_badMuonFilter2016 << " " << t.Flag_badMuonFilter2016_loose << " " << t.Flag_badChargedCandidateFilter2016 << " " << ratio << " " << 
                        t.recojet_leadingPFCandId->at(irj) << " " << (int)t.recojet_muFrac->at(irj) << " " << (int)t.genjet_muFrac->at(igj) << 
                        " " << t.genjet_pt->at(igj) << " " << t.genjet_eta->at(igj) << " " << t.genjet_phi->at(igj) << endl;
                }

                // // use gen-level flavor
                // if(t.genjet_flavour_cmssw->at(igj) == 5){
                //     JRTs_b->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                // }else if(t.genjet_flavour_cmssw->at(igj) >= 1){
                //     JRTs_l->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                // }                    
                // if(t.genjet_flavour_cmssw->at(igj) > 0)
                //     JRTs_a->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                
                // // use reco b-tagging
                // if(t.recojet_btagCSV->at(irj) > 0.8838)
                //     JRTs_b->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                // else
                //     JRTs_l->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                // JRTs_a->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);

                // use gen-flavor, but weight by efficiency corrected with btag SFs
                TH2D* eff_hist;
                BTagEntry::JetFlavor flavor;
                if(t.genjet_flavour_cmssw->at(igj) == 5){
                    eff_hist = h_btag_eff_b;
                    flavor = BTagEntry::FLAV_B;
                }else if(t.genjet_flavour_cmssw->at(igj) == 4){
                    eff_hist = h_btag_eff_c;
                    flavor = BTagEntry::FLAV_C;
                }else if(t.genjet_flavour_cmssw->at(igj) >= 1){
                    eff_hist = h_btag_eff_udsg;
                    flavor = BTagEntry::FLAV_UDSG;
                }
                if(t.genjet_flavour_cmssw->at(igj) >= 1){
                    double eff = getBtagEffFromHist(t.recojet_pt->at(irj), t.recojet_eta->at(irj), flavor, eff_hist);
                    float weight_cent = reader_btagsf->eval_auto_bounds("central", flavor, t.recojet_eta->at(irj), t.recojet_pt->at(irj));
                    float btagprob_data = weight_cent * eff;
                    JRTs_b->at(pt_bin)->at(eta_bin)->Fill(ratio, scale * btagprob_data);
                    JRTs_l->at(pt_bin)->at(eta_bin)->Fill(ratio, scale * (1-btagprob_data));
                    JRTs_a->at(pt_bin)->at(eta_bin)->Fill(ratio, scale);
                }


            }// match loop

            h_pileup->Fill(t.evt_nvertices, scale);

        }//event loop
    }//file loop

    // TFile *fout = new TFile("JetResponseTemplates.root","RECREATE");
    TFile *fout = new TFile(Form("looper_output/%s.root",argv[1]),"RECREATE");
    TDirectory *d_aux = fout->mkdir("auxiliary");
    d_aux->cd();
    h_pileup->Write();
    h_dr->Write();

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

bool inDeadEcalCell(float eta, float phi){

    for(uint i=0; i<NDEADCELLS; i++){
        if( eta > dead_cell_boundaries[i][0] &&
            eta < dead_cell_boundaries[i][2] &&
            phi > dead_cell_boundaries[i][1] &&
            phi < dead_cell_boundaries[i][3] )

            return true;
    }
    return false;

}

bool passTightPFJetID_2017(float eta, char chf, char nhf, char cef, char nef, int cm, int nm, int npfcands){
    if(fabs(eta) <= 2.4){
        if(chf == 0) return false;
        if(cm == 0) return false;
    }
    if(fabs(eta) <= 2.7){
        if(nef >= 90) return false;
        if(nhf >= 90) return false;
        if(npfcands <= 1) return false;
    }else if(fabs(eta) <= 3.0){
        if(nef <= 2 || nef >= 99) return false;
        if(nm <= 2) return false;
    }else{
        if(nhf <= 2) return false;
        if(nef >= 90) return false;
        if(nm <= 10) return false;
    }
    
    return true;
}

double getBtagEffFromHist(float pt, float eta, BTagEntry::JetFlavor flavor, TH2D* effs){
    float pt_cutoff;
    if(flavor == BTagEntry::FLAV_B)
        pt_cutoff = std::max(20.,std::min(599.,double(pt)));
    else
        pt_cutoff = std::max(20.,std::min(399.,double(pt)));

    int binx = effs->GetXaxis()->FindBin(pt_cutoff);
    int biny = effs->GetYaxis()->FindBin(fabs(eta));
    return effs->GetBinContent(binx, biny);
}
