#include <vector>
#include <iostream>
#include <algorithm>

#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTreeCache.h>
#include <TString.h>
#include <TDirectory.h>
#include <TLorentzVector.h>

#include "JetResponseTemplates/JetResponseTemplates/interface/JRTreader.h"
#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"
#include "JetResponseTemplates/JetResponseTemplates/interface/JRTrebalTree.h"

// #include "/home/users/bemarsh/CORE/Tools/hemJet.h"
// #include "/home/users/bemarsh/CORE/Tools/MT2/MT2.h"

#define NSMEARS 100
#define EVTRED 1

using namespace std;

void computeJetVars(string jettype, float met_pt, float met_phi, float *results);

float metCut = 0;
bool doGenJetMatching = false;
bool doRecoJetMatching = false;

JRTTree t;
vector<float> smearjet_pt;
vector<float> smearjet_eta;
vector<float> smearjet_phi;
vector<float> rebaljet_pt;
vector<float> rebaljet_eta;
vector<float> rebaljet_phi;
vector<float> RSjet_pt;
vector<float> RSjet_eta;
vector<float> RSjet_phi;

bool doRebalance = false;
JRTrebalTree rebalTree;
string rebalanceDir = "/home/users/bemarsh/analysis/mt2/current/MT2Analysis/RebalanceAndSmear/output/test";

int main(int argc, char* argv[]) 
{

    if(argc<5){
        cout << "USAGE: JRTsmeartest <template file> <input dir> <output dir> <sample>" << endl;
        return 1;
    }

    //setup chain
    TChain *chain = new TChain("Events");
    TChain *rebalChain = new TChain("rebalance");

    chain->Add(Form("%s/%s*.root",argv[2],argv[4]));
    if(doRebalance)
        rebalChain->Add(Form("%s/%s*.root",rebalanceDir.c_str(),argv[4]));

    // open output file and initialize histograms
    TFile *fout = new TFile(Form("%s/%s.root",argv[3],argv[4]), "RECREATE");

    string tags[] = {"reco","gen","gensmear","rebalance"};
    vector<TH1D*> h_nJet30_vec;
    vector<TH1D*> h_met_vec;
    vector<TH1D*> h_ht_vec;
    vector<TH1D*> h_mht_vec;
    vector<TH1D*> h_diffMetMhtOverMet_vec;
    vector<TH1D*> h_j1pt_vec;
    vector<TH1D*> h_j2pt_vec;
    vector<TH1D*> h_deltaPhiMin_vec;
    for(int i=0; i<4; i++){
        h_nJet30_vec.push_back(new TH1D(("h_nJet30_"+tags[i]).c_str(),";nJet30",10,0,10));
        h_met_vec.push_back(new TH1D(("h_met_"+tags[i]).c_str(),";E_{T}^{miss} [GeV]",150,0,1500));
        h_ht_vec.push_back(new TH1D(("h_ht_"+tags[i]).c_str(), ";H_{T} [GeV]",150,0,1500));
        h_mht_vec.push_back(new TH1D(("h_mht_"+tags[i]).c_str(), ";H_{T}^{miss} [GeV]",150,0,1500));
        h_diffMetMhtOverMet_vec.push_back(new TH1D(("h_diffMetMhtOverMet_"+tags[i]).c_str(), ";diffMetMhtOverMet",150,0,10));
        h_j1pt_vec.push_back(new TH1D(("h_j1pt_"+tags[i]).c_str(), ";p_{T}(jet1) [GeV]",150,0,1500));
        h_j2pt_vec.push_back(new TH1D(("h_j2pt_"+tags[i]).c_str(), ";p_{T}(jet2) [GeV]",150,0,1500));
        h_deltaPhiMin_vec.push_back(new TH1D(("h_deltaPhiMin_"+tags[i]).c_str(), ";deltaPhiMin",150,0,3.141593));
    }
    
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    TFile *rebalFile = 0;
    if(doRebalance){
        // FOR NOW: expecting only 1 rebalance file. Kinda hackish. Oh well.
        rebalFile = new TFile(((TFile*)TIter(rebalChain->GetListOfFiles()).Next())->GetTitle());
        if(rebalFile->IsZombie()){
            cout << "Couldn't find rebalance file!" << endl;
            return 1;
        } 
        rebalTree.Init((TTree*)rebalFile->Get("rebalance"));
    }

    cout << "Processing " << nEventsChain << " events." << endl;

    // reader to get random values from templates.
    JRTreader reader;
    reader.UseRawHistograms(false);
    reader.Init(argv[1]);

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
        
        int nBadEvents = 0;

        unsigned int nEventsTree = tree->GetEntriesFast();
        cout << " Events in tree: " << nEventsTree << endl;
        unsigned int maxEvents = nEventsTree;
        maxEvents = nEventsTree;
        for( unsigned int event = 0; event < maxEvents; ++event) {
    
            if(event%EVTRED != 0)
                continue;

            // Get Event Content
            Long64_t tentry = tree->LoadTree(event);
            t.GetEntry(tentry);

            if(doRebalance)
                rebalTree.GetEntry(nEventsTotal);

            ++nEventsTotal;

            scale *= EVTRED;

            //ANALYSIS CODE HERE
            
            // variables to store the met after subtracting recojets.
            // used to recompute the met with the smeared jets
            float met_sub_x = t.pfmet_pt * cos(t.pfmet_phi);
            float met_sub_y = t.pfmet_pt * sin(t.pfmet_phi);
            
            float jetVars[7]; //will hold ht, mht, diffMetMht/met

            // subtract reco jets from met to recompute using smeared jets
            for(int ir=0; ir<t.n_recojet; ir++){
                float pt = t.recojet_pt->at(ir);
                float phi = t.recojet_phi->at(ir);

                // subtract reco jets from met
                if(pt > 10){
                    met_sub_x += pt * cos(phi);
                    met_sub_y += pt * sin(phi);
                }
            }// recojet loop

            // COMPUTE RECOJET QUANTITIES
            computeJetVars("reco", t.pfmet_pt, t.pfmet_phi, jetVars);
            // computeJetVars(t.recojet_pt, t.recojet_eta, t.recojet_phi, NULL, t.pfmet_pt, t.pfmet_phi, jetVars);
            h_nJet30_vec.at(0)->Fill(jetVars[0], scale);
            if(jetVars[0] >= 2 && t.pfmet_pt>metCut){
                h_met_vec.at(0)->Fill(t.pfmet_pt, scale);
                h_ht_vec.at(0)->Fill(jetVars[1], scale);
                h_mht_vec.at(0)->Fill(jetVars[2], scale);
                h_diffMetMhtOverMet_vec.at(0)->Fill(jetVars[3], scale);
                h_j1pt_vec.at(0)->Fill(jetVars[4], scale);
                h_j2pt_vec.at(0)->Fill(jetVars[5], scale);
                h_deltaPhiMin_vec.at(0)->Fill(jetVars[6], scale);
            }

            // COMPUTE GENJET QUANTITIES //
            computeJetVars("gen", t.genmet_pt, t.genmet_phi, jetVars);
            h_nJet30_vec.at(1)->Fill(jetVars[0], scale);
            if(jetVars[0] >= 2){
                h_met_vec.at(1)->Fill(t.genmet_pt, scale);
                h_ht_vec.at(1)->Fill(jetVars[1], scale);
                h_mht_vec.at(1)->Fill(jetVars[2], scale);
                h_diffMetMhtOverMet_vec.at(1)->Fill(jetVars[3], scale);
                h_j1pt_vec.at(1)->Fill(jetVars[4], scale);
                h_j2pt_vec.at(1)->Fill(jetVars[5], scale);
                h_deltaPhiMin_vec.at(1)->Fill(jetVars[6], scale);
            }

            // SMEAR JETS AND COMPUTE QUANTITIES
            for(int ismear=0; ismear < NSMEARS; ismear++){
                smearjet_pt.clear();
                smearjet_eta.clear();
                smearjet_phi.clear();

                // weird corrupt babies causing crashes
                if(t.genjet_flavour_cmssw->size() != t.genjet_pt->size()){
                    nBadEvents++;
                    break;
                }

                float met_smear_x = t.genmet_pt * cos(t.genmet_phi);
                float met_smear_y = t.genmet_pt * sin(t.genmet_phi);

                for(int ig=0; ig<t.n_genjet; ig++){
                    int flavour = t.genjet_flavour_cmssw->at(ig);

                    //either a leptonic or unidentifiable jet
                    if(flavour < 1){
                        continue;
                        // flavour = 1;
                    }

                    bool isBjet = 0;
                    if(flavour==5)
                        isBjet = 1;

                    // get the smearing factor
                    float smearfact = reader.GetRandomResponse(t.genjet_pt->at(ig), t.genjet_eta->at(ig), isBjet);
                    float smearpt = t.genjet_pt->at(ig) * smearfact;
                    
                    //add the smeared jet to the list
                    smearjet_pt.push_back(smearpt);
                    smearjet_eta.push_back(t.genjet_eta->at(ig));
                    smearjet_phi.push_back(t.genjet_phi->at(ig));

                    met_smear_x += (t.genjet_pt->at(ig)-smearpt) * cos(t.genjet_phi->at(ig));
                    met_smear_y += (t.genjet_pt->at(ig)-smearpt) * sin(t.genjet_phi->at(ig));
                    
                }//genjet loop

 
                // //recompute met with smeared jets
                // float met_smear_x = met_sub_x;
                // float met_smear_y = met_sub_y;
                // for(unsigned int ij=0; ij<smearjet_pt.size(); ij++){
                //     float pt = smearjet_pt.at(ij);
                //     float phi = smearjet_phi.at(ij);

                //     if(pt > 10){
                //         met_smear_x -= pt * cos(phi);
                //         met_smear_y -= pt * sin(phi);
                //     }
                    
                // }// loop over smeared jets

                
                float met_smear = sqrt(met_smear_x*met_smear_x + met_smear_y*met_smear_y);
                float met_smear_phi = atan2(met_smear_y, met_smear_x);
                // COMPUTE SMEARJET QUANTITIES //
                computeJetVars("smear", met_smear, met_smear_phi, jetVars);
                h_nJet30_vec.at(2)->Fill(jetVars[0], scale/NSMEARS);
                // if(jetVars[0] >= 2 && met_smear>metCut){
                if(jetVars[0] >= 2 && met_smear>30 && jetVars[1]>1000){
                    h_met_vec.at(2)->Fill(met_smear, scale/NSMEARS);
                    h_ht_vec.at(2)->Fill(jetVars[1], scale/NSMEARS);
                    h_mht_vec.at(2)->Fill(jetVars[2], scale/NSMEARS);
                    h_diffMetMhtOverMet_vec.at(2)->Fill(jetVars[3], scale/NSMEARS);
                    h_j1pt_vec.at(2)->Fill(jetVars[4], scale/NSMEARS);
                    h_j2pt_vec.at(2)->Fill(jetVars[5], scale/NSMEARS);
                    h_deltaPhiMin_vec.at(2)->Fill(jetVars[6], scale/NSMEARS);
                }

            }// smearing loop


            if(!doRebalance)
                continue;

            rebaljet_pt.clear();
            rebaljet_eta.clear();
            rebaljet_phi.clear();

            if((unsigned int)t.n_recojet != rebalTree.useJet->size()){
                cout << "WARNING: mismatched event! " << event << " " << nEventsTotal-1 << endl;
                continue;
            }

            float rebal_met_pt = t.pfmet_pt;
            float rebal_met_phi = t.pfmet_phi;
            float rebal_met_x = rebal_met_pt * cos(rebal_met_phi);
            float rebal_met_y = rebal_met_pt * sin(rebal_met_phi);

            int ireb = -1;
            for(int i=0; i<t.n_recojet; i++){
                if(rebalTree.useJet->at(i) != 1)
                    continue;
                ireb += 1;

                float newpt = t.recojet_pt->at(i) * rebalTree.rebalanceFactors->at(ireb);

                rebal_met_x += (t.recojet_pt->at(i) - newpt) * cos(t.recojet_phi->at(i));
                rebal_met_y += (t.recojet_pt->at(i) - newpt) * sin(t.recojet_phi->at(i));

                rebaljet_pt.push_back(newpt);
                rebaljet_eta.push_back(t.recojet_eta->at(i));
                rebaljet_phi.push_back(t.recojet_phi->at(i));
            }
            rebal_met_pt = sqrt(rebal_met_x*rebal_met_x + rebal_met_y*rebal_met_y);
            rebal_met_phi = atan2(rebal_met_y, rebal_met_x);

            // COMPUTE REBALANCED JET QUANTITIES //
            computeJetVars("rebalance", rebal_met_pt, rebal_met_phi, jetVars);
            h_nJet30_vec.at(3)->Fill(jetVars[0], scale);
            if(jetVars[0] >= 2){
                h_met_vec.at(3)->Fill(rebal_met_pt, scale);
                h_ht_vec.at(3)->Fill(jetVars[1], scale);
                h_mht_vec.at(3)->Fill(jetVars[2], scale);
                h_diffMetMhtOverMet_vec.at(3)->Fill(jetVars[3], scale);
                h_j1pt_vec.at(3)->Fill(jetVars[4], scale);
                h_j2pt_vec.at(3)->Fill(jetVars[5], scale);
                h_deltaPhiMin_vec.at(3)->Fill(jetVars[6], scale);
            }

        }//event loop
        
        cout << " nBadEvents: " << nBadEvents << endl;

    }//file loop

    fout->cd();
    for(int i=0; i<4; i++){
        TDirectory *d = fout->mkdir(tags[i].c_str());
        d->cd();
        h_nJet30_vec.at(i)->Write();
        h_met_vec.at(i)->Write();
        h_ht_vec.at(i)->Write();
        h_mht_vec.at(i)->Write();
        h_diffMetMhtOverMet_vec.at(i)->Write();
        h_j1pt_vec.at(i)->Write();
        h_j2pt_vec.at(i)->Write();
        h_deltaPhiMin_vec.at(i)->Write();
    }



    fout->Close();

}

// return indices of jets passing pt>30, eta<2.5
vector<int> getGoodJets(vector<float> &pt, vector<float> &eta){
    
    vector<int> ind;
    for(unsigned int i=0; i<pt.size(); i++){
        if(pt.at(i)>30 && eta.at(i)<2.5)
            ind.push_back(i);
    }
    return ind;

}

bool pairCompare(pair<float, float> p1, pair<float,float> p2){
    return p1.first > p2.first;
}

float DeltaPhi(float phi1, float phi2){
    float dphi = fabs(phi1-phi2);
    if(dphi > 3.14159265){
        dphi = 6.28318530 - dphi;
    }
    return dphi;
}

float DeltaR(float eta1, float phi1, float eta2, float phi2){
    float dphi = fabs(phi1-phi2);
    const float PI = 3.14159265359;
    if(dphi > PI)
        dphi = 2*PI - dphi;
    return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
}

float overlapsGenJet(float eta, float phi){
    for(int i=0; i<t.n_genjet; i++){
        if(DeltaR(eta, phi, t.genjet_eta->at(i), t.genjet_phi->at(i)) < 0.3)
            return true;
    }
    return false;
}

float overlapsRecoJet(float eta, float phi){
    for(int i=0; i<t.n_recojet; i++){
        if(DeltaR(eta, phi, t.recojet_eta->at(i), t.recojet_phi->at(i)) < 0.3)
            if(t.recojet_pt->at(i) > 10)
                return true;
    }
    return false;
}

// compute nJet30, ht, mht, diffMetMht/met, jet1_pt, jet2_pt, deltaPhiMin
void computeJetVars(string jettype, float met_pt, float met_phi, float *results){
    vector<float> *pt = NULL;
    vector<float> *eta = NULL;
    vector<float> *phi = NULL;
    vector<int> *jetID = NULL;

    if(jettype=="reco"){
        pt = t.recojet_pt;
        eta = t.recojet_eta;
        phi = t.recojet_phi;
        jetID = t.recojet_isLoosePFJet;
    }
    if(jettype=="gen"){
        pt = t.genjet_pt;
        eta = t.genjet_eta;
        phi = t.genjet_phi;
        jetID = NULL;
    }
    if(jettype=="smear"){
        pt = &smearjet_pt;
        eta = &smearjet_eta;
        phi = &smearjet_phi;
        jetID = NULL;
    }
    if(jettype=="rebalance"){
        pt = &rebaljet_pt;
        eta = &rebaljet_eta;
        phi = &rebaljet_phi;
        jetID = NULL;
    }


    int njets = pt->size();

    float ht=0;
    float mht_x = 0;
    float mht_y = 0;
    
    float met_x = met_pt * cos(met_phi);
    float met_y = met_pt * sin(met_phi);

    vector< pair<float,float> > ptphi;
    int nJet30 = 0;
    int nJet30FailId = 0;

    for(int ij=0; ij<njets; ij++){
        if(pt->at(ij) > 30 && eta->at(ij) < 2.5){
            if(jettype=="reco" && jetID->at(ij)==0){
                nJet30FailId++;
                continue;                
            }

            if(jettype=="reco" && doGenJetMatching && !overlapsGenJet(eta->at(ij),phi->at(ij)))
                continue;

            if((jettype=="gen" || jettype=="smear") && doRecoJetMatching && !overlapsRecoJet(eta->at(ij),phi->at(ij)))
                continue;

            ht += pt->at(ij);
            mht_x -= pt->at(ij)*cos(phi->at(ij));
            mht_y -= pt->at(ij)*sin(phi->at(ij));
            ptphi.push_back(pair<float,float> (pt->at(ij), phi->at(ij)));
            nJet30++;


        }
    }

    float mht = sqrt(mht_x*mht_x + mht_y*mht_y);
    float diffMetMht = sqrt((met_x-mht_x)*(met_x-mht_x) + (met_y-mht_y)*(met_y*mht_y));
    
    sort(ptphi.begin(), ptphi.end(), pairCompare);

    results[0] = nJet30;
    results[1] = ht;
    results[2] = mht;
    results[3] = diffMetMht/met_pt;
    results[4] = 0;
    results[5] = 0;    
    if(nJet30 > 0)
        results[4] = ptphi.at(0).first;
    if(nJet30 > 1)
        results[5] = ptphi.at(1).first;

    results[6] = 999;
    for(int i=0; i<min(nJet30,4); i++){
        results[6] = min(results[6], DeltaPhi(met_phi, ptphi.at(i).second));
    }
    if(results[6]==999)
        results[6] = -1;


}
