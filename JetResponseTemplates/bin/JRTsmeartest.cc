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

#include "JetResponseTemplates/JetResponseTemplates/interface/JRTreader.h"
#include "JetResponseTemplates/JetResponseTemplates/interface/JRTTree.h"

#define NSMEARS 1

using namespace std;

void computeJetVars(vector<float> &pt, vector<float> &eta, vector<float> &phi, float met_pt, float met_phi, float *results);

int main(int argc, char* argv[]) 
{

    if(argc<3){
        cout << "USAGE: JRTsmeartest <template file> <files to run over>" << endl;
        return 1;
    }

    //setup chain
    TChain *chain = new TChain("Events");

    for(int i=2; i<argc; i++)
        chain->Add(argv[i]);

    // open output file and initialize histograms
    TFile *fout = new TFile("smeartest.root", "RECREATE");

    string tags[] = {"reco","gen","gensmear"};
    vector<TH1D*> h_met_vec;
    vector<TH1D*> h_ht_vec;
    vector<TH1D*> h_mht_vec;
    vector<TH1D*> h_diffMetMhtOverMet_vec;
    vector<TH1D*> h_j1pt_vec;
    vector<TH1D*> h_j2pt_vec;
    vector<TH1D*> h_deltaPhiMin_vec;
    for(int i=0; i<3; i++){
        h_met_vec.push_back(new TH1D(("h_met_"+tags[i]).c_str(),";E_{T}^{miss} [GeV]",150,0,1500));
        h_ht_vec.push_back(new TH1D(("h_ht_"+tags[i]).c_str(), ";H_{T} [GeV]",150,0,1500));
        h_mht_vec.push_back(new TH1D(("h_mht_"+tags[i]).c_str(), ";H_{T}^{miss} [GeV]",150,0,1500));
        h_diffMetMhtOverMet_vec.push_back(new TH1D(("h_diffMetMhtOverMet_"+tags[i]).c_str(), ";diffMetMhtOverMet",150,0,1500));
        h_j1pt_vec.push_back(new TH1D(("h_j1pt_"+tags[i]).c_str(), ";p_{T}(jet1) [GeV]",150,0,1500));
        h_j2pt_vec.push_back(new TH1D(("h_j2pt_"+tags[i]).c_str(), ";p_{T}(jet2) [GeV]",150,0,1500));
        h_deltaPhiMin_vec.push_back(new TH1D(("h_deltaPhiMin_"+tags[i]).c_str(), ";deltaPhiMin",150,0,3.141593));
    }
    
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    cout << "Processing " << nEventsChain << " events." << endl;

    JRTTree t;

    // reader to get random values from templates.
    JRTreader reader;
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

        unsigned int nEventsTree = tree->GetEntriesFast();
        unsigned int maxEvents = nEventsTree;
        // maxEvents = nEventsTree / 100;
        for( unsigned int event = 0; event < maxEvents; ++event) {
    
            if(event%100000 == 0)
                cout << "Processing event: " << event << endl;

            // Get Event Content
            Long64_t tentry = tree->LoadTree(event);
            t.GetEntry(tentry);
            ++nEventsTotal;

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
            computeJetVars(*t.recojet_pt, *t.recojet_eta, *t.recojet_phi, t.pfmet_pt, t.pfmet_phi, jetVars);
            if(jetVars[0] >= 2){
                h_met_vec.at(0)->Fill(t.pfmet_pt, scale);
                h_ht_vec.at(0)->Fill(jetVars[1], scale);
                h_mht_vec.at(0)->Fill(jetVars[2], scale);
                h_diffMetMhtOverMet_vec.at(0)->Fill(jetVars[3], scale);
                h_j1pt_vec.at(0)->Fill(jetVars[4], scale);
                h_j2pt_vec.at(0)->Fill(jetVars[5], scale);
                h_deltaPhiMin_vec.at(0)->Fill(jetVars[6], scale);
            }


            // COMPUTE GENJET QUANTITIES //
            computeJetVars(*t.genjet_pt, *t.genjet_eta, *t.genjet_phi, t.genmet_pt, t.genmet_phi, jetVars);
            if(jetVars[0] >= 2){
                h_met_vec.at(1)->Fill(t.genmet_pt, scale);
                h_ht_vec.at(1)->Fill(jetVars[1], scale);
                h_mht_vec.at(1)->Fill(jetVars[2], scale);
                h_diffMetMhtOverMet_vec.at(1)->Fill(jetVars[3], scale);
                h_j1pt_vec.at(1)->Fill(jetVars[4], scale);
                h_j2pt_vec.at(1)->Fill(jetVars[5], scale);
                h_deltaPhiMin_vec.at(1)->Fill(jetVars[6], scale);
            }

            //vectors to store kinematics of smeared jets
            vector<float> smearjet_pt;
            vector<float> smearjet_eta;
            vector<float> smearjet_phi;

            // SMEAR JETS AND COMPUTE QUANTITIES
            for(int ismear=0; ismear < NSMEARS; ismear++){
                smearjet_pt.clear();
                smearjet_eta.clear();
                smearjet_phi.clear();

                for(int ig=0; ig<t.n_genjet; ig++){
                    int flavour = t.genjet_flavour_bennett->at(ig);

                    //either a leptonic or unidentifiable jet
                    if(flavour < 1)
                        continue;

                    bool isBjet = 0;
                    if(flavour==5)
                        isBjet = 1;

                    // get the smearing factor
                    float smearfact = reader.GetResponse(t.genjet_pt->at(ig), t.genjet_eta->at(ig), isBjet);
                    float smearpt = t.genjet_pt->at(ig) * smearfact;
                    
                    //add the smeared jet to the list
                    smearjet_pt.push_back(smearpt);
                    smearjet_eta.push_back(t.genjet_eta->at(ig));
                    smearjet_phi.push_back(t.genjet_phi->at(ig));

                }//genjet loop

                //recompute met with smeared jets
                float met_smear_x = met_sub_x;
                float met_smear_y = met_sub_y;
                for(unsigned int ij=0; ij<smearjet_pt.size(); ij++){
                    float pt = smearjet_pt.at(ij);
                    float phi = smearjet_phi.at(ij);

                    if(pt > 10){
                        met_smear_x -= pt * cos(phi);
                        met_smear_y -= pt * sin(phi);
                    }
                    
                }// loop over smeared jets
                
                float met_smear = sqrt(met_smear_x*met_smear_x + met_smear_y*met_smear_y);
                float met_smear_phi = atan2(met_smear_y, met_smear_x);
                // COMPUTE GENJET QUANTITIES //
                computeJetVars(smearjet_pt, smearjet_eta, smearjet_phi, met_smear, met_smear_phi, jetVars);
                if(jetVars[0] >= 2){
                    h_met_vec.at(2)->Fill(met_smear, scale);
                    h_ht_vec.at(2)->Fill(jetVars[1], scale);
                    h_mht_vec.at(2)->Fill(jetVars[2], scale);
                    h_diffMetMhtOverMet_vec.at(2)->Fill(jetVars[3], scale);
                    h_j1pt_vec.at(2)->Fill(jetVars[4], scale);
                    h_j2pt_vec.at(2)->Fill(jetVars[5], scale);
                    h_deltaPhiMin_vec.at(2)->Fill(jetVars[6], scale);
                }

            }// smearing loop


        }//event loop
    }//file loop

    fout->cd();
    for(int i=0; i<3; i++){
        TDirectory *d = fout->mkdir(tags[i].c_str());
        d->cd();
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

// compute nJet30, ht, mht, diffMetMht/met, jet1_pt, jet2_pt, deltaPhiMin
void computeJetVars(vector<float> &pt, vector<float> &eta, vector<float> &phi, float met_pt, float met_phi, float *results){
    int njets = pt.size();

    float ht=0;
    float mht_x = 0;
    float mht_y = 0;
    
    float met_x = met_pt * cos(met_phi);
    float met_y = met_pt * sin(met_phi);

    vector< pair<float,float> > ptphi;
    int nJet30 = 0;

    for(int ij=0; ij<njets; ij++){
        if(pt.at(ij) > 30 && eta.at(ij) < 2.5){
            ht += pt.at(ij);
            mht_x -= pt.at(ij)*cos(phi.at(ij));
            mht_y -= pt.at(ij)*sin(phi.at(ij));
            ptphi.push_back(pair<float,float> (pt.at(ij), phi.at(ij)));
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
    if(nJet30 > 0)
        results[4] = ptphi.at(0).first;
    if(nJet30 > 1)
        results[5] = ptphi.at(1).first;

    results[6] = 999;
    for(int i=0; i<min(nJet30,4); i++){
        results[6] = min(results[6], DeltaPhi(met_phi, ptphi.at(i).second));
    }

}
