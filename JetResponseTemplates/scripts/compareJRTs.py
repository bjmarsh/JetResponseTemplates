#! /usr/bin/env python
import os
import ROOT
import pyRootPlotMaker as ppm

ROOT.gROOT.SetBatch(1)

# f1 = ROOT.TFile("JRTs/current.root")
# f2 = ROOT.TFile("JRTs/old.root")

# f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates.root")
# f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_80x_v2.root")
# f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_usedByJason.root")
# f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_2_8/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_92x_TGfit.root")
# f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_2_8/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_92x_smooth.root")
# f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_4_1/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x.root")
# f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_4_1/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_fixMatching3.root")
# f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_4_1/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_fixMatching4.root")
f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_4_7/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_JetID_PUID_BTagSFs_core2sigma.root")
f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_10_2_5/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_102x_JetID_PUID_BTagSFs_core2sigma.root")
# f2 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_80x_JetID_PUID_BTagSFs_core2sigma.root")
# f1 = ROOT.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_9_4_7/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_fixMatching5_ecalDeadCell.root")

outdir = "~/public_html/JRTplots/compare_fits_ptBinned_94x_vs_102x_JetID_PUID_BTagSFs_core2sigma"

os.system("mkdir -p {0}/nonbjets".format(outdir))
os.system("mkdir -p {0}/bjets".format(outdir))
os.system("cp ~/scripts/index.php {0}/nonbjets".format(outdir))
os.system("cp ~/scripts/index.php {0}/bjets".format(outdir))

n_pt_bins = 22
n_eta_bins = 12

def ZeroErrors(h):
    for i in range(h.GetNbinsX()):
        h.SetBinError(i,0)

for pt_bin in range(n_pt_bins):
    for eta_bin in range(n_eta_bins):

        # b jets
        # name1 = "h_b_JetAll_ResponsePt_Pt{0}_Eta{1}".format(pt_bin, eta_bin)
        name1 = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_bjet".format(pt_bin, eta_bin)
        name2 = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_bjet".format(pt_bin, eta_bin)
        found = True
        try:
            h1 = f1.Get(name2)
            h2 = f2.Get(name2)
        except:
            found = False
        # h1.Sumw2()
        # h2.Sumw2()
        if found and type(h1)!=ROOT.TObject and type(h2)!=ROOT.TObject:
            # ZeroErrors(h1)
            # ZeroErrors(h2)
            saveAs = outdir+"/bjets/pt{0:02d}_eta{1:02d}_bjets".format(pt_bin,eta_bin)
            for ext in [".pdf",".png"]:
                ppm.plotComparison(h1, h2, ratioTitle="New/Old", h1Title="Old", h2Title="New",
                                   saveAs=saveAs+ext, isLog=True, normalize=True, style=1, doOverflow=False,
                                   xAxisTitle="p_{T}^{reco}/p_{T}^{gen}", yRangeUser=(1e-5,10))

        # nonb jets
        name1 = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_nonbjet".format(pt_bin, eta_bin)
        name2 = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_nonbjet".format(pt_bin, eta_bin)
        found = True
        try:
            h1 = f1.Get(name2).Clone()
            h2 = f2.Get(name2).Clone()
        except:
            found = False
        # h1.Sumw2()
        # h2.Sumw2()
        if found and type(h1)!=ROOT.TObject and type(h2)!=ROOT.TObject:
            # ZeroErrors(h1)
            # ZeroErrors(h2)
            saveAs = outdir+"/nonbjets/pt{0:02d}_eta{1:02d}_nonbjets".format(pt_bin,eta_bin)
            for ext in [".pdf",".png"]:
                ppm.plotComparison(h1, h2, ratioTitle="New/Old", h1Title="Old", h2Title="New",
                                   saveAs=saveAs+ext, isLog=True, normalize=True, style=1, doOverflow=False,
                                   xAxisTitle="p_{T}^{reco}/p_{T}^{gen}", yRangeUser=(1e-5,10))

        # # all jets
        # name1 = "h_tot_JetAll_ResponsePt_Pt{0}_Eta{1}".format(pt_bin, eta_bin)
        # name2 = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_alljet".format(pt_bin, eta_bin)
        # h1 = f1.Get(name2)
        # h2 = f2.Get(name2)
        # h2.Sumw2()

        # saveAs = "JRTplots/compare201/justPT/alljets/pt{0:02d}_eta{1:02d}_alljets".format(pt_bin,eta_bin)
        # for ext in [".pdf",".png"]:
        #     ppm.plotComparison(h1, h2, ratioTitle="2016/2015", h1Title="2015", h2Title="2016",
        #                        saveAs=saveAs+ext, isLog=False, normalize=True, 
        #                        xAxisTitle="p_{T}^{reco}/p_{T}^{gen}")
