#! /usr/bin/env python

import os
import ROOT as r
import pyRootPlotMaker as ppm

r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(1)

# fid = r.TFile("~/analysis/mt2/current/MT2Analysis/RebalanceAndSmear/data/JetResponseTemplates.root")
# fid = r.TFile("~/analysis/mt2/current/MT2Analysis/babymaker/rebal/JetResponseTemplates.root")
# fid = r.TFile("/home/users/bemarsh/analysis/jet_response_templates/CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_usedByJason.root")
# fid = r.TFile("~/analysis/jet_response_templates/CMSSW_9_4_1/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_Bennett.root")
# fid = r.TFile("~/analysis/jet_response_templates/CMSSW_9_4_7/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_JECV32_JetID_PUID_BTagSFs_core2sigma.root")
fid = r.TFile("~/analysis/jet_response_templates/CMSSW_10_2_5/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_102x_JetID_PUID_BTagSFs_core2sigma.root")
# fid = r.TFile("~/analysis/jet_response_templates/CMSSW_9_4_7/src/JetResponseTemplates/JetResponseTemplates/output/JetResponseTemplates_ptBinned_94x_fixMatching5_TG.root")

# outdir = "~/public_html/JRTplots/fits_ptBinned_94x_fixMatching5_TG"
# outdir = "~/public_html/JRTplots/fits_ptBinned_94x_JECV32_JetID_PUID_BTagSFs_core2sigma"
outdir = "~/public_html/JRTplots/fits_ptBinned_102x_JetID_PUID_BTagSFs_core2sigma"
# outdir = "~/public_html/JRTplots/fits_usedByJason"

n_pt_bins = 22
# n_pt_bins = 5
n_eta_bins = 12

pt_bins = [10,20,30,50,80,120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500]
eta_bins = [0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0, "#infty"]

os.system("mkdir -p {0}/nonbjets_lin".format(outdir))
os.system("mkdir -p {0}/nonbjets_log".format(outdir))
os.system("mkdir -p {0}/bjets_lin".format(outdir))
os.system("mkdir -p {0}/bjets_log".format(outdir))
os.system("cp ~/scripts/index.php {0}/nonbjets_lin".format(outdir))
os.system("cp ~/scripts/index.php {0}/nonbjets_log".format(outdir))
os.system("cp ~/scripts/index.php {0}/bjets_lin".format(outdir))
os.system("cp ~/scripts/index.php {0}/bjets_log".format(outdir))

def DrawJRT(h_jrt, f_core, f_tail, f_all, outdir, pt_bin, eta_bin, isB=False, doLog=False):
    h = h_jrt.Clone()
    c = r.TCanvas()
    c.SetTickx()
    c.SetTicky()
    c.SetCanvasSize(800,700)
    c.SetTopMargin(0.06)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.08)
    c.SetRightMargin(0.05)
    c.SetLogy(doLog)
    h.SetTitle("")
    h.GetXaxis().SetTitle("p_{T}^{reco} / p_{T}^{gen}")
    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitleSize(0.04)
    h.GetXaxis().SetTitleOffset(1.20)
    h.GetYaxis().SetTitle("")
    h.GetYaxis().SetLabelSize(0.04)
    h.GetYaxis().SetLabelOffset(0.01)
    if doLog:
        h.SetMinimum(1e-5)
        h.SetMaximum(10)
    else:
        h.SetMaximum(h.GetMaximum()*1.1)
    h.Draw("PE")
    if type(f1)==type(r.TH1D()):
        f_all.SetLineColor(r.kGray+1)
        f_all.SetLineWidth(2)
        f_all.Draw("SAME L HIST")
        f_core.SetLineColor(r.kRed)
        f_core.SetLineWidth(2)
        f_core.Draw("SAME L HIST")
        f_tail.SetLineColor(r.kGreen+2)
        f_tail.SetLineWidth(2)
        f_tail.Draw("SAME L HIST")
    h.Draw("PE SAME")            
    text = r.TLatex()
    text.SetNDC(1)
    text.SetTextFont(62)
    text.SetTextColor(r.kBlue)
    text.SetTextSize(0.05)
    text.SetTextAlign(11)
    x = 0.58 if pt_bin <= 3 else 0.52 if pt_bin <= 12 else 0.48
    y = (0.87 if pt_bin<=1 else 0.82 if pt_bin<=3 else 0.72) if doLog else 0.72
    text.DrawLatex(x, y, "{0}b jets".format("" if isB else "non-"))
    text.DrawLatex(x, y-0.07, "{0} #leq p_{{T}} < {1} GeV".format(pt_bins[pt_bin], pt_bins[pt_bin+1]))
    text.DrawLatex(x, y-0.14, "{0} #leq |#eta| < {1}".format(eta_bins[eta_bin], eta_bins[eta_bin+1]))
    text.SetTextColor(r.kBlack)
    text.SetTextSize(0.045)
    text.DrawLatex(0.09, 0.95, "CMS Simulation")
    text.SetTextFont(42)
    text.SetTextAlign(31)
    year = 2018 if "102x" in outdir else 2017 if "94x" in outdir else 2016
    text.DrawLatex(0.94, 0.95, "13 TeV ({0})".format(year))

    saveAs = os.path.join(outdir, "pt{0:02d}_eta{1:02d}_{2}bjets".format(pt_bin,eta_bin,"" if isB else "non"))
    c.SaveAs(saveAs+".pdf")
    c.SaveAs(saveAs+".png")


for pt_bin in range(n_pt_bins):
    for eta_bin in range(n_eta_bins):

        # nonb jets
        hname = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_nonbjet".format(pt_bin, eta_bin)
        h = fid.Get(hname)
        fname1 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_nonbjet_core".format(pt_bin, eta_bin)
        fname2 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_nonbjet_tail".format(pt_bin, eta_bin)
        fname3 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_nonbjet".format(pt_bin, eta_bin)
        f1 = fid.Get(fname1)
        f2 = fid.Get(fname2)
        f3 = fid.Get(fname3)
        DrawJRT(h, f1, f2, f3, os.path.join(outdir, "nonbjets_lin"), pt_bin, eta_bin, isB=False, doLog=False)
        DrawJRT(h, f1, f2, f3, os.path.join(outdir, "nonbjets_log"), pt_bin, eta_bin, isB=False, doLog=True)


        # b jets
        hname = "pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_bjet".format(pt_bin, eta_bin)
        h = fid.Get(hname)
        fname1 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_bjet_core".format(pt_bin, eta_bin)
        fname2 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_bjet_tail".format(pt_bin, eta_bin)
        fname3 = "pt{0}/pt{0}_eta{1}/fit_pt{0}_eta{1}_bjet".format(pt_bin, eta_bin)
        f1 = fid.Get(fname1)
        f2 = fid.Get(fname2)
        f3 = fid.Get(fname3)
        DrawJRT(h, f1, f2, f3, os.path.join(outdir, "bjets_lin"), pt_bin, eta_bin, isB=True, doLog=False)
        DrawJRT(h, f1, f2, f3, os.path.join(outdir, "bjets_log"), pt_bin, eta_bin, isB=True, doLog=True)
        

