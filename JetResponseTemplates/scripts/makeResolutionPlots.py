import numpy as np
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)

jrt_file = r.TFile("../output/JetResponseTemplates_ptBinned_102x_JetID_PUID_BTagSFs_core2sigma.root")

pt_bins = np.array([10, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300], dtype=float)
eta_bins = [0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.3, 2.8, 3.2, 4.1, 5.0]

def get_mean_res(h_jrt):
    f = r.TF1("fit","gaus",0,3)
    max = h_jrt.GetBinContent(h_jrt.GetMaximumBin())
    mean = h_jrt.GetBinCenter(h_jrt.GetMaximumBin())
    rms = h_jrt.GetRMS()
    f.SetRange(mean-rms, mean+rms)
    f.SetParameters(max, mean, min(rms,0.1))
    # f.SetParLimits(0, max-0.003, max+0.003)
    f.SetParLimits(1, mean-0.02, mean+0.02)
    h_jrt.Fit(f, "QNR", "goff")
    return f.GetParameter(1), f.GetParameter(2)

def GetResPlot(eta_bin, isB, gr):
    ieta = 0
    for ipt in range(len(pt_bins)-1):
        jrt = jrt_file.Get("pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_{2}bjet".format(ipt,eta_bin,"" if isB else "non"))
        if jrt.GetEntries() < 100:
            continue
        mean, res = get_mean_res(jrt)
        gr.SetBinContent(ipt+1, res*100)
        gr.SetBinError(ipt+1, 0.00001)

c = r.TCanvas("c","c",800,600)
c.SetLogx()

hdummy = r.TH1F("h","",10,7,5000)
hdummy.GetYaxis().SetRangeUser(0,40)
hdummy.GetYaxis().SetTitle("Resolution (%)")
hdummy.GetYaxis().SetTitleSize(0.04)
hdummy.GetXaxis().SetTitle("p_{T}^{gen}")
hdummy.GetXaxis().SetTitleSize(0.04)
hdummy.GetXaxis().SetTitleOffset(1.10)
hdummy.Draw()

leg = r.TLegend(0.4, 0.6, 0.88, 0.88)
leg.SetBorderSize(0)
leg.SetNColumns(2)

grs_b = []
grs_nonb = []
colors = [r.kBlack, r.kRed, r.kBlue, r.kGreen+2]
for i,ieta in enumerate([0, 3, 5, 7]):
    grs_nonb.append(r.TH1D("h_nonb"+str(ieta), "", pt_bins.size-1, pt_bins))
    grs_b.append(r.TH1D("h_b"+str(ieta), "", pt_bins.size-1, pt_bins))
    GetResPlot(ieta, False, grs_nonb[-1])
    GetResPlot(ieta, True, grs_b[-1])

    grs_nonb[-1].SetMarkerStyle(20)
    grs_nonb[-1].SetMarkerSize(0.01)
    grs_nonb[-1].SetMarkerColor(colors[i])
    grs_nonb[-1].SetLineColor(colors[i])
    grs_nonb[-1].SetLineStyle(1)
    grs_nonb[-1].SetLineWidth(2)
    grs_nonb[-1].Draw("E SAME")
    grs_b[-1].SetMarkerStyle(20)
    grs_b[-1].SetMarkerSize(0.01)
    grs_b[-1].SetMarkerColorAlpha(colors[i],0.0)
    grs_b[-1].SetLineColor(colors[i])
    grs_b[-1].SetLineStyle(2)
    grs_b[-1].SetLineWidth(2)
    grs_b[-1].Draw("E SAME")

    leg.AddEntry(grs_nonb[-1], "{0} < |#eta| < {1}".format(eta_bins[ieta], eta_bins[ieta+1]), 'l')

leg.AddEntry(grs_nonb[0], "non-b jets", 'l')
leg.AddEntry(grs_b[0], "b jets", 'l')

leg.Draw()

text = r.TLatex()
text.SetNDC(1)
text.SetTextFont(62)
text.SetTextAlign(11)
text.SetTextSize(0.05)
text.DrawLatex(0.11,0.91,"CMS Simulation")
text.SetTextFont(42)
text.SetTextAlign(31)
text.SetTextSize(0.05)
text.DrawLatex(0.89,0.91,"13 TeV (2018)")

c.SaveAs("~/public_html/test.pdf")
#raw_input()





