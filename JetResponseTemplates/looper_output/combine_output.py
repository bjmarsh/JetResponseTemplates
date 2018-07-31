import ROOT as r
import glob

n_pt_bins = 23
n_eta_bins = 12

files = glob.glob("qcd_pt*.root")

hs = []
for ipt in range(n_pt_bins):
    hs.append([])
    for ieta in range(n_eta_bins):
        hs[-1].append([])
        hs[-1][-1].append(r.TH1D("JRT_pt{0}_eta{1}_bjet".format(ipt,ieta),"p_{T}^{reco}/p_{T}^{gen}",150,0,3))
        hs[-1][-1].append(r.TH1D("JRT_pt{0}_eta{1}_nonbjet".format(ipt,ieta),"p_{T}^{reco}/p_{T}^{gen}",150,0,3))
        hs[-1][-1].append(r.TH1D("JRT_pt{0}_eta{1}_alljet".format(ipt,ieta),"p_{T}^{reco}/p_{T}^{gen}",150,0,3))

        for f in files:
            extra_sf = 1.0

            # if f.startswith("qcd_pt800to1000"):
            #     extra_sf = 2.07386
            # if f.startswith("qcd_pt600to800"):
            #     extra_sf = 1.087379

            fid = r.TFile(f)
            h = fid.Get("pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_bjet".format(ipt,ieta))
            h.Scale(h.GetEntries() * extra_sf)
            hs[-1][-1][0].Add(h)

            h = fid.Get("pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_nonbjet".format(ipt,ieta))
            h.Scale(h.GetEntries() * extra_sf)
            hs[-1][-1][1].Add(h)

            h = fid.Get("pt{0}/pt{0}_eta{1}/JRT_pt{0}_eta{1}_alljet".format(ipt,ieta))
            h.Scale(h.GetEntries() * extra_sf)
            hs[-1][-1][2].Add(h)

            fid.Close()

fout = r.TFile("combined.root","RECREATE")
for ipt in range(n_pt_bins):
    ptdir = fout.mkdir("pt{0}".format(ipt))
    for ieta in range(n_eta_bins):
        etadir = ptdir.mkdir("pt{0}_eta{1}".format(ipt,ieta))
        etadir.cd()
        if hs[ipt][ieta][0].Integral(0,-1) > 0:
            hs[ipt][ieta][0].Scale(1.0 / hs[ipt][ieta][0].Integral(0,-1) / hs[ipt][ieta][0].GetBinWidth(1))
        if hs[ipt][ieta][1].Integral(0,-1) > 0:
            hs[ipt][ieta][1].Scale(1.0 / hs[ipt][ieta][1].Integral(0,-1) / hs[ipt][ieta][1].GetBinWidth(1))
        if hs[ipt][ieta][2].Integral(0,-1) > 0:
            hs[ipt][ieta][2].Scale(1.0 / hs[ipt][ieta][2].Integral(0,-1) / hs[ipt][ieta][2].GetBinWidth(1))
        hs[ipt][ieta][0].Write()
        hs[ipt][ieta][1].Write()
        hs[ipt][ieta][2].Write()
        
fout.Close()
