import ROOT as r
import numpy as np

r.gROOT.ProcessLine(".L ~/scripts/rootalias.C")

r.gStyle.SetOptStat(0)
r.gStyle.SetNumberContours(255)

fin = r.TFile("combined.root")
h = fin.Get("auxiliary/h_lowsmear")

c = r.TCanvas()
c.SetTickx(1)
c.SetTicky(1)
c.SetLogz(1)
h.SetTitle("eta/phi for jets w/ p_{T}(gen) > 300, p_{T}(reco/gen) < 0.5")
h.GetXaxis().SetTitle("eta")
h.GetYaxis().SetTitle("phi")

h.Draw("COLZ")
# h.ProjectionY().Draw()

# line = r.TLine()
# for eta in np.arange(-3,3.01,0.1):
#     line.DrawLine(eta,-3.10,eta,3.10)
# line.SetLineColor(r.kRed)
# for phi in np.arange(-3.2,3.21,0.1):
#     line.DrawLine(-2.90,phi,2.90,phi)

boundaries = [
    (-2.4,-1.6,-2.15,-1.3),
    (-2.15,-2.6,-1.95,-2.35),
    (-1.65,0.5,-1.5,0.7),
    (-1.35,-2.75,-1.25,-2.6),
    (-1.3,0.5,-1.2,0.6),
    (-1.15,0.5,-1.05,0.6),
    (-0.4,0.1,-0.2,0.3),
    (-1.0,2.9,-0.85,3.1),
    (0.85,2.8,1.0,2.9),
    (0.25,-0.8,0.35,-0.7),
    (0.6,-0.8,0.7,-0.7),
    (0.95,-0.8,1.05,-0.7),
    (1.35,-2.65,1.45,-2.5),
    (1.65,0.7,1.85,0.95),
    (1.65,-0.7,1.8,-0.5),
    (1.5,-1.6,1.65,-1.4),
    (1.7,-2.2,1.9,-1.95),
    (-1.55,2.6,-1.35,2.7)
]

box = r.TBox()
box.SetLineColor(r.kRed)
box.SetFillStyle(0)
for b in boundaries:
    box.DrawBox(*b)

tot = 0.0
tota = 0.0
for b in boundaries:
    tot += r.GetIntegral2D(h, b[0], b[2], b[1], b[3])
    tota += (b[2]-b[0])*(b[3]-b[1])
print "{0:.2f}%".format(tot / h.Integral(0,-1,0,-1) * 100)
print "{0:.2f}%".format(tota / (2*2.4*2*np.pi) * 100)

print "\n Paste into looper:"
print "const uint NDEADCELLS = {0};".format(len(boundaries))
print "float dead_cell_boundaries[NDEADCELLS][4] = {"
for i,b in enumerate(boundaries):
    print "    {{{0:5.2f},{1:5.2f},{2:5.2f},{3:5.2f}}}{4}".format(b[0], b[1], b[2], b[3], "," if i<len(boundaries)-1 else "")
print "};"
print ""

# c.SaveAs("~/public_html/mt2/RebalanceAndSmear/low_smears_eta_phi_afterFilter_afterManualRemoval.png")
raw_input()

