import numpy as np
import sys

from looks          import *
from methods        import *

import ROOT as rr
import root_numpy as rn

looks_minos()
rr.gStyle.SetPalette(1)
rr.gStyle.SetOptStat(1111)

f    = rr.TFile.Open("output.root","READ")
tree = f.Get("out_tree")
rec  = rn.tree2rec(tree)


c1 = rr.TCanvas()
c2 = rr.TCanvas()
c3 = rr.TCanvas()
c4 = rr.TCanvas()

tgX = rr.TGraph()
tgY = rr.TGraph()

thX = rr.TH1D("thX",";;",100,0,0)
thY = rr.TH1D("thY",";;",100,0,0)


Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])

rn.fill_graph(tgX,Xs)
rn.fill_graph(tgY,Ys)
rn.fill_hist(thX,rec['reco_X'] - rec['true_X'])
rn.fill_hist(thY,rec['reco_Y'] - rec['true_Y'])


c1.cd()
tgX.SetMarkerStyle(20)
tgX.Draw("AP")
setaxis(tgX,"Reco X [cm]", "True X [cm]")
c1.Update()
c1.Modified()

c2.cd()
tgY.SetMarkerStyle(20)
tgY.Draw("AP")
setaxis(tgY,"Reco Y [cm]", "True Y [cm]")

c2.Update()
c2.Modified()

c3.cd()

thX.Draw()
setaxis(thX,"Reco - True X [cm]","Count")
c3.Update()
c3.Modified()

c4.cd()
thY.Draw()
setaxis(thY,"Reco - True Y [cm]","Count")
c4.Update()
c4.Modified()
