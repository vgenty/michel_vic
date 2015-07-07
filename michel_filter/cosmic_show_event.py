import numpy as np
import sys

from looks          import *
from methods        import *

import ROOT as rr
import root_numpy as rn

looks_minos()
rr.gStyle.SetPalette(1)
rr.gStyle.SetOptStat(1111)


evt_num = int(sys.argv[1])

f    = rr.TFile.Open("cosmic_output.root","READ")
tree = f.Get("out_tree")
tree.GetEntry(evt_num)

tmg = rr.TMultiGraph()

tgAHITS = [rr.TGraph() for i in xrange(tree.ahits_X.size())]
tgTRUE  = [rr.TGraph() for i in xrange(tree.true_X.size()) ]

for i in xrange(len(tgAHITS)):
    for j in xrange(tree.ahits_X[i].size()):
        tgAHITS[i].SetPoint(j,tree.ahits_X[i][j],tree.ahits_Y[i][j])
        tgAHITS[i].SetMarkerColor(i)

for i in xrange(len(tgTRUE)):
    tgTRUE[i].SetPoint(0,tree.true_X[i],tree.true_Y[i])
    tgTRUE[i].SetMarkerStyle(34)
    tgTRUE[i].SetMarkerSize(3)
    #tgTRUE[i].SetMarkerColor(30)

c1 = rr.TCanvas()
c1.Divide(2)
c1.cd(1)
for a in tgAHITS:
    tmg.Add(a)

for b in tgTRUE:
    tmg.Add(b)


tmg.Draw("AP")


c1.cd(2)
tgORDER  = [rr.TGraph() for i in xrange(tree.ahits_X.size())]
tmgORDER = rr.TMultiGraph()


c = 0
for i in xrange(len(tgAHITS)):
    for j in tree.ordered_pts[i]:
        tgORDER[i].SetPoint(c,tree.ahits_X[i][j],tree.ahits_Y[i][j])
        c+=1
    c = 0
    tmgORDER.Add(tgORDER[i])

tmgORDER.Draw("AL")

c1.Update()
c1.Modified()


