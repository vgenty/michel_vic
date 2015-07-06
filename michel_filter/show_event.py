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

f    = rr.TFile.Open("output.root","READ")
tree = f.Get("out_tree")
rec  = rn.tree2rec(tree)

Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])

c1 = rr.TCanvas()
c1 = rr.TCanvas("c1")
c1.Divide(3,2)

c1.cd(1)

tmg    = rr.TMultiGraph()
tcc    = rr.TGraph()
ttrue  = rr.TGraph()
treco  = rr.TGraph()
    
ahits = np.array([[rec['ahits_X'][evt_num][i],rec['ahits_Y'][evt_num][i]] 
                  for i in xrange(rec['ahits_X'][evt_num].size)])

rn.fill_graph(tcc,ahits)

ttrue.SetPoint(0,Xs[evt_num][1],Ys[evt_num][1])
treco.SetPoint(0,Xs[evt_num][0],Ys[evt_num][0])

ttrue.SetMarkerStyle(34)
ttrue.SetMarkerSize(2)
ttrue.SetMarkerColor(30)

treco.SetMarkerStyle(23)
treco.SetMarkerSize(2)
treco.SetMarkerColor(40)

tmg.Add(tcc)
tmg.Add(ttrue)
tmg.Add(treco)

tmg.Draw("AP")
setaxis(tmg,"Wire [cm]","Time [cm]")

c1.cd(3)

ordered = rr.TGraph()
order = np.array([[rec['ahits_X'][evt_num][k],rec['ahits_Y'][evt_num][k]] 
                  for k in rec['ordered_pts'][evt_num]])

rn.fill_graph(ordered,order)
ordered.Draw("AL")
setaxis(ordered,"Wire [cm]","Time [cm]")

c1.cd(2)
meancharge = rr.TGraph()
thecharge =  np.array([[rec['s'][evt_num][k],rec['mean_charges'][evt_num][k]] 
                       for k in xrange(rec['mean_charges'][evt_num].size)])
rn.fill_graph(meancharge,thecharge)
meancharge.Draw("ALP")

setaxis(meancharge,"s [cm]","Q [ADC]")

c1.cd(4)
tdqds  = rr.TGraph()
dqds   =  np.array([[rec['s'][evt_num][k],rec['dqds'][evt_num][k]] 
                    for k in xrange(rec['dqds'][evt_num].size)])
rn.fill_graph(tdqds,dqds)
tdqds.Draw("ALP")
setaxis(tdqds,"s [cm]","dQ/ds [ADC/cm]")

c1.cd(5)
truecharge = rr.TGraph()

cc = [rec['charges'][evt_num][k] for k in rec['ordered_pts'][evt_num]]
#ss = [rec['s'][evt_num][k] for k in xrange(['ordered_pts'][evt_num].size)]

dd = np.array([[rec['s'][evt_num][i],cc[i]] for i in xrange(len(cc))])
rn.fill_graph(truecharge,dd)
truecharge.Draw("ALP")
setaxis(truecharge,"s [cm]","Charge [ADC]")


c1.Update()
c1.Modified()
# Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
# Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])


