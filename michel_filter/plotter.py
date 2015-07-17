import numpy as np
import sys

from looks          import *
from methods        import *

import ROOT as rr
import root_numpy as rn

looks_minos()
rr.gStyle.SetOptStat("emr")
rr.gStyle.SetPalette(1)
#rr.gStyle.SetOptStat(1111)

event = int(sys.argv[1])
evt_num = event
f    = rr.TFile.Open("output.root","READ")
tree = f.Get("out_tree")

rec  = rn.tree2rec(tree)

c9 = rr.TCanvas()   #1
c10 = rr.TCanvas()  #2
c11 = rr.TCanvas()  #3
c12 = rr.TCanvas()  #4

tgX = rr.TGraph()
tgY = rr.TGraph()

thX = rr.TH1D("thX",";;",250,-100,100)
thY = rr.TH1D("thY",";;",250,-100,100)

Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])

rn.fill_graph(tgX,Xs)
rn.fill_graph(tgY,Ys)
rn.fill_hist(thX,rec['reco_X'] - rec['true_X'])
rn.fill_hist(thY,rec['reco_Y'] - rec['true_Y'])

c9.cd()
c9.cd()
tgX.SetMarkerStyle(20)
tgX.Draw("AP")
setaxis(tgX,"Reco Z [cm]", "True Z [cm]")
c9.Update()
c9.Modified()


c10.cd()
tgY.SetMarkerStyle(20)
tgY.Draw("AP")
setaxis(tgY,"Reco X [cm]", "True X [cm]")
c10.Update()
c10.Modified()

c11.cd()
thX.Draw()
setaxis(thX,"Reco - True Z [cm]","Count/%.2f [cm]" % thX.GetXaxis().GetBinWidth(0))
c11.Update()
c11.Modified()

c12.cd()
thY.Draw()
setaxis(thY,"Reco - True X [cm]","Count/%.2f [cm]" % thY.GetXaxis().GetBinWidth(0))
c12.Update()
c12.Modified()

c1 = rr.TCanvas()
# c1.Divide(4,2)

c1.cd()

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

c1.Update()
c1.Modified()

c3 = rr.TCanvas()
c3.cd()

tmgs = rr.TMultiGraph()
ordered = rr.TGraph()
tstart  = rr.TGraph()
tend    = rr.TGraph()

order = np.array([[rec['ahits_X'][evt_num][k],rec['ahits_Y'][evt_num][k]] 
                  for k in rec['ordered_pts'][evt_num]])

ordered.SetMarkerSize(0)
rn.fill_graph(ordered,order)
tstart.SetPoint(0,rec['startX'][evt_num],rec['startY'][evt_num])
tend.SetPoint  (0,rec['endX'][evt_num]  ,rec['endY'][evt_num])

tstart.SetMarkerStyle(29)
tstart.SetMarkerSize(2)
tstart.SetMarkerColor(6)

tend.SetMarkerStyle(20)
tend.SetMarkerSize(2)
tend.SetMarkerColor(7)

tmgs.Add(ordered)
tmgs.Add(tstart)
tmgs.Add(tend)

tmgs.Draw("ALP")
setaxis(tmgs,"Wire [cm]","Time [cm]")


c3.Update()
c3.Modified()

c2 = rr.TCanvas()
c2.cd()

mcharges = rr.TMultiGraph()
meancharge = rr.TGraph()
meanchargepeaks = rr.TGraph()

thecharge =  np.array([[rec['s'][evt_num][k],rec['mean_charges'][evt_num][k]] 
                       for k in xrange(rec['mean_charges'][evt_num].size)])
rn.fill_graph(meancharge,thecharge)

we = 0;
for peak in rec['_the_tmean_max_peak'][evt_num] :
    meanchargepeaks.SetPoint(we,thecharge[peak][0],thecharge[peak][1])
    we += 1
    
mcharges.Add(meancharge)
mcharges.Add(meanchargepeaks)

mcharges.Draw("ALP")
setaxis(mcharges,"s [cm]","Truncated Q [ADC]")

c2.Update()
c2.Modified()


c4 = rr.TCanvas()
c4.cd()

tdqds  = rr.TGraph()
dqds   = np.array([[rec['s'][evt_num][k],rec['dqds'][evt_num][k]] 
                   for k in xrange(rec['dqds'][evt_num].size)])
rn.fill_graph(tdqds,dqds)
tdqds.Draw("ALP")
setaxis(tdqds,"s [cm]","Truncated dQ/ds [ADC/cm]")

c4.Update()
c4.Modified()


c5 = rr.TCanvas()
c5.cd()

truecharge = rr.TGraph()

cc = [rec['charges'][evt_num][k] for k in rec['ordered_pts'][evt_num]]

dd = np.array([[rec['s'][evt_num][i],cc[i]] for i in xrange(len(cc))])
rn.fill_graph(truecharge,dd)
truecharge.Draw("ALP")
setaxis(truecharge,"s [cm]","Q [ADC]")

c5.Update()
c5.Modified()

c6 = rr.TCanvas()
c6.cd()

allhits = rr.TGraph()
shohits = rr.TGraph()

ppp= rr.TMultiGraph()

zzz = np.array([[rec['_ALL_hits_p2_X'][evt_num][k],rec['_ALL_hits_p2_Y'][evt_num][k]] 
                for k in xrange(len(rec['_ALL_hits_p2_X'][evt_num]))])
qqq = np.array([[rec['_large_frac_shower_hits_X'][evt_num][k],rec['_large_frac_shower_hits_Y'][evt_num][k]] 
                for k in xrange(len(rec['_large_frac_shower_hits_X'][evt_num]))])


rn.fill_graph(allhits,zzz)
rn.fill_graph(shohits,qqq)

shohits.SetMarkerColor(2)

allhits.SetMarkerStyle(20)
shohits.SetMarkerStyle(20)

ppp.Add(allhits)
ppp.Add(shohits)

ppp.Draw("AP")
setaxis(ppp,"Wire [cm]","Time [cm]")

c6.Update()
c6.Modified()

c7 = rr.TCanvas()
c7.cd()


xmin = np.amin(order[:,0]) - 10.0
xmax = np.amax(order[:,0]) + 10.0
ymin = np.amin(order[:,1]) - 10.0
ymax = np.amax(order[:,1]) + 10.0

nbinsx = int(np.ceil((xmax-xmin)/0.3))
nbinsy = int(np.ceil((ymax-ymin)/0.08))

th2 = rr.TH2D("baka",";;",nbinsx,xmin,xmax,nbinsy,ymin,ymax)
th2.SetTitle(";Wire [cm]; Time [cm]")
for i in xrange(len(order)) :
    th2.Fill(order[i][0],order[i][1],rec['_chi2_copy'][evt_num][i])

th2.Draw("COLZ")

c7.Update()
c7.Modified()

c8 = rr.TCanvas()
c8.cd()

ccc      = rr.TMultiGraph()
chipeaks = rr.TGraph()
chiS     = rr.TGraph()

chhh = np.array([[rec['s'][evt_num][k],rec['_chi2_copy'][evt_num][k]] 
                 for k in xrange(rec['s'][evt_num].size)])
we = 0;
for peak in rec['_the_chi_max_peak'][evt_num] :
    chipeaks.SetPoint(we,chhh[peak][0],chhh[peak][1])
    we += 1
    
    
rn.fill_graph(chiS,chhh)
chipeaks.SetMarkerColor(4)
#chipeaks.SetMarkerSize(2)
chipeaks.SetMarkerStyle(20)
    
ccc.Add(chiS)
ccc.Add(chipeaks)
ccc.Draw("AP")

setaxis(ccc,"s [cm]","#chi^{2}/NDF")

# Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
# Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])


c8.Update()
c8.Modified()

c13 = rr.TCanvas()
c13.cd()

#simch_plane

simch = rr.TH1D("simch", "", 100, 0, 14000)
rrr = np.array([rec['_simch_plane_true_shower_E'][k] for k in xrange(rec['_simch_plane_true_shower_E'].size)])

rn.fill_hist(simch, rrr)
simch.Draw()

setaxis(simch ,"True Shower N electrons","Count")

c13.Update()
c13.Modified()

c14 = rr.TCanvas()
c14.cd()

lifetime = rr.TH1D("lifetime", "", 100, 0.9, 1.8)

sss = np.array([rec['_lifetime_correction'][k] for k in xrange(rec['_lifetime_correction'].size)])

rn.fill_hist(lifetime, sss)
lifetime.Draw()
setaxis(lifetime,"exp[t/ #tau]","")

c14.Update()
c14.Modified()

c15 = rr.TCanvas()
c15.cd()

true_michel = rr.TH1D("true_michel", "", 100, 0, 60)

lll = np.array([rec['_true_michel_Det'][k] for k in xrange(rec['_true_michel_Det'].size)])
rn.fill_hist(true_michel, lll)


true_michel.Draw()
setaxis(true_michel, "True Michel Energy [MeV]","Count/%.2f MeV" % true_michel.GetXaxis().GetBinWidth(0))

c15.Update()
c15.Modified()


c16 = rr.TCanvas()
c16.cd()

mcQfrac = rr.TH1D("mcQfrac", "", 100, 0, 1.1)

nnn = np.array([rec['_mcQ_frac'][k] for k in xrange(rec['_mcQ_frac'].size)])
rn.fill_hist(mcQfrac, nnn)
mcQfrac.Draw()
setaxis(mcQfrac,"Efficiency","Count/%.2f" % mcQfrac.GetXaxis().GetBinWidth(0))

c16.Update()
c16.Modified()

c17 = rr.TCanvas()
c17.cd()

E = rr.TGraph()

mmm = np.array([[rec['_simch_plane_true_shower_E'][i], rec['_michel_E'][i]] for i in xrange(rec['_michel_E'].size)])
   
rn.fill_graph(E, mmm)
setaxis(E,"True Shower Num e^{-}","Reco Michel Q [ADC]")

E.Draw("AP")

c17.Update()
c17.Modified()

    
raw_input('')

#sys.stdin.readline()
