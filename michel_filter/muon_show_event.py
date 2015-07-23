import numpy as np
import sys

from looks          import *
from methods        import *

import ROOT as rr
import root_numpy as rn

def graph(event):

    looks_minos()
    rr.gStyle.SetPalette(1)
    #rr.gStyle.SetOptStat(1111)

    event = int(event)
    evt_num = event
    #f    = rr.TFile.Open("output_muons_only.root","READ")
    f    = rr.TFile.Open(sys.argv[1],"READ")
    tree = f.Get("out_tree")

    rec  = rn.tree2rec(tree)

    Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
    Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])

    c1 = rr.TCanvas()
    c1 = rr.TCanvas("c1")
    c1.Divide(4,2)

    c1.cd(1)

    tmg    = rr.TMultiGraph()
    tcc    = rr.TGraph()
    #ttrue  = rr.TGraph()
    treco  = rr.TGraph()

    ahits = np.array([[rec['ahits_X'][evt_num][i],rec['ahits_Y'][evt_num][i]] 
                      for i in xrange(rec['ahits_X'][evt_num].size)])

    rn.fill_graph(tcc,ahits)

    #ttrue.SetPoint(0,Xs[evt_num][1],Ys[evt_num][1])
    treco.SetPoint(0,Xs[evt_num][0],Ys[evt_num][0])

    # ttrue.SetMarkerStyle(34)
    # ttrue.SetMarkerSize(2)
    # ttrue.SetMarkerColor(30)

    treco.SetMarkerStyle(23)
    treco.SetMarkerSize(2)
    treco.SetMarkerColor(40)

    tmg.Add(tcc)
    #tmg.Add(ttrue)
    tmg.Add(treco)

    tmg.Draw("AP")
    setaxis(tmg,"Wire [cm]","Time [cm]")

    c1.cd(3)

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
    
    c1.cd(2)
    
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
    if(meanchargepeaks.GetN() > 0):
        mcharges.Add(meanchargepeaks)

    meanchargepeaks.SetMarkerColor(5)
    meanchargepeaks.SetMarkerStyle(20)


    mcharges.Draw("AP")
    setaxis(mcharges,"s [cm]","Q [ADC]")
    
    c1.cd(4)
    mdqds     = rr.TMultiGraph()
    tdqds     = rr.TGraph()
    dqdspeaks = rr.TGraph()

    dqds   = np.array([[rec['s'][evt_num][k],rec['dqds'][evt_num][k]] 
                       for k in xrange(rec['dqds'][evt_num].size)])
    
    rn.fill_graph(tdqds,dqds)

    we = 0;
    for peak in rec['_the_tdqds_min_peak'][evt_num] :
        dqdspeaks.SetPoint(we,dqds[peak][0],dqds[peak][1])
        we += 1
    
    
    mdqds.Add(tdqds)
    mdqds.Add(dqdspeaks)
    
    dqdspeaks.SetMarkerColor(6)
    dqdspeaks.SetMarkerStyle(20)

    mdqds.Draw("AP")
    setaxis(mdqds,"s [cm]","dQ/ds [ADC/cm]")

    c1.cd(5)
    truecharge = rr.TGraph()

    cc = [rec['charges'][evt_num][k] for k in rec['ordered_pts'][evt_num]]
    
    dd = np.array([[rec['s'][evt_num][i],cc[i]] for i in xrange(len(cc))])
    rn.fill_graph(truecharge,dd)
    truecharge.Draw("ALP")
    setaxis(truecharge,"s [cm]","Charge [ADC]")

    
    c1.cd(6)
    
    allhits = rr.TGraph()
    shohits = rr.TGraph()
    
    ppp= rr.TMultiGraph()

    zzz = np.array([[rec['_ALL_hits_p2_X'][evt_num][k],rec['_ALL_hits_p2_Y'][evt_num][k]] 
                    for k in xrange(len(rec['_ALL_hits_p2_X'][evt_num]))])
    # qqq = np.array([[rec['_large_frac_shower_hits_X'][evt_num][k],rec['_large_frac_shower_hits_Y'][evt_num][k]] 
    #                 for k in xrange(len(rec['_large_frac_shower_hits_X'][evt_num]))])
    
    
    rn.fill_graph(allhits,zzz)
    #rn.fill_graph(shohits,qqq)
    
    #shohits.SetMarkerColor(2)

    allhits.SetMarkerStyle(20)
    #shohits.SetMarkerStyle(20)
    
    ppp.Add(allhits)
    #ppp.Add(shohits)

    ppp.Draw("AP")
    setaxis(ppp,"Z","X")


    
    c1.cd(7)
    
    
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
    
    
    c1.cd(8)
    
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
    
    setaxis(ccc,"s [cm]","Chi2/NDF")

    # Xs = np.array([[rec['reco_X'][i],rec['true_X'][i]] for i in xrange(rec['reco_X'].size)])
    # Ys = np.array([[rec['reco_Y'][i],rec['true_Y'][i]] for i in xrange(rec['reco_Y'].size)])

    
    c1.Update()
    c1.Modified()

    
    raw_input('')

    #sys.stdin.readline()


event = 0
current = event
while True:
    print "Which event? (enter a number to choose or return to move to next event)"
    print "give me -1 to quit, + to move forward and - to move back"
    current = event
    event = raw_input()
    
    if event  == '+':
        print "vic says poo"
        current += 1
        event = current

        print "you gave me event: " + str(event)
       # current = event;
       # print current
        graph(event)

    if event == '-':
        current -= 1
        event = current

        print "you gave me event: " + str(event)
       # current = event;
           # print current
        graph(event)
        
    else:
        event = int(event)
        if event == -1 :
            break;


        if event >= 0 :
            print "you gave me event: " + str(event)
            graph(event)
        




