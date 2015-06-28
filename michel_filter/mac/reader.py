import numpy as np
from a_cluster import ACluster

from cluster_merger import *
from windows        import *
from looks          import *
from methods        import *

import ROOT as rr
import root_numpy as rn

looks_minos()
rr.gStyle.SetPalette(1)


def extract_clusters(evt,recofile) :
    
    recofile.cluster_fuzzycluster_tree.GetEntry(evt)
    recofile.hit_gaushit_tree.GetEntry(evt)
    recofile.ass_fuzzycluster_tree.GetEntry(evt)

    hrr = recofile.hit_gaushit_tree.hit_gaushit_branch
    crr = recofile.cluster_fuzzycluster_tree.cluster_fuzzycluster_branch
    ass = recofile.ass_fuzzycluster_tree.ass_fuzzycluster_branch


    cluster_to_hit_ass = ass.association(crr.id(),hrr.id())


    # for c in xrange(len(crr)):
    #     if(crr[c].View() == 2):
            
    
    clusterQ = 0.0
    count    = 0

    clusters = { c : [] for c in xrange(len(crr)) if crr[c].View() == 2 }

    for hit_indicies in cluster_to_hit_ass :
        if(crr[count].View() == 2):
            clusters[count] = np.array(hit_indicies)
        count += 1
        
    
    return [clusters,crr,hrr]
 
        # for hit_index in hit_indicies :
        #     count +=1 
        #     clusterQ += hrr[hit_index].Integral()
        
        

def make_clusters(c) :
    clusters, crr, hrr = c
    aclusters = []
    
    
    # Make the usual clusters
    for c in clusters:
        aclusters.append(ACluster(c,crr,hrr,clusters))


    # merge them
    
    aclusters = check_boundaries(aclusters)

    return aclusters


def read_true_michels(f,i,geo) :
    
    found_events = {}
    counter = 0;

    #nevents = f.mcshower_mcreco_tree.GetEntries()
    #for i in xrange(nevents):
    f.mcshower_mcreco_tree.GetEntry(i)
    b = f.mcshower_mcreco_tree.mcshower_mcreco_branch
    for s in xrange(len(b)):
        if(b[s].Process() == "muMinusCaptureAtRest" and
           b[s].Charge(2) > 1.0  #and
           #b[s].PdgCode() == 11  and
           #b[s].MotherPdgCode() == 13
       ):
            found_events[counter] = { "file"       : f,
                                      "event"      : i,
                                      "charge"     : b[s].Charge(2),
                                      "startx"     : b[s].Start().X(),
                                      "starty"     : b[s].Start().Y(),
                                      "startz"     : b[s].Start().Z(),
                                      "energy"     : b[s].Start().E(),
                                      "motherendX" : b[s].MotherEnd().X(),
                                      "motherendY" : b[s].MotherEnd().Y(),
                                      "motherendZ" : b[s].MotherEnd().Z()
            }
            
            counter += 1
            
            print "Found a michel"
            print "The XYZ is ("  + str(b[s].Start().X()) + " , " + str(b[s].Start().Y())  + " , " + str(b[s].Start().Z()) + ")"
            
            loc =  geo.Get2DPointProjectionVIC(b[s].Start().X(),
                                               b[s].Start().Y(),
                                               b[s].Start().Z(),
                                                   2);
            
            print "( " + str(loc[0]) + "," + str(loc[1]) + " )"
            
            print "event: " + str(i)
            
    return found_events
            
            
def EvtDisplay(clusters,michels,geo) :

    c1 = rr.TCanvas("c1")
    c1.Divide(2,2)
    c1.cd(1)

    
    tmg  = rr.TMultiGraph()
    tgs     = [rr.TGraph() for j in xrange(len(clusters))]

    for i in xrange(len(clusters)) :
        for j in xrange(clusters[i].hits_xy[:,0].size) :
            # print "i : " + str(i) + " j: " + str(j) + " ~ " + str(clusters[i].hits_xy[j])
            tgs[i].SetPoint(j,
                            clusters[i].hits_xy[j][0],
                            clusters[i].hits_xy[j][1])
            
    u = 2
    for tk in tgs:
        if(u == 0 or u == 10):
            u += 1
        tk.SetMarkerColor(u)
        tk.SetMarkerStyle(20)
        tmg.Add(tk)
        u += 1

    tTruth    = [rr.TGraph() for i in xrange(len(michels))]

    if(len(michels) > 0) :
        cntr = 0
        for key in michels:
            loc    =  geo.Get2DPointProjectionVIC(michels[key]["startx"],
                                                  michels[key]["starty"],
                                                  michels[key]["startz"],
                                                  2);

            tTruth[cntr].SetPoint(0,loc[0]*0.3,loc[1]*0.0802814)
            tTruth[cntr].SetMarkerStyle(34)
            tTruth[cntr].SetMarkerSize(2)
            tTruth[cntr].SetMarkerColor(30)
            tmg.Add(tTruth[cntr])
            cntr += 1
    
        
        
    tmg.Draw("AP")
    setaxis(tmg,"Wire [cm]","Time [cm]")

    

    c1.cd(2)
    tgmeans = rr.TGraph()
    tc = clusters[-1]
    
    charges = [h.Integral() for h in tc.hits[tc.ordered_pts] ]
               
    meancharge = windowed_means(charges,25,0,25)
    meancharge.pop(0); meancharge.pop(-1)
    meancharge.pop(0); meancharge.pop(-1)

    g = np.asarray([[tc.s[i],meancharge[i]] for i in xrange(len(meancharge))])

    rn.fill_graph(tgmeans,g)
    tgmeans.Draw("ALP")
    setaxis(tgmeans,"s [cm]","Q [ADC]")


    c1.cd(4)
    tgderive  = rr.TGraph()

    s = 2
    baka = []
    for i in xrange(s,len(g)-s+1) :
        baka.append(smooth_derive(g[i-s : i+s][:,1],g[i-s : i+s][:,0],2*s+1))
    
    baka = [[i,baka[i]] for i in xrange(len(baka))]
    rn.fill_graph(tgderive,baka)
    setaxis(tgderive,"s [cm]","#frac{dQ}{ds}")

    tgderive.Draw("ALP")


    c1.cd(3)
    tgorder  = rr.TGraph()
    aho = tc.hits_xy[tc.ordered_pts]
    rn.fill_graph(tgorder,aho)
    setaxis(tgorder,"Wire [cm]","Time [cm]")
    
    tgorder.Draw("AL")

    c1.Update()
    c1.Modified()

    
    
    raw_input('')

    c1.Clear()
