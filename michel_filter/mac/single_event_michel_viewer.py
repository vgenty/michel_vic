#!/usr/bin/python
import ROOT as rr
import root_numpy as rn
import numpy as np
import sys
from scipy.stats.stats import pearsonr

rr.gSystem.Load("libLArLite_DataFormat.so")
rr.gSystem.Load("libLArLite_LArUtil.so")

#def Get2DPointProjection(x,y,z,plane):
#TODO???

def near(p1,p2,dx,dy):
    if(abs(p2[0] - p1[0]) < 1.5 and 
       abs(p2[1] - p1[1]) < 3.0):
        return True
    return False
    
def is_val_in_lists(w,lis):
    return w in [item for sublist in lis for item in sublist]

def is_val_in_list(w,lis):
    return w in lis;

def list_overlap(l1,l2):
    return any(xx in l1 for xx in l2)
    
def find_michel_events(TF):
    found_events = {}
    counter = 0;
    for f in TF :
        nevents = f.mcshower_mcreco_tree.GetEntries()
        for i in xrange(nevents):
            f.mcshower_mcreco_tree.GetEntry(i)
            b = f.mcshower_mcreco_tree.mcshower_mcreco_branch
            for s in xrange(len(b)):
                if(b[s].Process() == "muMinusCaptureAtRest" and
                   b[s].Charge(2) > 50.0 and
                   b[s].PdgCode() == 11  and
                   b[s].MotherPdgCode() == 13):
                    found_events[counter] = { "file"   : f,
                                              "event"  : i,
                                              "charge" : b[s].Charge(2),
                                              "startx" : b[s].Start().X(),
                                              "starty" : b[s].Start().Y(),
                                              "startz" : b[s].Start().Z(),
                                              "energy" : b[s].Start().E()}
                    counter += 1



                    #charge b[s].Charge(0)/100000.0
                    #energy b[s].Start().E()
                    
    
    # Start().X/Y/Z() is useless since I can't get it's
    # projection onto the Y plane -.-
    return found_events


def extract_hits(evt,recofile):
    # Get the event you asked for, next set the branch as br
    recofile.hit_gaushit_tree.GetEntry(evt)
    br = recofile.hit_gaushit_tree.hit_gaushit_branch
    
    hits_xy = []
    
    for h in xrange(len(br)): #lol br not iterable???
        if(br[h].View() == 2):
            x = br[h].WireID().Wire * 0.3       # from larsoft
            y = br[h].PeakTime()    * 0.0802814 # from larsoft
        
            hits_xy.append([x,y])

    #make into numpy array
    hits_xy = np.asarray(hits_xy)
    return hits_xy

def AhoCluster(hits_xy):
    
    # idx will hold indicies of nearby hits
    # for example idx[5] = [6,7,45]  means the 5'th indexed
    # hit (hits_xy[5]) has neighbors hits_xy[6,7,45]
    nhits = len(hits_xy)
    idx = [ [j for j in xrange(nhits) if near(hits_xy[i],hits_xy[j],1.5,1.5) and i != j] for i in xrange(nhits) ]

    
    # we will use pearsonsr as a way to find regions of interest
    pvals = np.zeros(nhits)
    for i in xrange(nhits): 
        pvals[i] = abs(pearsonr(hits_xy[idx[i],:][:,0],hits_xy[idx[i],:][:,1])[0])     
       
    #now idx is technically a list of clusters so lets start
    #ignoring that idx index matches hits_xy

    idx = [idx[i] + [i] for i in xrange(nhits)]
    
    pthresh = 0.9
    
    #clusters have bad hits in them (pvals < pthresh), remove them
    u = 0
    bad = True
    while(bad):
        for idxs in idx:
            for i in idxs:
                if(pvals[i] < pthresh):
                    idxs.remove(i)
                    u += 1
        if(u == 0): bad = False
        u = 0
    
    #there will now be some empty lists in idx, remove them
    idx = [x for x in idx if x != []]
    
    #you will also find that some lists are duplicated, this 
    #following lines helps remove some not all i don't understand why
    idx = sorted(idx)
    idx = [idx[i] for i in range(len(idx)) if i == 0 or sorted(idx[i]) != sorted(idx[i-1])]

    #here is the meat, find overlaping clusters, confusing...
    isoverlap = True
    qq = 0
    while(isoverlap):
        for c in idx:
            for z in idx:
                if( set(c) != set(z) ):
                    if(any(xx in c for xx in z)):
                        c.extend([o for o in z if o not in c])
                        idx.remove(z)
                        qq += 1
        if(qq == 0): isoverlap = False
        qq = 0
        
    return [idx,pvals]
    
def EvtDisplay(clusters,hits_xy,p) :
        #Draw something useful...

    c2 = rr.TCanvas("c2")
    c2.cd()
    
    xmin = np.amin(hits_xy[:,0]) - 10.0
    xmax = np.amax(hits_xy[:,0]) + 10.0
    ymin = np.amin(hits_xy[:,1]) - 10.0
    ymax = np.amax(hits_xy[:,1]) + 10.0
    
    nbinsx = int(np.ceil((xmax-xmin)/0.3))
    nbinsy = int(np.ceil((ymax-ymin)/0.08))

    th2 = rr.TH2D("baka",";;",nbinsx,xmin,xmax,nbinsy,ymin,ymax)
        
    for i in xrange(len(hits_xy)) :
        th2.Fill(hits_xy[i][0],hits_xy[i][1],p[i])
    
    th2.Draw("COLZ")
    c2.Update()
    c2.Modified()


    c1 = rr.TCanvas("c1")
    c1.cd()
    
    tmg  = rr.TMultiGraph()
    rest = rr.TGraph()
    tgs  = [rr.TGraph() for j in xrange(len(clusters))]
    for i in xrange(len(clusters)):
        for j in xrange(len(clusters[i])):
            tgs[i].SetPoint(j,hits_xy[clusters[i][j]][0],hits_xy[clusters[i][j]][1])
        
    u = 2
    for tk in tgs:
        tk.SetMarkerColor(u)
        tk.SetMarkerStyle(20)
        tmg.Add(tk)
        u += 1

        
    rn.fill_graph(rest,hits_xy[np.where(p < 0.9)])
    rest.SetMarkerStyle(20)

    tmg.Add(rest)
    tmg.Draw("AP")
    tmg.GetXaxis().SetTitle("Wire [cm]")
    tmg.GetYaxis().SetTitle("Time [cm]")
    c1.Update()
    c1.Modified()
    
    print "Just showed Event %d. Hit enter to go next event..." % user_input_evt_no
    raw_input('')
    c1.Clear()
    c2.Clear()
    
if __name__ == '__main__':
    print "Parsing mc_info and reco files..."
    # 20 is a hard number here
    truefiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_%d/larlite_mcinfo.root" % g,"READ")  for g in xrange(20)]
    recofiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_%d/larlite_reco2d.root" % g,"READ")  for g in xrange(20)]

    k = find_michel_events(truefiles)
    print "Found michel events"
    
    while True:
        try:
            user_input_evt_no = input('Hit Enter to continue to next evt, or type in an event number to jump to that event:')
        except SyntaxError:
            user_input_evt_no = user_input_evt_no + 1
            

        print "Extracting hits..."
        hits_xy = extract_hits(int(user_input_evt_no),recofiles[0])
        print "Ahoclustering..."
        cc = AhoCluster(hits_xy)
        clusters = cc[0]; p = cc[1]
        print "Drawing..."
        
        EvtDisplay(clusters,hits_xy,p)


