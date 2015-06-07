#!/usr/bin/python
import ROOT as rr
import root_numpy as rn
import numpy as np
import sys

from cluster import Cluster
from scipy.stats.stats import pearsonr

rr.gSystem.Load("libLArLite_DataFormat.so")
rr.gSystem.Load("libLArLite_LArUtil.so")

#def Get2DPointProjection(x,y,z,plane):
#TODO???


def distance(p1,p2):
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2

def near(p1,p2,dx,dy):
    # if(abs(p2[0] - p1[0]) < 1.5 and 
    #    abs(p2[1] - p1[1]) < 5.65):   #(0.3/0.08)
    if(abs(p2[0] - p1[0]) < 1.0 and 
       abs(p2[1] - p1[1]) < 3.75):   #(0.3/0.08)
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
        if(np.isnan(pvals[i])) : pvals[i] = 0.0
       
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
    
    idx = Do_Grouping(idx)

    return [idx,pvals]

def Do_Grouping(idx):
    
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
        
    return idx

def EvtDisplay(cobjects,hits_xy,p) :
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
    #rest = rr.TGraph()
    tgs  = [rr.TGraph() for j in xrange(len(cobjects))]

    for i in xrange(len(cobjects)):
        for j in xrange(cobjects[i].size):
            tgs[i].SetPoint(j,hits_xy[cobjects[i].idxs[j]][0],hits_xy[cobjects[i].idxs[j]][1])
        
    u = 2
    for tk in tgs:
        if(u == 0 or u == 10):
            u += 1
        tk.SetMarkerColor(u)
        tk.SetMarkerStyle(20)
        tmg.Add(tk)
        u += 1

        
    #rn.fill_graph(rest,hits_xy[np.where(p < 0.8)])
    #rest.SetMarkerStyle(20)

    #tmg.Add(rest)
    tmg.Draw("AP")
    tmg.GetXaxis().SetTitle("Wire [cm]")
    tmg.GetYaxis().SetTitle("Time [cm]")
    c1.Update()
    c1.Modified()
    
    raw_input('')
    c1.Clear()
    c2.Clear()
    
if __name__ == '__main__':
    print "Parsing mc_info and reco files..."
    # 20 is a hard number here
    truefiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_%d/larlite_mcinfo.root" % g,"READ")  for g in xrange(20)]
    recofiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_%d/larlite_reco2d.root" % g,"READ")  for g in xrange(20)]

    pthresh = 0.9
    kk = find_michel_events(truefiles)
    print "Found michel events"

    hits_xy = extract_hits(int(sys.argv[1]),recofiles[0])
    cc = AhoCluster(hits_xy)
    clusters = cc[0]; p = cc[1]
    
    if(len(clusters) == 1): 
        print "Found only one cluster so who cares..."
        #continue;
        
    
    # convert clusters[i] to np.array and sort based on wire
    # clusters  = [np.array(sorted(c,key=lambda wire : hits_xy[wire][0]))
    #              for c in clusters]


    
    #some clusters are cut in half or more, try and add the indicies
    
    left_overs = [ i for i in xrange(len(hits_xy)) if not is_val_in_lists(i,clusters) ] 
    left_over_nears = {k : [] for k in left_overs}
    for key in left_over_nears:
        for i in left_overs:
            if(i != key):
                if(near(hits_xy[key],hits_xy[i],1,1)):
                    left_over_nears[key].append(i)
    
    left = [left_over_nears[key] + [key] for key in left_over_nears]
    left = Do_Grouping(left)
    
    # left holds the other clusters
    # cluster holds the main muon one

    # create cluster objects
    cobjects =  [ Cluster(xx,
                          p[xx].mean(),
                          p[xx].std(),
                          sorted(xx, 
                                 key = lambda id: hits_xy[id][0])[0],
                          sorted(xx, 
                                 key = lambda id: hits_xy[id][0])[-1],
                          p,
                          True) 
                  for xx in clusters ]
    for l in left:
        cobjects.append(Cluster(l,
                                p[l].mean(),
                                p[l].std(),
                                sorted(l, 
                                       key = lambda id: hits_xy[id][0])[0],
                                sorted(l, 
                                       key = lambda id: hits_xy[id][0])[-1],
                                p,
                                False))
        
    

    #look through clusters and find ones that are "inside" others....
    poop = 0
    yes = True
    final_objects = []
    deleted = []

    while yes:
        for c in cobjects:
            for k in cobjects:
                if(c != k ):
                    if(hits_xy[k.start][0] < hits_xy[c.start][0]
                       and 
                       hits_xy[k.end][0]   > hits_xy[c.end][0]):
                        final_objects.append(c + k)
                        cobjects.remove(c)
                        cobjects.remove(k)
                        poop += 1
        if(poop == 0): 
            yes = False
            final_objects += cobjects
        poop = 0
        
    # import copy
    # cobjects = copy.copy(final_objects)
    # final_objects = []
    # #look through clusters and find ones that are between two high pavg clusters
    # yes = True
    # while yes:
    #     for c in cobjects:
    #         for k in cobjects:
    #             for z in cobjects:
    #                 if(c in cobjects and 
    #                    k in cobjects and 
    #                    z in cobjects and
    #                    c != k and c != z and k != z 
    #                    and z.large == 0
    #                    and c.large == 1
    #                    and k.large == 1
    #                    and z.size < c.size
    #                    and z.size < k.size):
    #                     if(near(hits_xy[z.start],hits_xy[c.end],1,1)
    #                        and 
    #                        near(hits_xy[z.end],hits_xy[k.start],1,1)):
    #                         final_objects.append((c + k) + z)
    #                         print str(hits_xy[z.start]) + " near " + str(hits_xy[c.end]) + " and " + str(hits_xy[z.end]) + " near " + str(hits_xy[z.start])
    #                         print "z.large " + str(z.large) + " c.large " + str(c.large) + " k.large " + str(k.large)
    #                         cobjects.remove(c)
    #                         cobjects.remove(k)
    #                         cobjects.remove(z)
    #                         poop += 1

    #     if(poop == 0): 
    #         yes = False
    #         final_objects += cobjects
    #     poop = 0
    
        
    # cobjects = final_objects
    # xx.avg   = p[xx.idxs].mean()
    # xx.std   = p[xx.idxs].std()
    # xsorted  = sorted(xx.idxs, key = lambda id: hits_xy[id][0])
    # xx.start = xsorted[0]
    # xx.end   = xsorted[-1]


    #     min_keys = sorted(nears, key=lambda key: len(nears[key]))
    #     cstart   = min_keys.pop(0)
    #     cend     = 0.0

    #     #find the next key that is not "near"
    #     for i in min_keys:
    #         if (near(hits_xy[i],hits_xy[cstart],1,1) is False):
    #             cend = i
    #             break
       
    #     cobjects[o].start = cstart
    #     cobjects[o].end   = cend
    #     o += 1
 
            
        
            
    # For event display

    # while True:
    #     try:
    #         user_input_evt_no = input('Hit Enter to continue to next evt, or type in an event number to jump to that event:')
    #     except SyntaxError:
    #         user_input_evt_no = user_input_evt_no + 1
            

    #     print "Extracting hits..."
    #     hits_xy = extract_hits(int(user_input_evt_no),recofiles[0])
    #     print "Ahoclustering..."
    #     cc = AhoCluster(hits_xy)
    #     clusters = cc[0]; p = cc[1]
    #     print "Drawing..."
        
    
    
    
    EvtDisplay(cobjects,hits_xy,p)


