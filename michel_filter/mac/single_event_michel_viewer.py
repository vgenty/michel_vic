#!/usr/bin/python
import ROOT as rr
from ROOT import TLatex
import root_numpy as rn
import numpy as np
import sys
from looks import *

from windows import *

from cluster import Cluster
from scipy.stats.stats import pearsonr

rr.gSystem.Load("libLArLite_DataFormat.so")
rr.gSystem.Load("libLArLite_LArUtil.so")
geo = rr.larutil.GeometryUtilities.GetME()  

#def Get2DPointProjection(x,y,z,plane):
#TODO???

looks_minos()
rr.gStyle.SetPalette(1)


def distance(p1,p2):
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2

def near(p1,p2,dx,dy):
    # if(abs(p2[0] - p1[0]) < 1.5 and 
    #    abs(p2[1] - p1[1]) < 5.65):   #(0.3/0.08)
    if(abs(p2[0] - p1[0]) < 1.0 and 
       abs(p2[1] - p1[1]) < 3.75):   #(0.3/0.08)
        return True
        
    #if(abs(p2[0] - p1[0]) < 0.61 and 
    #   abs(p2[1] - p1[1]) < 2.29):   #(0.3/0.08)
    #    return True
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
                    #charge b[s].Charge(0)/100000.0
                    #energy b[s].Start().E()
                    
    
    # Start().X/Y/Z() is useless since I can't get it's
    # projection onto the Y plane -.-
    return found_events


def extract_hits(evt,recofile):
    # Get the event you asked for, next set the branch as br
    recofile.hit_gaushit_tree.GetEntry(evt)
    br = recofile.hit_gaushit_tree.hit_gaushit_branch
    
    hits_xy     = []
    hits_xy_err = []
    charge      = []

    for h in xrange(len(br)): #lol br not iterable???
        if(br[h].View() == 2):
            x = br[h].WireID().Wire * 0.3       # from larsoft
            y = br[h].PeakTime()    * 0.0802814 # from larsoft

            
            #i made this shit up
            hits_xy_err.append(50.0/(br[h].Integral()))
            hits_xy.append([x,y])
            charge.append(br[h].Integral())
            #print "charge : " + str(br[h].Integral()) + "sigma charge : " + str(br[h].SigmaIntegral())

    #make into numpy array
    hits_xy = np.asarray(hits_xy)
    hits_xy_err = np.asarray(hits_xy_err)
    charge      = np.asarray(charge)
    return [hits_xy,hits_xy_err,charge]
    

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

def EvtDisplay(cobjects,hits_xy,p,kk,evt,charge) :
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
    th2.SetTitle(";Wire [cm]; Time [cm]")

    for i in xrange(len(hits_xy)) :
        th2.Fill(hits_xy[i][0],hits_xy[i][1],p[i])
    
    th2.Draw("COLZ")
    c2.Update()
    c2.Modified()


    c1 = rr.TCanvas("c1")
    c1.cd()
    
    tmg  = rr.TMultiGraph()
    #rest = rr.TGraph()
    tgs     = [rr.TGraph() for j in xrange(len(cobjects))]
    tgspts  = [rr.TGraph() for j in xrange(len(cobjects))]

    for i in xrange(len(cobjects)):
        for j in xrange(cobjects[i].size):
            tgs[i].SetPoint(j,hits_xy[cobjects[i].idxs[j]][0],hits_xy[cobjects[i].idxs[j]][1])
            
    for t in xrange(len(tgspts)):
        tgspts[t].SetPoint(0,hits_xy[cobjects[t].start][0],
                           hits_xy[cobjects[t].start][1])
        tgspts[t].SetPoint(1,hits_xy[cobjects[t].end][0],
                           hits_xy[cobjects[t].end][1])

    u = 2
    for tk in tgs:
        if(u == 0 or u == 10):
            u += 1
        tk.SetMarkerColor(u)
        tk.SetMarkerStyle(20)
        tmg.Add(tk)
        u += 1
    u = 2
    for tk in tgspts:
        if(u == 0 or u == 10):
            u += 1
        tk.SetMarkerColor(1)
        tk.SetMarkerStyle(29)
        tk.SetMarkerSize(2)
        tmg.Add(tk)
        u += 1

    
    the_key = 999

    for key in kk:
        if(kk[key]["event"] == int(evt)):
            the_key = key
            break

    tTruth    = rr.TGraph()
    tTruthMom = rr.TGraph()
    
    if(the_key != 999):
        loc    =  geo.Get2DPointProjectionVIC(kk[key]["startx"],
                                              kk[key]["starty"],
                                              kk[key]["startz"],
                                              2);

        tTruth.SetPoint(0,loc[0]*0.3,loc[1]*0.0802814)
        tTruth.SetMarkerStyle(34)
        tTruth.SetMarkerSize(2)
        tTruth.SetMarkerColor(30)
        tmg.Add(tTruth)
    
        
        
    #rn.fill_graph(rest,hits_xy[np.where(p < 0.8)])
    #rest.SetMarkerStyle(20)

    #tmg.Add(rest)
    tmg.Draw("AP")
    tmg.GetXaxis().SetTitle("Wire [cm]")
    tmg.GetYaxis().SetTitle("Time [cm]")
    c1.Update()
    c1.Modified()

    l  = 0
    ll = 0

    h = 0
    for c in cobjects:
        if (c.size > h) :
            h  = c.size
            ll = l
        l+=1
    
    c3 = rr.TCanvas()
    tdqdx = rr.TGraph()
    rn.fill_graph(tdqdx,np.array(cobjects[ll].dqdx))
    tdqdx.GetXaxis().SetTitle("s [cm]")
    tdqdx.GetYaxis().SetTitle("#frac{dq}{ds} [ADC/cm]")
    tdqdx.Draw("AL")
    
    
    c3.Update()
    c3.Modified()
    
    
    c4 = rr.TCanvas()
    tgc = rr.TGraph()
    rn.fill_graph(tgc,hits_xy[cobjects[ll].ordered_pts])
    tgc.GetXaxis().SetTitle("wire")
    tgc.GetYaxis().SetTitle("time")
    tgc.Draw("AL")
    
    c4.Update()
    c4.Modified()
    
    
    c5 = rr.TCanvas()
    tgcharge = rr.TGraph()
    ch = 0.0
    dd = 10.0
    print len(cobjects[ll].ordered_pts)
    #the_charge = [[ooo,charge[cobjects[ll].ordered_pts[ooo]]] for ooo in xrange(len(cobjects[ll].ordered_pts))]
    
    # the_charge = [charge[cobjects[ll].ordered_pts[ooo]] for ooo in xrange(len(cobjects[ll].ordered_pts))]
    # from itertools import count, tee, izip, islice
    # the_charge = map(np.median, izip(*(islice(it,i,None,10)
    #                                  for i, it in enumerate(tee(the_charge, 10)))))
    
    the_charge = []
    gggg = 0.0
    for w in xrange(len(cobjects[ll].ordered_pts) - 1):
        tgcharge.SetPoint(w,gggg,charge[cobjects[ll].ordered_pts[w]])
        gggg += cobjects[ll].ds[w]
        
        
    # for u in xrange(len(the_charge)):
    #     tgcharge.SetPoint(u,u,the_charge[u])
        
    # rn.fill_graph(tgcharge,the_charge)
    
    tgcharge.SetMarkerStyle(20)
    tgcharge.GetXaxis().SetTitle("s [cm]")
    tgcharge.GetYaxis().SetTitle("Charge [ADC]")
    tgcharge.Draw("AP")
    
    
    c5.Update()
    c5.Modified()

    c6 = rr.TCanvas()
    tgacharge = rr.TGraph()
    ch = 0.0
    dd = 10.0
    the_chargey = [charge[cobjects[ll].ordered_pts[ooo]] for ooo in xrange(len(cobjects[ll].ordered_pts))]

    a = 4
    b = 10

    # a = 5
    # b = 10
    
    from itertools import count, tee, izip, islice
    the_chargey = map(np.median, izip(*(islice(it,i,None,a)
                                        for i, it in enumerate(tee(the_chargey, b)))))
    # the_x = map(np.mean, izip(*(islice(it,i,None,a)
    #                         for i, it in enumerate(tee(hits_xy[cobjects[ll].ordered_pts][:,0]), b))))
    # the_y = map(np.mean, izip(*(islice(it,i,None,a)
    #                             for i, it in enumerate(tee(hits_xy[cobjects[ll].ordered_pts][:,1]), b))))
    
    for u in xrange(len(the_chargey)):
        tgacharge.SetPoint(u,u,the_chargey[u])
         
        
    tgacharge.SetMarkerStyle(20)
    tgacharge.GetXaxis().SetTitle("z [arbitrary]")
    tgacharge.GetYaxis().SetTitle("Median Charge [ADC]")
    tgacharge.Draw("AP")
    
    
    c6.Update()
    c6.Modified()

    
    c7     = rr.TCanvas()
    thchrage = rr.TH1D("ch",";;",100,0,0)
    rn.fill_hist(thchrage,charge[cobjects[ll].ordered_pts])
    # the_dirs = [k for k in cobjects[ll]]
    # gggg = 0.0
    # for w in xrange(len(cobjects[ll].ordered_pts)-1):
    #     tgcharge.SetPoint(w,gggg,XXXXXXXXXX)
    #     gggg += cobjects[ll].ds[w
    thchrage.GetXaxis().SetTitle("charge")
    thchrage.GetYaxis().SetTitle("count")
    
    thchrage.Draw()
    c7.Update()
    c7.Modified()


    
    c9 = rr.TCanvas(); tgmeans = rr.TGraph()
    
    meancharge = windowed_means(charge[cobjects[ll].ordered_pts],13,0,20)
    g = np.asarray([[i,meancharge[i]] for i in xrange(len(meancharge))])

    rn.fill_graph(tgmeans,g)
    tgmeans.Draw("ALP")
    c9.Update()
    c9.Modified()

    
    c8        = rr.TCanvas()
    c8.Update()
    c8.Modified()


    raw_input('')

    c1.Clear()
    c2.Clear()
    c3.Clear()
    c4.Clear()
    c5.Clear()
    c6.Clear()
    c7.Clear()
    c8.Clear()
    c9.Clear()

def remove_inside(cobjects,hits_xy):
    
    # check inside bounding rectangle

    yes = True
    before = 0
    after  = 0
    poop = 0
    while yes:
        final_objects = []
        poop = 0
        #print " before " + str(len(cobjects))
        #before = len(cobjects)
        for c in cobjects:
            for k in cobjects:
                if(c in cobjects and k in cobjects and c != k):
                    if(k.inside(c)):
                        final_objects.append(c + k)
                        cobjects.remove(c)
                        cobjects.remove(k)
                        poop += 1
        
        #print " after " + str(len(cobjects))
        #after = len(cobjects)
        cobjects += final_objects
        if(poop == 0): 
            yes = False
            #final_objects += cobjects
              
        
    #cobjects = final_objects
    return cobjects

    
def check_boundaries(cobjects,hits_xy):
    
    yes = True
    poop = 0
    while yes:
        poop = 0
        final_objects = []
        #for c, k, z in zip(cobjects,cobjects,cobjects):
        for c in cobjects:
            for k in cobjects:
                for z in cobjects:
                    if(c != k and c != z and k != z and
                       c in cobjects and k in cobjects and z in cobjects):
                        if(z.close_to(k) and z.close_to(c) and
                           c.size > 15 and k.size > 15): #this here could use some improvement
                            final_objects.append((c + k) + z)
                            cobjects.remove(c)
                            cobjects.remove(k)
                            cobjects.remove(z)
                            poop += 1

        cobjects += final_objects
        if(poop == 0): 
            yes = False

    return cobjects



def mostly_contained(cobjects,hits_xy):
    
    yes = True
    poop = 0
    while yes:
        final_objects = []
        poop = 0
        for c in cobjects:
            for k in cobjects:
                if(c in cobjects and k in cobjects 
                   and c != k and c.size > k.size):
                    if(k.mostly_contained(c,0.5)):
                        final_objects.append(c + k)
                        cobjects.remove(c)
                        cobjects.remove(k)
                        poop += 1
        
        #print " after " + str(len(cobjects))
        #after = len(cobjects)
        cobjects += final_objects
        if(poop == 0): 
            yes = False
            #final_objects += cobjects
    
    
    
    return cobjects

def stitch(cobjects,hits_xy):
    

    yes = True
    poop = 0
    while yes:
        poop = 0
        final_objects = []
        for c in cobjects:
            for k in cobjects:
                for z in cobjects:
                    if(c in cobjects and 
                       k in cobjects and 
                       z in cobjects and
                       c != k and c != z and k != z 
                       and z.size < c.size
                       and z.size < k.size):
                       #and z.pears < c.pears
                       #and z.pears < k.pears):
                        if(z.touching(c) and z.touching(k)):
                            final_objects.append((c + k) + z)
                            cobjects.remove(c)
                            cobjects.remove(k)
                            cobjects.remove(z)
                            poop += 1
        cobjects += final_objects
        if(poop == 0): 
            yes = False


    return cobjects

def dir(xy,err):
    invERR2 = 1.0/np.square(err)

    X  = (xy[:,0]*invERR2).sum()
    Y  = (xy[:,1]*invERR2).sum()
    XY = ((xy[:,0]*xy[:,1]*invERR2).sum())*(invERR2.sum())
    X2 = ((xy[:,0]*xy[:,0]*invERR2).sum())*(invERR2.sum())
    denom = X*X-X2
    return (-1.0)*(X*Y - XY)/(denom) #(-1.0 return (a,b)\cdot(x,y) for unit vector...

def calc_direction(c,hits_xy,hits_xy_err):
    nearbys = {k: 0.0 for k in c.idxs}
        
    for key in nearbys:
        close = [ l for l in c.idxs if near(hits_xy[key],hits_xy[l],1.0,1.0) and key != l] + [key]
        if(len(close) > 1):
            nearbys[key] = dir(hits_xy[close],
                               hits_xy_err[close])
        
    c.dirs = nearbys
    
def calculate_directions(cobjects,hits_xy,hits_xy_err):
    for c in cobjects:
        calc_direction(c,hits_xy,hits_xy_err)

    return cobjects
                            
def merge_final(cobjects):
    return cobjects
    
def order_points_dqdx(cobjects,hits_xy,charge):

    for c in cobjects : 
        c.order_points()
        c.fill_dqdx(charge)
    
    return cobjects

def add_kinks(cobjects,hits_xy) :

    yes = True
    poop = 0
    while yes:
        final_objects = []
        poop = 0
        for c in cobjects:
            for k in cobjects:
                if(c in cobjects and k in cobjects 
                   and c != k and c.size > k.size):
                    if(k.n_shared_boundaries(c) >= 2):
                        final_objects.append(c + k)
                        cobjects.remove(c)
                        cobjects.remove(k)
                        poop += 1
                        
        #print " after " + str(len(cobjects))
        #after = len(cobjects)
        cobjects += final_objects
        if(poop == 0): 
            yes = False
            #final_objects += cobjects
    
    
    
    return cobjects

if __name__ == '__main__':

    print "Parsing mc_info and reco files..."
    # 20 is a hard number here
    
    #For now let's just get the first folder.
    
    #truefiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_%d/larlite_mcinfo.root" % g,"READ")  for g in xrange(20)]
    #recofiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_%d/larlite_reco2d.root" % g,"READ")  for g in xrange(20)]

    
    
    truefiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_%d/larlite_mcinfo.root" % g,"READ")  for g in xrange(2)]
    recofiles = [rr.TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_%d/larlite_reco2d.root" % g,"READ")  for g in xrange(2)]

    pthresh = 0.9
    kk = find_michel_events(truefiles)
    print "Found michel events"

    hhh = extract_hits(int(sys.argv[1]),recofiles[0])
    hits_xy     = hhh[0]
    hits_xy_err = hhh[1]
    charge      = hhh[2]

        
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
                          sorted(xx, 
                                 key = lambda id: hits_xy[id][0])[0],
                          sorted(xx, 
                                 key = lambda id: hits_xy[id][0])[-1],
                          hits_xy)
                  for xx in clusters ]
    for l in left:
        cobjects.append(Cluster(l,
                                sorted(l, 
                                       key = lambda id: hits_xy[id][0])[0],
                                sorted(l, 
                                       key = lambda id: hits_xy[id][0])[-1],
                                hits_xy))
                        
                        
                        
    # look through clusters and find ones that are "inside" others....
    # currently this loop doesn't catch "overlapping ones"
    
    
    #cobjects = remove_inside(cobjects,hits_xy)
    cobjects = [c for c in cobjects if c.size > 1]
    
    #group clusters mostly contained in others
    cobjects = mostly_contained(cobjects,hits_xy)
    
    #pick three clusters, if one is touching two distinct big ones group
    cobjects = check_boundaries(cobjects,hits_xy)
    
    
    #group clusters mostly contained in others...
    cobjects = mostly_contained(cobjects,hits_xy)

    # #add "kinks" so if one cluster shares two boundary points with another then something is fishy, combine
    cobjects = add_kinks(cobjects,hits_xy)
    
    # #group clusters mostly contained in others one final time..................................
    cobjects = mostly_contained(cobjects,hits_xy)
    
    #cobjects = remove_inside(cobjects,hits_xy)
    #temporarily disable stitch
    #cobjects = stitch(cobjects,hits_xy)
    
    
    cobjects  = calculate_directions(cobjects,hits_xy,hits_xy_err)    
    
    cobjects  = order_points_dqdx(cobjects,hits_xy,charge)
    cobjects  = merge_final(cobjects)
    
    
    
    
    # at this point we are left with long individual clusters broken up by
    # somewhat linear ones, we need a connector algorithm
    
    
    

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
    
    ## probably as good as we are going to get, lets not try
    ## to add more clusters geometrically 
    
        
            
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
        
        
                
    EvtDisplay(cobjects,hits_xy,p,kk,sys.argv[1],charge)
           

