
# coding: utf-8

# In[71]:

#!/usr/bin/python -i
from scipy.stats.stats import pearsonr
from numpy import array, zeros, empty, append,asarray
import numpy as np
from ROOT import *
from root_numpy import *
gSystem.Load("libLArLite_DataFormat.so")

ff = [TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_%d/larlite_reco2d.root" % g,"READ")  for g in xrange(20)]


def near(p1,p2):
    if(abs(p2[0] - p1[0]) < 2.0 and
       abs(p2[1] - p1[1]) < 2.0):
        return True;
        
    return False;

evt = 273
ff[0].hit_gaushit_tree.GetEntry(evt)
b = ff[0].hit_gaushit_tree.hit_gaushit_branch

#k = zeros([len(b),2])
k = []

for h in xrange(len(b)):
    if(b[h].View() == 0) :
        
        x = b[h].WireID().Wire * 0.3;
        y = b[h].PeakTime()    * 0.0802814;
        
        k.append([x,y])
    
k = asarray(k)
idx = []

#event hits are now in k
#find the indicies in k that are near

for i in xrange(len(k)):
    local_nears = []
    for j in xrange(len(k)):
        if(i != j):
            if(near(k[i],k[j])):
                local_nears.append(j)
    idx.append(local_nears)


# we will use pearsonsr as a way to find regions of interest

p = zeros(len(k))
for i in xrange(len(k)) :
    p[i] = abs(pearsonr(k[idx[i],:][:,0],k[idx[i],:][:,1])[0])




# c = 0
# for i in xrange(len(idx)):
#     for j in xrange(len(idx)):
#         if(i != j):
#             if (i in idx[j]):
#                 c += 1
#     print "Found " + str(i) + ", " + str(c) + " times"
#     c = 0

cluster = []

 
c1 = TCanvas()
tg = TGraph()

fill_graph(tg,k)
tg.SetMarkerStyle(20)
tg.Draw("AP")


# for f in ff :
#     nevents = f.hit_gaushit_tree.GetEntries()
#     nevents = 1
#     for i in xrange(nevents):
#         f.hit_gaushit_tree.GetEntry(i)
#         b = f.hit_gaushit_tree.hit_gaushit_branch
        
        

# th1d = TH1D("xx",";;",100,0,100)

# for f in ff :
#     nevents = f.mcshower_mcreco_tree.GetEntries()
#     for i in xrange(nevents):
#         f.mcshower_mcreco_tree.GetEntry(i)
#         b = f.mcshower_mcreco_tree.mcshower_mcreco_branch
#         for s in xrange(len(b)):
#             if(b[s].Process() == "muMinusCaptureAtRest" and 
#                b[s].Charge(0) > 50.0 and
#                b[s].PdgCode() == 11  and
#                b[s].MotherPdgCode() == 13):
#                 print "On event: %d in shower %d" % (i,s) 
#                 th1d.Fill(b[s].Charge(0)/100000.0)
                
# th1d.Draw()

    
c1.Update()
c1.Modified()

xmin = np.amin(k[:,0]) - 10.0
xmax = np.amax(k[:,0]) + 10.0
ymin = np.amin(k[:,1]) - 10.0
ymax = np.amax(k[:,1]) + 10.0

nbinsx = int(np.ceil((xmax-xmin)/0.3))
nbinsy = int(np.ceil((ymax-ymin)/0.08))

th2 = TH2D("xx",";;",nbinsx,xmin,xmax,nbinsy,ymin,ymax)

for i in xrange(len(k)) :
    th2.Fill(k[i][0],k[i][1],p[i])
    
c2 = TCanvas()    
th2.Draw("COLZ")
c2.Update()
c2.Modified()


# In[65]:

def is_inn(w,lis):
    return w in [item for sublist in lis for item in sublist]

def is_in(w,lis):
    #return any(w in sublist for sublist in [lis])
    # if(w in [item for sublist in lis for item in sublist]):
    #     return True
    # return False
    return w in lis;


# In[76]:

cluster = []; clusters = []

clusters.append([])
# In[77]:
bad = []
for i in np.argsort(k[:,0]):
    if(p[i] < 0.9):
        if(cluster): clusters.append(cluster)
        cluster = []
        bad.append(i)
        continue
        
    for j in idx[i]:
        if(p[j] > 0.9 and not is_in(j,cluster)
           and not is_inn(j,clusters)):
            cluster.append(j)
        

    if(i == len(idx) - 1 and cluster):
        clusters.append(cluster)
   

clusters.pop(0)
# # In[78]:

clusters


# In[75]:

np.where(p < 0.9)


# In[102]:

len(clusters);



# In[96]:

import ROOT as root


# In[97]:

tmg = root.TMultiGraph()


# In[98]:

c3 = TCanvas()
c3.cd()

# In[99]:

tgs = [TGraph() for j in xrange(len(clusters))]


# In[104]:

for i in xrange(len(clusters)):
    for j in xrange(len(clusters[i])):
        tgs[i].SetPoint(j,k[clusters[i][j]][0],k[clusters[i][j]][1])


# In[105]:

tgs


# In[109]:

c=1
for tk in tgs:
    tk.SetMarkerColor(c+1)
    tk.SetMarkerStyle(20)
    c+=1
    tmg.Add(tk)
    
rest = TGraph()
fill_graph(rest,k[np.where(p < 0.9)])
rest.SetMarkerStyle(20)

tmg.Add(rest)
# 

# In[117]:

tmg.Draw("AP")
c3.Update()
c3.Modified()

# In[ ]:




# In[ ]:



