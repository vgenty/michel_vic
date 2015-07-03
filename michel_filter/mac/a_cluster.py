import numpy as np
from pyhull.convex_hull import ConvexHull;
import matplotlib.path as mplPath;
import ROOT

class ACluster:    
    def __init__(self,idx,clusters,hits,ass):
        self.cluster = [clusters[idx]]
        self.hits    = np.array([ hits.at(ROOT.Long(h)) for h in ass[idx] ])

        _xy = lambda h : [h.WireID().Wire * 0.3, h.PeakTime() * 0.0802814]
        self.hits_xy      = np.array([ _xy(h) for h in self.hits ])
        
        self.left2right  = sorted([i for i in xrange(self.hits.size)], 
                                  key = lambda id: self.hits_xy[id][0])
        self.start       = self.left2right[ 0]
        self.end         = self.left2right[-1]
        
        self.ds          = []
        self.s           = []
        self.ordered_pts = []

        self.order_points()

    
    
    def __add__(self, other):
        self.cluster += other.cluster
        self.hits     = np.concatenate((self.hits, other.hits), 
                                       axis=0)
        self.hits_xy  = np.concatenate((self.hits_xy,other.hits_xy),
                                       axis=0)
        
        # python allows only one constructor so I have to redo this
        self.left2right  = sorted([i for i in xrange(self.hits.size)], 
                                  key = lambda id: self.hits_xy[id][0])

        self.start       = self.left2right[ 0]
        self.end         = self.left2right[-1]
        
        self.ds          = []
        self.s           = []
        self.ordered_pts = []

        self.order_points()
        
        
        return self

    def distance(self,p1,p2):
        return np.sqrt(  (self.hits_xy[p1][0] - self.hits_xy[p2][0])**2 
                       + (self.hits_xy[p1][1] - self.hits_xy[p2][1])**2 )

    def order_points(self):
        self.ordered_pts.append(self.start)
        self.s.append(0)
        idxx = self.left2right[:]

        aho = True
        closest = 9999.0
        zz = 0.0
        cnt = 0
        stot = 0
        idxholder = 0
        
        idxx.remove(self.start)
        
        j = 0
        while aho:
            for i in idxx:
                zz = self.distance(self.ordered_pts[cnt],i)
                if(zz < closest and zz < 0.3*6):
                    idxholder = i
                    closest   = zz
                    j = 1
                    

            
            if(j == 1):        
                idxx.remove(idxholder)            
                self.ordered_pts.append(idxholder)
                self.ds.append(closest)
                stot += closest
                self.s.append(stot)
                
                idxholder = 0
                closest = 9999.0
                zz   = 0.0
                cnt += 1
            
            if(len(idxx) == 0 or j == 0) :
                aho = False
                
            j = 0

    def near(self,p1,p2,dx,dy):
        if(abs(p2[0] - p1[0]) < 1.0 and 
           abs(p2[1] - p1[1]) < 3.75):   #(0.3/0.08)
            return True
        return False

    def touching(self,other) :
        if(self.near(self.hits_xy[self.start],other.hits_xy[other.end],1,1)
           or
           self.near(self.hits_xy[self.end],other.hits_xy[other.start],1,1)):
            return True
        else:
            return False

        
        
        
        
        
        
