import numpy as np
from scipy.stats.stats import pearsonr
from pyhull.convex_hull import ConvexHull;

class Cluster:
    def __init__(self,idx,start,end,hits_xy):
        self.idxs   = idx
        self.size   = len(idx)
        self.start  = start
        self.end    = end
        self.pears  = abs(pearsonr(hits_xy[idx,:][:,0],hits_xy[idx,:][:,1])[0])
        self.hitsxy = hits_xy
        self.boundary = []
        if(self.size > 2):
            hull = ConvexHull(self.hitsxy[self.idxs])
            for v in hull.vertices:
                if(v[0] not in self.boundary): self.boundary.append(v[0])

        
    def __add__(self, other):
        start = float(0.0)
        end   = float(0.0)
    
        if   (self.start < other.start) and (self.end < other.end):
            start = self.start
            end   = other.end
        elif (self.start < other.start) and (self.end > other.end):
            start = self.start
            end   = self.end
        elif (self.start > other.start) and (self.end > other.end):
            start = other.start
            end   = self.end
        elif (self.start > other.start) and (self.end < other.end):
            start = other.start
            end   = other.end

        return Cluster(self.idxs + other.idxs,
                       start,
                       end,
                       self.hitsxy)

    def inside(self,other):
        if((self.hitsxy[self.start][0] >= self.hitsxy[other.start][0]
            and 
            self.hitsxy[self.end][0]   <= self.hitsxy[other.end][0])
           and
           ((self.hitsxy[self.start][1] <= self.hitsxy[other.start][1]
             and 
             self.hitsxy[self.end][1]   >= self.hitsxy[other.end][1])
            or
            (self.hitsxy[self.start][1] >= self.hitsxy[other.start][1]
             and 
             self.hitsxy[self.end][1]   <= self.hitsxy[other.end][1]))):
            return True
        else:
            return False
            
    def near(self,p1,p2,dx,dy):
    # if(abs(p2[0] - p1[0]) < 1.5 and 
    #    abs(p2[1] - p1[1]) < 5.65):   #(0.3/0.08)
        if(abs(p2[0] - p1[0]) < 1.0 and 
           abs(p2[1] - p1[1]) < 3.75):   #(0.3/0.08)
            return True
        return False

    def touching(self,other):
        if(self.near(self.hitsxy[self.start],self.hitsxy[other.end],1,1)
           or
           self.near(self.hitsxy[self.end],self.hitsxy[other.start],1,1)):
            return True
        else:
            return False

    def close_to(self,other):
        for b1 in self.boundary:
            for b2 in other.boundary:
                #print "b1: " + str(b1) + " b2: " + str(b2)
                if(self.near(self.hitsxy[b1],self.hitsxy[b2],1,1)):
                    #print "is near...."
                    return True

        return False
