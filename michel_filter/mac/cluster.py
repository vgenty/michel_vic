import numpy as np
from scipy.stats.stats import pearsonr
from pyhull.convex_hull import ConvexHull;
import matplotlib.path as mplPath;

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
            hull  = ConvexHull(self.hitsxy[self.idxs])
            vert  = hull.vertices
            first = vert.pop(0)
            self.boundary.append(first[0])
            self.boundary.append(first[1])
            while (len(vert) > 0):
                for l in xrange(len(vert)):
                    if(self.boundary[-1] == vert[l][0]):
                        self.boundary.append(vert[l][1])
                        vert.remove(vert[l])
                        break
        
            self.boundary.pop()
        self.boundary = [self.idxs[x] for x in self.boundary]
        self.path = mplPath.Path(np.array(self.hitsxy[self.boundary]))
        self.dirs = {}

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

    def mostly_contained(self,other,frac):
        tot_points = 0
        for i in self.idxs:
            if(other.path.contains_point(self.hitsxy[i])):
                tot_points += 1
                    
        if (float(frac)*self.size < tot_points):
            return True
        
        return False
