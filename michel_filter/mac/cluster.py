import numpy as np


class Cluster:
    def __init__(self,idx,mean,std,start,end,p,large):
        self.idxs   = idx
        self.size   = len(idx)
        self.start  = start
        self.end    = end
        self.pavg   = mean
        self.pstd   = std
        self.p      = p
        self.new    = large
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
                       self.p[self.idxs + other.idxs].mean(),
                       self.p[self.idxs + other.idxs].std(),
                       start,
                       end,
                       self.p,
                       bool(self.new) + bool(other.new))
            
