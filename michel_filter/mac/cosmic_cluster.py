#!/usr/bin/python
from ROOT import TLatex
import ROOT       as rr
import root_numpy as rn
import numpy      as np

import sys

from looks   import *
from windows import *
from reader  import *

rr.gSystem.Load("libLArLite_DataFormat.so")
rr.gSystem.Load("libLArLite_LArUtil.so")
geo = rr.larutil.GeometryUtilities.GetME()  



if __name__ == '__main__':
    truefile = rr.TFile("/Users/vgenty/git/data/cosmics/larlite_mcinfo.root","READ")
    recofile = rr.TFile("/Users/vgenty/git/data/cosmics/larlite_reco2d.root","READ")
    
    evt = int(sys.argv[1])
    michels  = read_true_michels(truefile,evt,geo)

    clusters = make_clusters(extract_clusters(evt,recofile))
    
    EvtDisplay(clusters,michels,geo)
    
    
    
    
