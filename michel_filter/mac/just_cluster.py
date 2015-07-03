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
    truefile = rr.TFile("~/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_7/larlite_mcinfo.root","READ")
    recofile = rr.TFile("~/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791010_7/larlite_reco2d.root","READ")
    
    evt = int(sys.argv[1])
    michels  = read_true_michels(truefile,evt,geo)

    clusters = make_clusters(extract_clusters(evt,recofile))
    
    EvtDisplay(clusters,michels,geo)
    
    
    
    
