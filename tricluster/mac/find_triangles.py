#!/usr/bin/python

import sys, os

from larlite import larlite
import ROOT

# Load my shared object
ROOT.gSystem.Load("libmichel_vic_tricluster.so")

the_file = sys.argv[1]

# My shared object knows about larlite hopefully if it's in same namespace
proc = larlite.ana_processor()
proc.set_io_mode(larlite.storage_manager.kREAD)
proc.add_input_file(the_file)

tri_ana = larlite.TriClusterAna()
proc.add_process(tri_ana)
#proc.set_ana_output_file("the_dataa.root");

#process one event?
proc.process_event(int(sys.argv[2]))



artist = tri_ana.get_artist()
p = artist.get_chosen();

hits     = artist.get_hits(p)
hits_err = artist.get_hits_err(p)

triangle = artist.get_triangle(p)

c1 = ROOT.TCanvas()
tg = ROOT.TGraphErrors()
for i in xrange(hits.size()) :
    tg.SetPoint     (i,hits[i].first,     hits[i].second)
    tg.SetPointError(i,hits_err[i].second,hits_err[i].second)
    
up   = ROOT.TF1("up","%f*x + %f" % (triangle.line_one  ().first,triangle.line_one  ().second),0,1100)
down = ROOT.TF1("up","%f*x + %f" % (triangle.line_two  ().first,triangle.line_two  ().second),0,1100)
hyp  = ROOT.TF1("up","%f*x + %f" % (triangle.line_three().first,triangle.line_three().second),0,1100)
    
tg.Draw("AP")
up.Draw("SAMES")
down.Draw("SAMES")
hyp.Draw("SAMES")



#call finalize explicitly
#proc.finalize()

#proc.run()

#sys.exit(0)
