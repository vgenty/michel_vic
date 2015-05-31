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

while True:
    
    try:
        user_input_evt_no = input('Hit Enter to continue to next evt, or type in an event number to jump to that event:')
    except SyntaxError:
        user_input_evt_no = user_input_evt_no + 1

    proc.process_event(int(user_input_evt_no))

    artist = tri_ana.get_artist()
    p = artist.get_chosen();

    hits     = artist.get_hits(p)
    hits_err = artist.get_hits_err(p)

    triangle = artist.get_triangle(p)

    c1    = ROOT.TCanvas()
    tmg   = ROOT.TMultiGraph()
    tgIN  = ROOT.TGraphErrors()
    tgOUT = ROOT.TGraphErrors()

    for i in xrange(triangle.get_fHit().size()) :
        tgIN.SetPoint     (i,hits    [triangle.get_fHit()[i]].first, hits    [triangle.get_fHit()[i]].second)
        tgIN.SetPointError(i,hits_err[triangle.get_fHit()[i]].second,hits_err[triangle.get_fHit()[i]].second)
    for i in xrange(hits.size()) :
        tgOUT.SetPoint     (i,hits    [i].first, hits    [i].second)
        tgOUT.SetPointError(i,hits_err[i].second,hits_err[i].second)

        
    xmin_one   = 0.0
    xmax_one   = 0.0

    xmin_two   = 0.0
    xmax_two   = 0.0

    xmin_thr = 0.0
    xmax_thr = 0.0

    x1 = (triangle.line_one().second -  triangle.line_two().second)/(triangle.line_one().first - triangle.line_two().first)
    x2 = (triangle.line_one().second -  triangle.line_three().second)/(triangle.line_one().first - triangle.line_three().first)
    x3 = (triangle.line_three().second - triangle.line_two().second)/(triangle.line_three().first - triangle.line_two().first)

    
    if(x1 < 0.0): x1 = -1.0*x1
    if(x2 < 0.0): x2 = -1.0*x2
    if(x3 < 0.0): x3 = -1.0*x3
    
    

    if(triangle.is_right()) :
        xmin_one   = x1
        xmin_two   = x1
        xmax_one   = x2
        xmax_two   = x3

    else:
        xmax_one = x1
        xmax_one = x1
        xmin_one = x2
        xmin_two = x3

    if(triangle.line_three().first > 0):
        xmax_thr = xmax_one
        xmin_thr = xmax_two
    else:
        xmax_thr = xmax_two
        xmin_thr = xmax_one

    

    up   = ROOT.TF1("up","%f*x + %f" % (triangle.line_one  ().first,triangle.line_one  ().second),xmin_one,xmax_one)
    down = ROOT.TF1("up","%f*x + %f" % (triangle.line_two  ().first,triangle.line_two  ().second),xmin_two,xmax_two)
    hyp  = ROOT.TF1("up","%f*x + %f" % (triangle.line_three().first,triangle.line_three().second),xmin_thr,xmax_thr)
    
    tgIN.SetMarkerColor(4)
    tgIN.SetMarkerStyle(20)
    tgOUT.SetMarkerColor(3)
    tgOUT.SetMarkerStyle(20)

    tmg.Add(tgOUT)
    tmg.Add(tgIN)

    tmg.Draw("AP")
    up.Draw("SAMES")
    down.Draw("SAMES")
    hyp.Draw("SAMES")

    c1.Update()
    c1.Modified()
    
    print "Hit enter to choose the next event..."
    sys.stdin.readline()
    
    del c1
    del tgIN
    del tgOUT
    del tmg
    del up
    del down
    del hyp



#call finalize explicitly
#proc.finalize()

#proc.run()

#sys.exit(0)
