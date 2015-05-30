#!/usr/bin/python

import sys, os

from larlite import larlite
import ROOT

# Load my shared object
ROOT.gSystem.Load("libmichel_vic_triangles.so")

the_file = sys.argv[1]


# My shared object knows about larlite hopefully if it's in same namespace
proc = larlite.ana_processor()
proc.set_io_mode(larlite.storage_manager.kREAD)
proc.add_input_file(the_file)

tri_ana = larlite.Triangles()
proc.add_process(tri_ana)
proc.set_ana_output_file("the_dataa.root");

#process one event?
#proc.process_event(int(sys.argv[2]))

#call finalize explicitly
#proc.finalize()

proc.run()

sys.exit(0)
