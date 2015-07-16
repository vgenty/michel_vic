import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk
import ROOT
ROOT.gSystem.Load("libmichel_vic_michel_filter.so")

# Create ana_processor instance
my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kREAD)
my_proc.set_ana_output_file("output.root");
my_proc.enable_filter(True)

the_filter = fmwk.MichelFilter()
the_ana    = fmwk.Michel2DAna("fuzzycluster");

#~~tune-able parameters, they have defaults...
the_ana.set_min_merge_cluster_size(25)
the_ana.set_min_proto_cluster_size(4)
the_ana.set_n_window_size(15)      
the_ana.set_window_cutoff(0.20)
the_ana.set_truncated_shave(3)
the_ana.set_min_rad(10)
the_ana.set_threshold(0); #min value to open window
the_ana.set_rise(5); #number sigma above baseline
the_ana.set_fall(10); #number sigma below baseline

#you absolutely must set these two until Reco2D becomes singelton :(
nwires = float(4.0)
the_ana.set_near_X(nwires*0.3)      # cm ~~ 1.0
the_ana.set_near_Y(nwires*0.3/0.08) # cm ~~ 3.75
the_ana.set_d_cutoff(20.0*0.3)      # cm ~~ 20*0.3

my_proc.add_process(the_filter)
my_proc.add_process(the_ana)

#~~Let's run it.
my_proc.run();

#~~Or process single event
# my_proc.process_event(302)
# the_ana.finalize()

sys.exit(0)
