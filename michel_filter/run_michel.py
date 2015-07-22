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
#the_filter = fmwk.RemoveMichel()

#~~tune-able parameters, they have defaults...
the_ana.set_min_merge_cluster_size(25)
the_ana.set_min_proto_cluster_size(4)
the_ana.set_n_window_size(15)      
the_ana.set_window_cutoff(0.25)
the_ana.set_truncated_shave(3)
the_ana.set_min_rad(10)

the_ana.set_chi2_threshold(0); #min value to open window
the_ana.set_chi2_rise(5);      #number sigma above baseline
the_ana.set_chi2_fall(6);     #number sigma below baseline

the_ana.set_tmean_threshold(0);
the_ana.set_tmean_rise(10);       
the_ana.set_tmean_fall(10); 

the_ana.set_tdqds_threshold(0);
the_ana.set_tdqds_rise(10);       
the_ana.set_tdqds_fall(10);      

#you absolutely must set these two until Reco2D becomes singelton :(
nwires = float(4.0)
the_ana.set_near_X(nwires*0.3)      # cm ~~ 1.0
the_ana.set_near_Y(nwires*0.3/0.08) # cm ~~ 3.75
the_ana.set_d_cutoff(20.0*0.3)      # cm ~~ 20*0.3

my_proc.add_process(the_filter)
my_proc.add_process(the_ana)

my_proc.set_verbosity(fmwk.msg.kDEBUG)

#~~Let's run it.
my_proc.run()

#~~Or process single event
# my_proc.process_event(302)
# the_ana.finalize()

sys.exit(0)
