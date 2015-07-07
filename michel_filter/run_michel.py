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
#the_ana    = fmwk.Michel2DAna("fuzzycluster");
the_ana    = fmwk.Michel2DCosmics("fuzzycluster");

my_proc.add_process(the_filter)
my_proc.add_process(the_ana)

# Let's run it.
my_proc.run();
# my_proc.process_event(302)
# the_ana.finalize()
# done!

sys.exit(0)
