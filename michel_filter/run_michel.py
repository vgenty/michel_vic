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

the_ana = fmwk.Michel2DAna("fuzzycluster");

my_proc.add_process(the_ana)
#my_proc.add_process(fmwk.josh())

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
#my_proc.run();
my_proc.process_event(1)
# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
