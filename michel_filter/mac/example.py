import sys
from ROOT import gSystem
gSystem.Load("libmichel_vic_michel_filter")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing michel_filter..."

sys.exit(0)

