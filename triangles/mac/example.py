import sys
from ROOT import gSystem
gSystem.Load("libmichel_vic_triangles")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing triangles..."

sys.exit(0)

