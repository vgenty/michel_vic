from ROOT import *

f = TFile("the_dataa.root","READ")
f.tgraph_plane_2.Draw("AP")
f.tfit_plane_2.Draw("SAMES")
f.ONE.Draw("SAMES")
f.TWO.Draw("SAMES")
f.THREE.Draw("SAMES")

