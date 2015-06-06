#!/usr/bin/python -i

from ROOT import *
gSystem.Load("libLArLite_DataFormat.so")

ff = [TFile("/Users/vgenty/git/data/prod_muminus_0.1-2.0GeV_isotropic_uboone/1791009_%d/larlite_mcinfo.root" % g,"READ")  for g in xrange(20)]

c1 = TCanvas()
thcharge = TH1D("xx",";;",50,0,100)
thenergy = TH1D("xx2",";;",50,0,100)


for f in ff :
    nevents = f.mcshower_mcreco_tree.GetEntries()
    for i in xrange(nevents):
        f.mcshower_mcreco_tree.GetEntry(i)
        b = f.mcshower_mcreco_tree.mcshower_mcreco_branch
        for s in xrange(len(b)):
            if(b[s].Process() == "muMinusCaptureAtRest" and 
               b[s].Charge(0) > 50.0 and
               b[s].PdgCode() == 11  and
               b[s].MotherPdgCode() == 13):
                print "On event: %d in shower %d" % (i,s) 
                thcharge.Fill(b[s].Charge(0)/100000.0)
                thenergy.Fill( b[s].Start().E())

thcharge.SetFillColor(4)
thenergy.SetFillColor(3)
thcharge.SetFillStyle(3004)
thenergy.SetFillStyle(3005)

thcharge.Draw()
thenergy.Draw("SAMES")

thcharge.GetXaxis().SetTitle("Charge/100000? and Energy")
c1.Update()
c1.Modified()



