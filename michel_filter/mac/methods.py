from ROOT import TGraph, TMultiGraph, TLegend, TPaveText

def tgs(tgraphs,title) :

    tmg = TMultiGraph()
    tle = TLegend(0.6,0.4,0.9,0.8)
    tit = TPaveText(0.2988506,0.934322,0.6997126,1,"nbNDC")

    # Make title
    tit.SetFillColor(0)
    tit.SetFillStyle(0)
    tit.SetLineColor(0)
    tit.AddText(title)
    
    for i in xrange(len(tgraphs)):
        tgraphs[i].SetLineWidth(2)
        tgraphs[i].SetLineColor(i+1)

        tmg.Add(tgraphs[i])
        tle.AddEntry(tgraphs[i],tgraphs[i].GetName(),"l")

    
    return [tmg,tle,tit]

def setaxis(obj,X,Y):
    
    try:
        # Method exists, and was used.
        obj.GetXaxis().SetTitle(X)
        obj.GetYaxis().SetTitle(Y)
        obj.GetXaxis().CenterTitle()
        obj.GetYaxis().CenterTitle()

    except AttributeError:
        # Method does not exist.
        print "Object " + str(obj) + " can not access TAxis."
