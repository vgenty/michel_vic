import ROOT
from ROOT import gStyle, gROOT, TStyle

def fix_hist(j):
    j.GetXaxis().CenterTitle()
    j.GetYaxis().CenterTitle()

def looks():
    gROOT.ProcessLine('TGaxis::SetMaxDigits(3)')
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameFillColor(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetStatColor(0)
    gStyle.SetTitleFont(42,'')
    gStyle.SetLabelFont(42,'x')
    gStyle.SetTitleFont(42,'x')
    gStyle.SetLabelFont(42,'y')
    gStyle.SetTitleFont(42,'y')
    gStyle.SetLabelFont(42,'z')
    gStyle.SetTitleFont(42,'z')
    gStyle.SetLabelSize(0.048,'x')
    gStyle.SetTitleSize(0.048,'x')
    gStyle.SetLabelSize(0.048,'y')
    gStyle.SetTitleSize(0.048,'y')
    gStyle.SetLabelSize(0.048,'z')
    gStyle.SetTitleSize(0.048,'z')
    gStyle.SetTitleSize(0.048,'')
    gStyle.SetTitleAlign(23)
    gStyle.SetTitleX(0.5)
    gStyle.SetTitleBorderSize(0)
    gStyle.SetTitleFillColor(0)
    gStyle.SetTitleStyle(0)
    gStyle.SetOptStat('emr')
    gStyle.SetOptStat(0)
    gStyle.SetStatBorderSize(0)
    gStyle.SetStatFont(42)
    gStyle.SetStatFontSize(0.048)
    gStyle.SetStatY(0.9)
    gStyle.SetStatX(0.86)
    gStyle.SetHistLineColor(ROOT.kBlue+2)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendFillColor(0)
    gStyle.SetFuncWidth(2)
    gStyle.SetFuncColor(2)
    gStyle.SetPadTopMargin(0.08)
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.06)  
    gStyle.SetCanvasDefX(400)
    gStyle.SetCanvasDefY(20)
    gStyle.SetCanvasDefH(420)
    gStyle.SetCanvasDefW(610)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetHistLineWidth(2)
    gStyle.SetTitleOffset(1.16,'y')
    gStyle.SetTitleOffset(1.20,'x')

def looks_minos():
    minosStyle = TStyle("minosStyle", "MINOS Style");

    # Colors
    #set the background color to white
    minosStyle.SetFillColor(10);
    minosStyle.SetFrameFillColor(10);
    minosStyle.SetCanvasColor(10);
    minosStyle.SetPadColor(10);
    minosStyle.SetTitleFillColor(0);
    minosStyle.SetStatColor(10);

    #dont put a colored frame around the plots
    minosStyle.SetFrameBorderMode(0);
    minosStyle.SetCanvasBorderMode(0);
    minosStyle.SetPadBorderMode(0);

    #use the primary color palette
    #minosStyle.SetPalette(1,0);

    #set the default line color for a histogram to be black
    minosStyle.SetHistLineColor(ROOT.kBlack);
    
    #set the default line color for a fit function to be red
    minosStyle.SetFuncColor(ROOT.kRed);

    #make the axis labels black
    minosStyle.SetLabelColor(ROOT.kBlack,"xyz");
    
    #set the default title color to be black
    minosStyle.SetTitleColor(ROOT.kBlack);
    
    # Sizes
    
    #set the margins
    minosStyle.SetPadBottomMargin(0.2);
    minosStyle.SetPadTopMargin(0.075);
    minosStyle.SetPadLeftMargin(0.15);
    
    #set axis label and title text sizes
    minosStyle.SetLabelSize(0.05,"xyz");
    minosStyle.SetTitleSize(0.06,"xyz");
    minosStyle.SetTitleOffset(0.9,"x");
    minosStyle.SetTitleOffset(0.8,"yz");
    minosStyle.SetStatFontSize(0.05);
    minosStyle.SetTextSize(0.06);
    minosStyle.SetTitleBorderSize(0);
    minosStyle.SetStatBorderSize(0);
    
    #set line widths
    minosStyle.SetHistLineWidth(2);
    minosStyle.SetFrameLineWidth(2);
    minosStyle.SetFuncWidth(2);
    
    # Misc

    #align the titles to be centered
    minosStyle.SetTitleAlign(22);
    
    #set the number of divisions to show
    minosStyle.SetNdivisions(506, "xy");

    #turn off xy grids
    minosStyle.SetPadGridX(0);
    minosStyle.SetPadGridY(0);
    
    #set the tick mark style
    #minosStyle.SetPadTickX(1);
    #minosStyle.SetPadTickY(1);
    
    #show the fit parameters in a box
    minosStyle.SetOptFit(1111);
    
    #turn off all other stats
    minosStyle.SetOptStat(0000000);
    
    #marker settings
    minosStyle.SetMarkerStyle(8);
    minosStyle.SetMarkerSize(0.9);
    
    # Fonts
    
    kMinosFont = int(42);
    
    minosStyle.SetStatFont(kMinosFont);
    minosStyle.SetLabelFont(kMinosFont,"xyz");
    minosStyle.SetTitleFont(kMinosFont,"xyz");
    minosStyle.SetTextFont(kMinosFont);

    #done
    minosStyle.cd();

    gROOT.ForceStyle();
    gStyle.ls();


def looks_atlas():
    atlasStyle = TStyle("ATLAS","Atlas style");
    
    # use plain black on white colors
    icol=int(0); # WHITE
    atlasStyle.SetFrameBorderMode(icol);
    atlasStyle.SetFrameFillColor(icol);
    atlasStyle.SetCanvasBorderMode(icol);
    atlasStyle.SetCanvasColor(icol);
    atlasStyle.SetPadBorderMode(icol);
    atlasStyle.SetPadColor(icol);
    atlasStyle.SetStatColor(icol);
    #atlasStyle.SetFillColor(icol); # don't use: white fill color for *all* objects
    
    # set the paper & margin sizes
    atlasStyle.SetPaperSize(20,26);
    
    # set margin sizes
    atlasStyle.SetPadTopMargin(0.05);
    atlasStyle.SetPadRightMargin(0.05);
    atlasStyle.SetPadBottomMargin(0.16);
    atlasStyle.SetPadLeftMargin(0.16);
    
    # set title offsets (for axis label)
    atlasStyle.SetTitleXOffset(1);
    atlasStyle.SetTitleYOffset(1);
    
    # use large fonts
    #Int_t font=72; # Helvetica italics
    font=int(42); # Helvetica
    tsize=int(0.05);
    atlasStyle.SetTextFont(font);
    
    atlasStyle.SetTextSize(tsize);
    atlasStyle.SetLabelFont(font,"x");
    atlasStyle.SetTitleFont(font,"x");
    atlasStyle.SetLabelFont(font,"y");
    atlasStyle.SetTitleFont(font,"y");
    atlasStyle.SetLabelFont(font,"z");
    atlasStyle.SetTitleFont(font,"z");
    
    atlasStyle.SetLabelSize(tsize,"x");
    atlasStyle.SetTitleSize(tsize,"x");
    atlasStyle.SetLabelSize(tsize,"y");
    atlasStyle.SetTitleSize(tsize,"y");
    atlasStyle.SetLabelSize(tsize,"z");
    atlasStyle.SetTitleSize(tsize,"z");
    
    # use bold lines and markers
    atlasStyle.SetMarkerStyle(20);
    atlasStyle.SetMarkerSize(1.2);
    atlasStyle.SetHistLineWidth(2);
    atlasStyle.SetLineStyleString(2,"[12 12]"); # postscript dashes
    
    # get rid of X error bars 
    #atlasStyle.SetErrorX(0.001);
    # get rid of error bar caps
    atlasStyle.SetEndErrorSize(0.);
    
    # do not display any of the standard histogram decorations
    atlasStyle.SetOptTitle(0);
    #atlasStyle.SetOptStat(1111);
    atlasStyle.SetOptStat(0);
    #atlasStyle.SetOptFit(1111)
    atlasStyle.SetOptFit(0);
    
    # put tick marks on top and RHS of plots
    atlasStyle.SetPadTickX(1);
    atlasStyle.SetPadTickY(1);
    atlasStyle.cd();
    gROOT.ForceStyle();
    gStyle.ls();


def looks_atlas2():

    atlasStyle = TStyle("ATLAS","Atlas style");
    icol=0
    atlasStyle.SetFrameBorderMode(icol)
    atlasStyle.SetFrameFillColor(icol)
    atlasStyle.SetCanvasBorderMode(icol)
    atlasStyle.SetCanvasColor(icol)
    atlasStyle.SetPadBorderMode(icol)
    atlasStyle.SetPadColor(icol)
    atlasStyle.SetStatColor(icol)
    atlasStyle.SetPaperSize(20,26)

    atlasStyle.SetPadTopMargin(0.05)
    atlasStyle.SetPadRightMargin(0.05)
    atlasStyle.SetPadBottomMargin(0.16)
    atlasStyle.SetPadLeftMargin(0.16)
    
    atlasStyle.SetTitleXOffset(1.4)
    atlasStyle.SetTitleYOffset(1.4)

    font=42
    tsize=0.05
    atlasStyle.SetTextFont(font)

    atlasStyle.SetTextSize(tsize)
    atlasStyle.SetLabelFont(font,'x')
    atlasStyle.SetTitleFont(font,'x')
    atlasStyle.SetLabelFont(font,'y')
    atlasStyle.SetTitleFont(font,'y')
    atlasStyle.SetLabelFont(font,'z')
    atlasStyle.SetTitleFont(font,'z')

    atlasStyle.SetLabelSize(tsize,'x')
    atlasStyle.SetTitleSize(tsize,'x')
    atlasStyle.SetLabelSize(tsize,'y')
    atlasStyle.SetTitleSize(tsize,'y')
    atlasStyle.SetLabelSize(tsize,'z')
    atlasStyle.SetTitleSize(tsize,'z')

    atlasStyle.SetMarkerStyle(20)
    atlasStyle.SetMarkerSize(1.2)
    atlasStyle.SetHistLineWidth(2)
    atlasStyle.SetLineStyleString(2,'[12 12]')

    atlasStyle.SetEndErrorSize(0.)
    

    atlasStyle.SetOptTitle(0)
    atlasStyle.SetOptStat(0)
    atlasStyle.SetOptFit(0)

    atlasStyle.SetPadTickX(1)
    atlasStyle.SetPadTickY(1)

    atlasStyle.cd()
    gROOT.ForceStyle()
    atlasStyle.ls()
    
    
def looks_CmsTDR():
    tdrStyle = TStyle("tdrStyle","Style for P-TDR");
        
    # For the canvas:
    tdrStyle.SetCanvasBorderMode(0);
    tdrStyle.SetCanvasColor(ROOT.kWhite);
    tdrStyle.SetCanvasDefH(600); #Height of canvas
    tdrStyle.SetCanvasDefW(600); #Width of canvas
    tdrStyle.SetCanvasDefX(0);   #POsition on screen
    tdrStyle.SetCanvasDefY(0);
    
    # For the Pad:
    tdrStyle.SetPadBorderMode(0);
    # tdrStyle.SetPadBorderSize(Width_t size = 1);
    tdrStyle.SetPadColor(ROOT.kWhite);
    #tdrStyle.SetPadGridX(false);
    #tdrStyle.SetPadGridY(false);
    tdrStyle.SetGridColor(0);
    tdrStyle.SetGridStyle(3);
    tdrStyle.SetGridWidth(1);
    
    # For the frame:
    tdrStyle.SetFrameBorderMode(0);
    tdrStyle.SetFrameBorderSize(1);
    tdrStyle.SetFrameFillColor(0);
    tdrStyle.SetFrameFillStyle(0);
    tdrStyle.SetFrameLineColor(1);
    tdrStyle.SetFrameLineStyle(1);
    tdrStyle.SetFrameLineWidth(1);
    
    # For the histo:
    # tdrStyle.SetHistFillColor(1);
    # tdrStyle.SetHistFillStyle(0);
    tdrStyle.SetHistLineColor(1);
    tdrStyle.SetHistLineStyle(0);
    tdrStyle.SetHistLineWidth(1);
    # tdrStyle.SetLegoInnerR(Float_t rad = 0.5);
    # tdrStyle.SetNumberContours(Int_t number = 20);
    
    tdrStyle.SetEndErrorSize(2);
    #tdrStyle.SetErrorMarker(20);  # Seems to give an error
    tdrStyle.SetErrorX(0.);
    
    tdrStyle.SetMarkerStyle(20);
    
    #For the fit/function:
    tdrStyle.SetOptFit(1);
    tdrStyle.SetFitFormat("5.4g");
    tdrStyle.SetFuncColor(2);
    tdrStyle.SetFuncStyle(1);
    tdrStyle.SetFuncWidth(1);
    
    #For the date:
    tdrStyle.SetOptDate(0);
    # tdrStyle.SetDateX(Float_t x = 0.01);
    # tdrStyle.SetDateY(Float_t y = 0.01);
    
    # For the statistics box:
    tdrStyle.SetOptFile(0);
    tdrStyle.SetOptStat(0); # To display the mean and RMS:   SetOptStat("mr");
    tdrStyle.SetStatColor(ROOT.kWhite);
    tdrStyle.SetStatFont(42);
    tdrStyle.SetStatFontSize(0.025);
    tdrStyle.SetStatTextColor(1);
    tdrStyle.SetStatFormat("6.4g");
    tdrStyle.SetStatBorderSize(1);
    tdrStyle.SetStatH(0.1);
    tdrStyle.SetStatW(0.15);
    # tdrStyle.SetStatStyle(Style_t style = 1001);
    # tdrStyle.SetStatX(Float_t x = 0);
    # tdrStyle.SetStatY(Float_t y = 0);
    
    # Margins:
    tdrStyle.SetPadTopMargin(0.05);
    tdrStyle.SetPadBottomMargin(0.13);
    tdrStyle.SetPadLeftMargin(0.16);
    tdrStyle.SetPadRightMargin(0.02);
    
    # For the Global title:
    tdrStyle.SetOptTitle(1);    # 0=No Title
    tdrStyle.SetTitleFont(42);
    tdrStyle.SetTitleColor(1);
    tdrStyle.SetTitleTextColor(1);
    tdrStyle.SetTitleFillColor(10);
    tdrStyle.SetTitleFontSize(0.05);
    # tdrStyle.SetTitleH(0); # Set the height of the title box
    # tdrStyle.SetTitleW(0); # Set the width of the title box
    # tdrStyle.SetTitleX(0); # Set the position of the title box
    # tdrStyle.SetTitleY(0.985); # Set the position of the title box
    # tdrStyle.SetTitleStyle(Style_t style = 1001);
    # tdrStyle.SetTitleBorderSize(2);
    
    # For the axis titles:
    tdrStyle.SetTitleColor(1, "XYZ");
    tdrStyle.SetTitleFont(42, "XYZ");
    tdrStyle.SetTitleSize(0.06, "XYZ");
    # tdrStyle.SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
    # tdrStyle.SetTitleYSize(Float_t size = 0.02);
    tdrStyle.SetTitleXOffset(0.9);
    tdrStyle.SetTitleYOffset(1.25);
    # tdrStyle.SetTitleOffset(1.1, "Y"); # Another way to set the Offset
    
    # For the axis labels:
    tdrStyle.SetLabelColor(1, "XYZ");
    tdrStyle.SetLabelFont(42, "XYZ");
    tdrStyle.SetLabelOffset(0.007, "XYZ");
    tdrStyle.SetLabelSize(0.05, "XYZ");
    
    # For the axis:
    tdrStyle.SetAxisColor(1, "XYZ");
    tdrStyle.SetStripDecimals(ROOT.kTRUE);
    tdrStyle.SetTickLength(0.03, "XYZ");
    tdrStyle.SetNdivisions(510, "XYZ");
    tdrStyle.SetPadTickX(0);  # 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom
    tdrStyle.SetPadTickY(1);
    
    # Change for log plots:
    tdrStyle.SetOptLogx(0);
    tdrStyle.SetOptLogy(0);
    tdrStyle.SetOptLogz(0);
    
    # Postscript options:
    tdrStyle.SetPaperSize(20.,20.);
    # tdrStyle.SetLineScalePS(Float_t scale = 3);
    # tdrStyle.SetLineStyleString(Int_t i, const char* text);
    # tdrStyle.SetHeaderPS(const char* header);
    # tdrStyle.SetTitlePS(const char* pstitle);
    
    # tdrStyle.SetBarOffset(Float_t baroff = 0.5);
    # tdrStyle.SetBarWidth(Float_t barwidth = 0.5);
    # tdrStyle.SetPaintTextFormat(const char* format = "g");
    # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    # tdrStyle.SetTimeOffset(Double_t toffset);
    # tdrStyle.SetHistMinimumZero(kTRUE);
    
    
    tdrStyle.cd();
    gROOT.ForceStyle();  # Try this if stuff doesn't work right
    gStyle.ls();
    
