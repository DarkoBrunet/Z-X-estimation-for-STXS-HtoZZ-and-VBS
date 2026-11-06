import ROOT
import math
import os
from ROOT import TColor
from analyzeZX import setCavasAndStyles, m4l_plot_config, VBS
_EOS_USER = "/eos/user/b/brunetd"
HTOZZ_EOS_DIR = os.path.join(_EOS_USER, "ZplusXpython_HtoZZ2024")
VBS_EOS_DIR = os.path.join(_EOS_USER, "ZpXVBS2024_results")
ANALYSIS_EOS_DIR = VBS_EOS_DIR if VBS else HTOZZ_EOS_DIR
if VBS:
    # Existing VBS archive: hists + estimates in ROOT_FILES/, plots at top level.
    STAGE1_DIR = os.path.join(ANALYSIS_EOS_DIR, "ROOT_FILES")
    STAGE2_DIR = STAGE1_DIR
    CR_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "CR_plots")
    SR_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "SR_plots")
    Z1L_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "Z1L_plots")
    MET_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "mEt_plots")
    FR_PTL3_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "FR_plots")
    FR_MZ1_PLOT_DIR = os.path.join(ANALYSIS_EOS_DIR, "FR_plots")
else:
    STAGE1_DIR = os.path.join(ANALYSIS_EOS_DIR, "stage1")
    STAGE2_DIR = os.path.join(ANALYSIS_EOS_DIR, "stage2")
    PLOT_BASE_DIR = os.path.join(ANALYSIS_EOS_DIR, "plots")
    CR_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "CR_plots")
    SR_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "SR_plots")
    Z1L_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "Z1L_plots")
    MET_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "mEt_plots")
    FR_PTL3_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "FR_ptl3_plots")
    FR_MZ1_PLOT_DIR = os.path.join(PLOT_BASE_DIR, "FR_mZ1_plots")


def _m4l_ytitle(cfg):
    return "Events / %g GeV" % cfg["bin_width"]


def _style_m4l_axis(hist, cfg):
    hist.GetXaxis().SetTitle("m_{4#font[12]{l}} (GeV)")
    hist.GetXaxis().SetRangeUser(cfg["xmin"], cfg["xmax"])
    hist.GetYaxis().SetTitle(_m4l_ytitle(cfg))


def _make_stack_frame(name, cfg):
    frame = ROOT.TH1F(name, "Proba", cfg["nbins"], cfg["xmin"], cfg["xmax"])
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetLineColor(TColor.GetColor("#000099"))
    frame.SetLineStyle(0)
    frame.SetMarkerStyle(20)
    for axis in (frame.GetXaxis(), frame.GetYaxis(), frame.GetZaxis()):
        axis.SetLabelFont(42)
        axis.SetLabelOffset(0.007)
        axis.SetLabelSize(0.05 if axis != frame.GetZaxis() else 0.05)
        axis.SetTitleSize(0.06 if axis != frame.GetZaxis() else 0.06)
        axis.SetTitleFont(42)
    frame.GetXaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.GetZaxis().SetLabelOffset(0.007)
    return frame

def CRPlot():
    m4l_cfg = m4l_plot_config()
    histName = ["2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu","3P1F_2mu2e"]
    StackTitle = "Proba"

    # --- Active: match main_FR_CR.py → Hist_MC_<var>_<Nickname>.root ---
    #ProcessNames = ["Data", "WZ"]
   # fileSize = 2
   # Fillcolors = ["#cc0099"]
   # Linecolors = ["#990066"]

    Fillcolors = ["#99ccff", "#cc0099", "#996666", "#669966"]
    Linecolors = ["#000099", "#990066", "#5f3f3f", "#003300"]
    ProcessNames = ["Data", "ZZ", "WZ", "TT", "DY50"]
    LabelNames = ["Data", "#Z\\gamma^*,ZZ#", "#WZ#", "#t\\bar{t}+jets#", "#Z+jets#"]
    fileSize = 5
    
    plotSize = 10
    variableN = "ptl3"
    for p in range(plotSize):
        print( "Ulazi u plotting" ) 
        c1 = ROOT.TCanvas("c1_"+histName[p]+"_process", "myPlots",0,67,600,600)
        leg = ROOT.TLegend(0.45,0.5,1.05,0.85)
        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
    
        if (p>4):
            entry =leg.AddEntry("NULL","Control Region "+ histName[p] ,"h")
        else:
            entry = leg.AddEntry("NULL","Control Region "+ histName[p] ,"h")

        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        entry=leg.AddEntry("NULL","Data","LP")
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(20)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","Z + jets","F")
        
        ci = TColor.GetColor("#669966")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#003300")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","t#bar{t} + jets","F")
        
        ci = TColor.GetColor("#996666")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#5f3f3f")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","WZ","F")
        
        ci = TColor.GetColor("#cc0099")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#990066")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","Z#gamma*,ZZ","F")
        
        ci = TColor.GetColor("#99ccff")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        if (p>4):

            ci = TColor.GetColor("#000099")
            entry.SetLineColor(ci)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)
            entry=leg.AddEntry("NULL","2P2F contribution","F")
            entry.SetFillStyle(4000)
            entry.SetLineColor(6)
            entry.SetLineStyle(1)
            entry.SetLineWidth(2)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)        
        
        
        
        A = ROOT.THStack(StackTitle,StackTitle)
        
        Proba_stack_5_stack_1 = _make_stack_frame("Proba_stack_5_stack_1", m4l_cfg)
        A.SetHistogram(Proba_stack_5_stack_1)
        
        _file = []
        histPlot = []
        for i in range(fileSize):

            if (i == 0):
                Title = "Data"
            else:
                Title = "MC"

            _file.append(ROOT.TFile("Hist_"+Title+"_"+variableN+"_"+ProcessNames[i]+".root"))
            
            if (i==0) :
                dataPlot = _file[i].Get("h1D_m4l_"+histName[p])
                histPlot.append(dataPlot)

                print( "Data " + histName[p]+ " := "  + str(dataPlot.Integral()) )
                
            
            else:
                histPlot.append(_file[i].Get("h1D_m4l_"+histName[p]))
                histPlot[i].SetFillColor(TColor.GetColor(Fillcolors[i-1]))
                histPlot[i].SetLineColor(TColor.GetColor(Linecolors[i-1]))
                histPlot[i].SetFillStyle(1001)
                
                histPlot[i].SetLineStyle(0)
                histPlot[i].SetMarkerStyle(20)
                histPlot[i].GetXaxis().SetLabelFont(42)
                histPlot[i].GetXaxis().SetLabelOffset(0.007)
                histPlot[i].GetXaxis().SetLabelSize(0.05)
                histPlot[i].GetXaxis().SetTitleSize(0.06)
                histPlot[i].GetXaxis().SetTitleOffset(0.9)
                histPlot[i].GetXaxis().SetTitleFont(42)
                histPlot[i].GetYaxis().SetLabelFont(42)
                histPlot[i].GetYaxis().SetLabelOffset(0.007)
                histPlot[i].GetYaxis().SetLabelSize(0.05)
                histPlot[i].GetYaxis().SetTitleSize(0.06)
                histPlot[i].GetYaxis().SetTitleOffset(1.25)
                histPlot[i].GetYaxis().SetTitleFont(42)
                histPlot[i].GetZaxis().SetLabelFont(42)
                histPlot[i].GetZaxis().SetLabelOffset(0.007)
                histPlot[i].GetZaxis().SetLabelSize(0.05)
                histPlot[i].GetZaxis().SetTitleSize(0.06)
                histPlot[i].GetZaxis().SetTitleFont(42)
                
                print( "MC: "+ histName[p] + " := "+ str( histPlot[i].Integral() ) )
                A.Add(histPlot[i])

        print ("Data status " + str(dataPlot))    
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        
        dataPlot.SetFillColor(TColor.GetColor("#000000"))
        dataPlot.SetLineColor(TColor.GetColor("#000000"))
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetRange(1, dataPlot.GetNbinsX())
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        _style_m4l_axis(dataPlot, m4l_cfg)
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.25)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(42)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        
        dataPlot.GetYaxis().SetRangeUser(0, 1.3 * dataPlot.GetMaximum())
        dataPlot.Draw("e1 goff")
        
        A.Draw("hist same")
        dataPlot.Draw("same e1 goff")
        if ( p>4 ):
            _file2p2f = ROOT.TFile(_resolve_estimate_file("estimateZXData"), "READ")
            histPlot2p2f = _file2p2f.Get("h1D_m4l_Add_"+histName[p-5])
            
            print("2P2F contr:= " + str(histPlot2p2f.Integral()))
            
            histPlot2p2f.Add(histPlot[1])  # ZZ
            histPlot2p2f.Add(histPlot[2])  # WZ
            histPlot2p2f.SetFillColor(TColor.GetColor("#ffffff"))
            histPlot2p2f.SetLineColor(TColor.GetColor("#ff00ff"))
            histPlot2p2f.SetLineWidth(2)
            histPlot2p2f.SetFillStyle(4000)
            histPlot2p2f.Smooth()
            histPlot2p2f.Draw("hist E1 same goff")

        leg.Draw("goff")
        
        tex1 = ROOT.TLatex(0.95,0.96,"109.82 fb^{-1} (13.6 TeV)")
        tex1.SetNDC()
        tex1.SetTextAlign(31)
        tex1.SetTextFont(42)
        tex1.SetTextSize(0.03)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.15,0.96,"CMS")
        tex2.SetNDC()
        tex2.SetTextFont(61)
        tex2.SetTextSize(0.0375)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        tex3 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex3.SetNDC()
        tex3.SetTextFont(52)
        tex3.SetTextSize(0.0285)
        tex3.SetLineWidth(2)
        tex3.Draw("goff")
        
        c1.SaveAs(CR_PLOT_DIR + "/" + histName[p] + "_auto.pdf")
        c1.SaveAs(CR_PLOT_DIR + "/m" + histName[p] + "_auto.C")
        c1.SaveAs(CR_PLOT_DIR + "/" + histName[p] + "_auto.root")



def CRMatchPlot():
    histName = ["2P2F","2P2F4e","2P2F4mu","2P2F2e2mu","2P2F2mu2e","3P1F","3P1F4e","3P1F4mu","3P1F2e2mu", "3P1F2mu2e"]
    dataName = ["2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu", "3P1F_2mu2e"]
    Fillcolors = ["#AF0096","#91FFFF","#0000ff","#669965","#3CB371"]
    Linecolors = ["#AF0096","#91FFFF","#0000ff","#669965","#3CB371"]
    StackTitle = "Proba"
    # (ProcessNames / fileSize are not used in CRMatchPlot loops; kept as documentation.)
    ProcessNames = ["Data", "WZ"]
    fileSize = 2
    # ProcessNames = ["Data", "WZ", "ZZ", "TT", "DY50"]
    # fileSize = 5
    plotSize = 10
    variableN = "ptl3"
    Title = "Data"
    _file1 = ROOT.TFile("Hist_"+Title+"_"+variableN+"_"+Title+".root", "READ")
    for p in range(plotSize):
        # Canvas first: TLegend dummy "NULL" entries are dropped if no pad exists.
        c1 = ROOT.TCanvas("c1_"+histName[p]+"_match", "myPlots",0,67,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)

        leg_xl = 0.45
        leg_xr = 0.95
        leg_yt = 0.85
        leg_yb = 0.48 if p > 4 else 0.52
        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        entry = leg.AddEntry("NULL", "Control Region " + dataName[p], "h")
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        A = ROOT.THStack(StackTitle,StackTitle)

        dataPlot = _file1.Get("h1D_m4l_"+dataName[p])
        dataPlot.SetDirectory(0)
        m4l_cfg = {
            "nbins": dataPlot.GetNbinsX(),
            "xmin": dataPlot.GetXaxis().GetXmin(),
            "xmax": dataPlot.GetXaxis().GetXmax(),
            "bin_width": dataPlot.GetBinWidth(1),
        }
        leg.AddEntry(dataPlot,Title,"LP")
        
        _fileZZ = ROOT.TFile("Hist_ID_"+variableN+"_ZZ.root", "READ")
        _fileWZ = ROOT.TFile("Hist_ID_"+variableN+"_WZ.root", "READ")
        _fileTT = ROOT.TFile("Hist_ID_"+variableN+"_TT.root", "READ")
        _fileDY = ROOT.TFile("Hist_ID_"+variableN+"_DY50.root", "READ")

        histPlotC = _fileZZ.Get("h1D_m4l_Conv_"+histName[p]).Clone("h1D_m4l_Conv_"+histName[p]+"_sum")
        histPlotC.SetDirectory(0)
        histPlotC.Add(_fileWZ.Get("h1D_m4l_Conv_"+histName[p]))
        histPlotC.Add(_fileTT.Get("h1D_m4l_Conv_"+histName[p]))
        histPlotC.Add(_fileDY.Get("h1D_m4l_Conv_"+histName[p]))
        histPlotC.SetFillColor(30)
        histPlotC.SetLineColor(1)
        histPlotC.SetFillStyle(1001)
        
        histPlotC.SetLineStyle(0)
        histPlotC.SetMarkerStyle(20)
        
        leg.AddEntry(histPlotC,"Conversions","F")
        A.Add(histPlotC)
        
        histPlotP = _fileZZ.Get("h1D_m4l_Prom_"+histName[p]).Clone("h1D_m4l_Prom_"+histName[p]+"_sum")
        histPlotP.SetDirectory(0)
        histPlotP.Add(_fileWZ.Get("h1D_m4l_Prom_"+histName[p]))
        histPlotP.Add(_fileTT.Get("h1D_m4l_Prom_"+histName[p]))
        histPlotP.Add(_fileDY.Get("h1D_m4l_Prom_"+histName[p]))
        histPlotP.SetFillColor(38)
        histPlotP.SetLineColor(1)
        histPlotP.SetFillStyle(1001)
        
        histPlotP.SetLineStyle(0)
        histPlotP.SetMarkerStyle(20)
        leg.AddEntry(histPlotP,"Prompt","F")
        A.Add(histPlotP)
        
        histPlotF = _fileZZ.Get("h1D_m4l_Fake_"+histName[p]).Clone("h1D_m4l_Fake_"+histName[p]+"_sum")
        histPlotF.SetDirectory(0)
        histPlotF.Add(_fileWZ.Get("h1D_m4l_Fake_"+histName[p]))
        histPlotF.Add(_fileTT.Get("h1D_m4l_Fake_"+histName[p]))
        histPlotF.Add(_fileDY.Get("h1D_m4l_Fake_"+histName[p]))
        histPlotF.SetFillColor(46)
        histPlotF.SetLineColor(1)
        histPlotF.SetFillStyle(1001)
        
        histPlotF.SetLineStyle(0)
        histPlotF.SetMarkerStyle(20)
        leg.AddEntry(histPlotF,"LF Fakes","F")
        A.Add(histPlotF)
        
        histPlotBD = _fileZZ.Get("h1D_m4l_BDfake_"+histName[p]).Clone("h1D_m4l_BDfake_"+histName[p]+"_sum")
        histPlotBD.SetDirectory(0)
        histPlotBD.Add(_fileWZ.Get("h1D_m4l_BDfake_"+histName[p]))
        histPlotBD.Add(_fileTT.Get("h1D_m4l_BDfake_"+histName[p]))
        histPlotBD.Add(_fileDY.Get("h1D_m4l_BDfake_"+histName[p]))
        histPlotBD.SetFillColor(49)
        histPlotBD.SetLineColor(1)
        histPlotBD.SetFillStyle(1001)
        histPlotBD.SetLineStyle(0)
        histPlotBD.SetMarkerStyle(20)
        leg.AddEntry(histPlotBD,"HF fakes","F")
        A.Add(histPlotBD)

        dataPlot.SetFillColor(1)
        dataPlot.SetLineColor(1)
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetRange(1, dataPlot.GetNbinsX())
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        _style_m4l_axis(dataPlot, m4l_cfg)
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.25)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(42)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        
        dataPlot.GetYaxis().SetRangeUser(0, 1.3 * dataPlot.GetMaximum())
        dataPlot.Draw("e1 goff")
        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        if (p > 4):
        
            _fileD = ROOT.TFile(_resolve_estimate_file("estimateZXData"), "READ")
            histPlot1 = _fileD.Get("h1D_m4l_Add_"+dataName[p-5]).Clone("h1D_m4l_Add_"+dataName[p-5]+"_match")
            histPlot1.SetDirectory(0)
            print ("h1D_m4l_Add_"+dataName[p-5])
            print (histPlot1.Integral())
            leg.AddEntry(histPlot1,"2P2F contr smooth","F")
            histPlot1.SetFillColor(0)
            histPlot1.SetLineColor(6)
            histPlot1.SetLineWidth(2)
            histPlot1.SetFillStyle(4000)
            histPlot1.Smooth()
            histPlot1.Draw("hist same goff")
            _fileD.Close()
            

        leg.Draw("goff")
        
        tex1 = ROOT.TLatex(0.95,0.96,"109.82 fb^{-1} (13.6 TeV)")
        tex1.SetNDC()
        tex1.SetTextAlign(31)
        tex1.SetTextFont(42)
        tex1.SetTextSize(0.03)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.15,0.96,"CMS")
        tex2.SetNDC()
        tex2.SetTextFont(61)
        tex2.SetTextSize(0.0375)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        tex3 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex3.SetNDC()
        tex3.SetTextFont(52)
        tex3.SetTextSize(0.0285)
        tex3.SetLineWidth(2)
        tex3.Draw("goff")
        
        c1.SaveAs(CR_PLOT_DIR + "/" + histName[p] + "_match_auto.pdf")
        c1.SaveAs(CR_PLOT_DIR + "/" + histName[p] + "_match_auto.C")

        _fileZZ.Close()
        _fileWZ.Close()
        _fileTT.Close()
        _fileDY.Close()
    _file1.Close()


def SRMatchPlot():
    """Truth-origin decomposition for 4P signal-region candidates.

    Uses the same Hist_ID naming convention as CRMatchPlot:
    h1D_m4l_{Prom,Fake,BDfake,Conv}_4P...
    """
    m4l_cfg = m4l_plot_config()
    histName = ["4P", "4P4e", "4P4mu", "4P2e2mu", "4P2mu2e"]
    dataName = ["4P", "4P_4e", "4P_4mu", "4P_2e2mu", "4P_2mu2e"]
    plotSize = len(histName)
    variableN = "ptl3"
    Title = "Data"

    _file1 = ROOT.TFile("Hist_" + Title + "_" + variableN + "_" + Title + ".root")
    _fileZZ = ROOT.TFile("Hist_ID_" + variableN + "_ZZ.root")
    _fileWZ = ROOT.TFile("Hist_ID_" + variableN + "_WZ.root")
    _fileTT = ROOT.TFile("Hist_ID_" + variableN + "_TT.root")
    _fileDY = ROOT.TFile("Hist_ID_" + variableN + "_DY50.root")
    mc_files = [_fileZZ, _fileWZ, _fileTT, _fileDY]

    def _sum_origin(origin, section):
        first = mc_files[0].Get("h1D_m4l_" + origin + "_" + section)
        if not first:
            return None
        hist = first.Clone("h1D_m4l_" + origin + "_" + section + "_4p_sum")
        hist.SetDirectory(0)
        for f in mc_files[1:]:
            h = f.Get("h1D_m4l_" + origin + "_" + section)
            if not h:
                return None
            hist.Add(h)
        return hist

    for p in range(plotSize):
        dataPlot = _file1.Get("h1D_m4l_" + dataName[p])
        histPlotC = _sum_origin("Conv", histName[p])
        histPlotP = _sum_origin("Prom", histName[p])
        histPlotF = _sum_origin("Fake", histName[p])
        histPlotBD = _sum_origin("BDfake", histName[p])
        if not dataPlot or not histPlotC or not histPlotP or not histPlotF or not histPlotBD:
            print("Skipping 4P match plot for %s: rerun stage 1 to produce Hist_ID 4P histograms" % histName[p])
            continue

        leg = ROOT.TLegend(0.45, 0.55, 0.95, 0.85)
        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)

        stack = ROOT.THStack("Proba", "Proba")
        for hist, label, color in (
            (histPlotC, "Conversions", 30),
            (histPlotP, "Prompt", 38),
            (histPlotF, "LF Fakes", 46),
            (histPlotBD, "HF fakes", 49),
        ):
            hist.SetFillColor(color)
            hist.SetLineColor(1)
            hist.SetFillStyle(1001)
            hist.SetLineStyle(0)
            hist.SetMarkerStyle(20)
            leg.AddEntry(hist, label, "F")
            stack.Add(hist)

        c1 = ROOT.TCanvas("c1_" + histName[p] + "_4p_match", "myPlots", 0, 67, 600, 600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)

        dataPlot.SetFillColor(1)
        dataPlot.SetLineColor(1)
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetRange(1, dataPlot.GetNbinsX())
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        _style_m4l_axis(dataPlot, m4l_cfg)
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.25)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetYaxis().SetRangeUser(0, 1.3 * max(dataPlot.GetMaximum(), 1.0))

        dataPlot.Draw("e1 goff")
        stack.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        leg.AddEntry(dataPlot, Title, "LP")
        leg.Draw("goff")

        tex1 = ROOT.TLatex(0.95, 0.96, "109.82 fb^{-1} (13.6 TeV)")
        tex1.SetNDC(); tex1.SetTextAlign(31); tex1.SetTextFont(42); tex1.SetTextSize(0.03); tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.15, 0.96, "CMS")
        tex2.SetNDC(); tex2.SetTextFont(61); tex2.SetTextSize(0.0375); tex2.SetLineWidth(2)
        tex2.Draw("goff")
        tex3 = ROOT.TLatex(0.23, 0.96, "Preliminary")
        tex3.SetNDC(); tex3.SetTextFont(52); tex3.SetTextSize(0.0285); tex3.SetLineWidth(2)
        tex3.Draw("goff")

        c1.SaveAs(SR_PLOT_DIR + "/" + histName[p] + "_match_auto.pdf")
        c1.SaveAs(SR_PLOT_DIR + "/" + histName[p] + "_match_auto.C")
        c1.SaveAs(SR_PLOT_DIR + "/" + histName[p] + "_match_auto.root")

    _file1.Close()
    _fileZZ.Close()
    _fileWZ.Close()
    _fileTT.Close()
    _fileDY.Close()


def _resolve_estimate_file(basename):
    for filename in (basename, basename + ".root"):
        path = os.path.join(STAGE2_DIR, filename)
        if os.path.isfile(path):
            return path
    raise RuntimeError("Missing estimate file: %s" % basename)


# Fake-rate variants for SR Z+X estimate plots (stage-2 output tags).
FR_SR_TAGS = ("nominal", "ge1j", "ge2j") if VBS else ("nominal",)
FR_SR_METHOD_LABELS = {
    "nominal": "",
    "ge1j": "#geq1 jet FR",
    "ge2j": "#geq2 jet FR",
}
FR_SR_ZX_LABELS = {
    "nominal": "Z+X estimate",
    "ge1j": "Z+X estimate (#geq1 jet FR)",
    "ge2j": "Z+X estimate (#geq2 jet FR)",
}
FR_SR_TOTAL_LABELS = {
    "nominal": "Estimated from data",
    "ge1j": "Estimated from data (#geq1 jet FR)",
    "ge2j": "Estimated from data (#geq2 jet FR)",
}


def _fr_sr_file_suffix(fr_tag):
    return "" if fr_tag == "nominal" else "_%s" % fr_tag


def _fr_sr_estimate_basenames(fr_tag, var_name="ptl3"):
    if var_name == "ptl3":
        data_base = "estimateZXData"
        zz_base = "estimateZXZZ"
    else:
        data_base = "estimateZXData_%s" % var_name
        zz_base = "estimateZXZZ_%s" % var_name
    if fr_tag != "nominal":
        data_base = "%s_%s" % (data_base, fr_tag)
        zz_base = "%s_%s" % (zz_base, fr_tag)
    return data_base, zz_base


def _fr_sr_estimate_available(fr_tag, var_name="ptl3"):
    data_base, zz_base = _fr_sr_estimate_basenames(fr_tag, var_name)
    for basename in (data_base, zz_base):
        if not (
            os.path.isfile(os.path.join(STAGE2_DIR, basename))
            or os.path.isfile(os.path.join(STAGE2_DIR, basename + ".root"))
        ):
            return False
    return True


def _build_sr_tot_histograms(estimate_basename="estimateZXData", zz_basename="estimateZXZZ"):
    """SR_tot = SR_3P1F - SR_2P2F (bin-by-bin), same logic as estimate_final_numbers_macro.py."""
    estimate_path = _resolve_estimate_file(estimate_basename)
    zz_path = _resolve_estimate_file(zz_basename)

    f_est = ROOT.TFile.Open(estimate_path, "READ")
    f_zz = ROOT.TFile.Open(zz_path, "READ")

    flavor_suffixes = ["_4e", "_4mu", "_2e2mu", "_2mu2e"]
    sr_3p1f = {}
    sr_2p2f = {}
    sr_tot = {}

    for suffix in [""] + flavor_suffixes:
        h3 = f_est.Get("h1D_m4l_SR_3P1F" + suffix)
        h2 = f_est.Get("h1D_m4l_SR_2P2F" + suffix)
        if h3 is None or h2 is None:
            raise RuntimeError("Missing SR histograms for suffix '%s'" % suffix)
        sr_3p1f[suffix] = h3.Clone("sr3p1f" + suffix)
        sr_3p1f[suffix].SetDirectory(0)
        sr_2p2f[suffix] = h2.Clone("sr2p2f" + suffix)
        sr_2p2f[suffix].SetDirectory(0)
        sr_tot[suffix] = sr_2p2f[suffix].Clone("sr_tot" + suffix)
        sr_tot[suffix].SetDirectory(0)

        hzz = f_zz.Get("h1D_m4l_SR_3P1F" + suffix)
        if hzz is not None:
            sr_3p1f[suffix].Add(hzz, -1)

    var_nBins = sr_2p2f[""].GetNbinsX()
    for suffix in flavor_suffixes:
        h3 = sr_3p1f[suffix]
        h2 = sr_2p2f[suffix]
        htot = sr_tot[suffix]
        for c in range(1, var_nBins + 1):
            if h3.GetBinContent(c) <= 2 * h2.GetBinContent(c):
                htot.SetBinContent(c, h2.GetBinContent(c))
                htot.SetBinError(c, h2.GetBinError(c))
            else:
                htot.SetBinContent(c, h3.GetBinContent(c) - h2.GetBinContent(c))
                htot.SetBinError(
                    c,
                    math.sqrt(h2.GetBinError(c) ** 2 + h3.GetBinError(c) ** 2),
                )

    sr_tot[""].Reset()
    for suffix in flavor_suffixes:
        sr_tot[""].Add(sr_tot[suffix])

    f_est.Close()
    f_zz.Close()
    return sr_tot


def _build_sr_input_histograms(estimate_basename="estimateZXData", zz_basename="estimateZXZZ"):
    """Return weighted SR 3P1F and 2P2F inputs used to build the final Z+X estimate."""
    estimate_path = _resolve_estimate_file(estimate_basename)
    zz_path = _resolve_estimate_file(zz_basename)

    f_est = ROOT.TFile.Open(estimate_path, "READ")
    f_zz = ROOT.TFile.Open(zz_path, "READ")

    flavor_suffixes = ["", "_4e", "_4mu", "_2e2mu", "_2mu2e"]
    sr_3p1f = {}
    sr_2p2f = {}

    for suffix in flavor_suffixes:
        h3 = f_est.Get("h1D_m4l_SR_3P1F" + suffix)
        h2 = f_est.Get("h1D_m4l_SR_2P2F" + suffix)
        if h3 is None or h2 is None:
            raise RuntimeError("Missing SR input histograms for suffix '%s'" % suffix)
        sr_3p1f[suffix] = h3.Clone("sr3p1f_input" + suffix)
        sr_3p1f[suffix].SetDirectory(0)
        sr_2p2f[suffix] = h2.Clone("sr2p2f_input" + suffix)
        sr_2p2f[suffix].SetDirectory(0)

        hzz = f_zz.Get("h1D_m4l_SR_3P1F" + suffix)
        if hzz is not None:
            sr_3p1f[suffix].Add(hzz, -1)

    f_est.Close()
    f_zz.Close()
    return sr_3p1f, sr_2p2f


def SRInputPlot():
    """Overlay the weighted SR 3P1F and 2P2F inputs before the final bin-by-bin prescription."""
    m4l_cfg = m4l_plot_config()
    plot_specs = [
        ("SR_inputs", ""),
        ("SR_inputs_4e", "_4e"),
        ("SR_inputs_4mu", "_4mu"),
        ("SR_inputs_2e2mu", "_2e2mu"),
        ("SR_inputs_2mu2e", "_2mu2e"),
    ]

    for fr_tag in FR_SR_TAGS:
        if not _fr_sr_estimate_available(fr_tag):
            print("SRInputPlot: skipping %s, missing estimateZX* files" % fr_tag)
            continue
        data_base, zz_base = _fr_sr_estimate_basenames(fr_tag)
        fr_suffix = _fr_sr_file_suffix(fr_tag)
        method_label = FR_SR_METHOD_LABELS[fr_tag]
        sr_3p1f, sr_2p2f = _build_sr_input_histograms(data_base, zz_base)

        for plot_label, suffix in plot_specs:
            out_label = plot_label + fr_suffix
            print("SRInputPlot: " + out_label)
            h3 = sr_3p1f[suffix].Clone("plot_3p1f_" + out_label)
            h2 = sr_2p2f[suffix].Clone("plot_2p2f_" + out_label)
            h3.SetDirectory(0)
            h2.SetDirectory(0)

            c1 = ROOT.TCanvas("c1_" + out_label, "myPlots", 0, 67, 600, 600)
            leg = ROOT.TLegend(0.42, 0.55, 0.95, 0.85)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.032)

            flavor = plot_label.replace("SR_inputs", "").strip("_")
            title = "SR weighted inputs"
            if flavor:
                title += " " + flavor
            if method_label:
                title += " (" + method_label + ")"
            entry = leg.AddEntry("NULL", title, "h")
            entry.SetTextFont(62)

            ROOT.gStyle.SetOptFit(1)
            ROOT.gStyle.SetOptStat(0)
            ROOT.gStyle.SetOptTitle(0)
            c1.SetFillColor(0)
            c1.SetBorderMode(0)
            c1.SetBorderSize(2)
            c1.SetTickx(1)
            c1.SetTicky(1)
            c1.SetLeftMargin(0.15)
            c1.SetRightMargin(0.05)
            c1.SetTopMargin(0.05)
            c1.SetBottomMargin(0.13)
            c1.SetFrameFillStyle(0)
            c1.SetFrameBorderMode(0)

            h3.SetLineColor(TColor.GetColor("#D55E00"))
            h3.SetMarkerColor(TColor.GetColor("#D55E00"))
            h3.SetMarkerStyle(20)
            h3.SetLineWidth(2)
            h3.SetFillStyle(0)
            h2.SetLineColor(TColor.GetColor("#0072B2"))
            h2.SetMarkerColor(TColor.GetColor("#0072B2"))
            h2.SetMarkerStyle(24)
            h2.SetLineWidth(2)
            h2.SetFillStyle(0)

            _style_m4l_axis(h3, m4l_cfg)
            h3.GetYaxis().SetLabelFont(42)
            h3.GetYaxis().SetLabelSize(0.05)
            h3.GetYaxis().SetTitleSize(0.06)
            h3.GetYaxis().SetTitleOffset(1.25)
            h3.GetYaxis().SetTitleFont(42)
            h3.GetYaxis().SetRangeUser(0, 1.3 * max(h3.GetMaximum(), h2.GetMaximum(), 1.0))

            h3.Draw("e1 goff")
            h2.Draw("same e1 goff")
            leg.AddEntry(h3, "3P1F - ZZ prompt", "lep")
            leg.AddEntry(h2, "2P2F", "lep")
            leg.Draw("goff")

            tex1 = ROOT.TLatex(0.95, 0.96, "109.82 fb^{-1} (13.6 TeV)")
            tex1.SetNDC(); tex1.SetTextAlign(31); tex1.SetTextFont(42); tex1.SetTextSize(0.03); tex1.SetLineWidth(2)
            tex1.Draw("goff")
            tex2 = ROOT.TLatex(0.15, 0.96, "CMS")
            tex2.SetNDC(); tex2.SetTextFont(61); tex2.SetTextSize(0.0375); tex2.SetLineWidth(2)
            tex2.Draw("goff")
            tex3 = ROOT.TLatex(0.23, 0.96, "Preliminary")
            tex3.SetNDC(); tex3.SetTextFont(52); tex3.SetTextSize(0.0285); tex3.SetLineWidth(2)
            tex3.Draw("goff")

            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.pdf")
            c1.SaveAs(SR_PLOT_DIR + "/m" + out_label + "_auto.C")
            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.root")


def SRPlot():
    """4P signal-region data vs prompt MC + final Z+X estimate (SR_3P1F - SR_2P2F)."""
    m4l_cfg = m4l_plot_config()
    plot_specs = [
        ("4P", "4P", ""),
        ("4P4e", "4P_4e", "_4e"),
        ("4P4mu", "4P_4mu", "_4mu"),
        ("4P2e2mu", "4P_2e2mu", "_2e2mu"),
        ("4P2mu2e", "4P_2mu2e", "_2mu2e"),
    ]

    extra_processes = [] if VBS else ["ggZZ", "HZZ"]
    Fillcolors = ["#4E79A7", "#76B7B2", "#E15759"]
    Linecolors = ["#1F4E79", "#2F6F6A", "#9C2F31"]
    ProcessNames = ["Data", "ZZ"] + extra_processes
    ProcessLabels = {
        "ZZ": "Z#gamma*,ZZ",
        "ggZZ": "ggZZ",
        "HZZ": "H#rightarrow ZZ",
    }
    fileSize = len(ProcessNames)
    variableN = "ptl3"
    StackTitle = "Proba"

    for fr_tag in FR_SR_TAGS:
        if not _fr_sr_estimate_available(fr_tag):
            print("SRPlot: skipping %s, missing estimateZX* files" % fr_tag)
            continue
        data_base, zz_base = _fr_sr_estimate_basenames(fr_tag)
        fr_suffix = _fr_sr_file_suffix(fr_tag)
        zx_label = FR_SR_ZX_LABELS[fr_tag]
        sr_tot = _build_sr_tot_histograms(data_base, zz_base)

        for plot_label, data_suffix, sr_suffix in plot_specs:
            out_label = plot_label + fr_suffix
            print("SRPlot: " + out_label)

            c1 = ROOT.TCanvas("c1_" + out_label + "_SR", "myPlots", 0, 67, 600, 600)
            leg = ROOT.TLegend(0.42, 0.48, 1.05, 0.85)
            leg.SetBorderSize(0)
            leg.SetLineColor(1)
            leg.SetLineStyle(1)
            leg.SetLineWidth(1)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)

            entry = leg.AddEntry("NULL", "Signal Region " + plot_label, "h")
            entry.SetLineColor(1)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(62)
            entry = leg.AddEntry("NULL", "Data", "LP")
            entry.SetLineColor(1)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(20)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)
            entry = leg.AddEntry("NULL", zx_label, "F")
            ci = TColor.GetColor("#59A14F")
            entry.SetFillColor(ci)
            entry.SetFillStyle(1001)
            ci = TColor.GetColor("#2E6B2E")
            entry.SetLineColor(ci)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)
            for i, process in enumerate(ProcessNames[1:]):
                entry = leg.AddEntry("NULL", ProcessLabels.get(process, process), "F")
                ci = TColor.GetColor(Fillcolors[i])
                entry.SetFillColor(ci)
                entry.SetFillStyle(1001)
                ci = TColor.GetColor(Linecolors[i])
                entry.SetLineColor(ci)
                entry.SetLineStyle(1)
                entry.SetLineWidth(1)
                entry.SetMarkerColor(1)
                entry.SetMarkerStyle(21)
                entry.SetMarkerSize(1)
                entry.SetTextFont(42)
            A = ROOT.THStack(StackTitle, StackTitle)
            Proba_stack_5_stack_1 = _make_stack_frame(
                "Proba_stack_5_stack_1_SR_" + out_label, m4l_cfg
            )
            A.SetHistogram(Proba_stack_5_stack_1)

            data_name = "h1D_m4l_" + data_suffix
            _file_data = ROOT.TFile("Hist_Data_" + variableN + "_Data.root")
            dataPlot = _file_data.Get(data_name)
            if not dataPlot:
                print(
                    "Skipping %s: missing %s in Data (re-run stage-1 merge after analyzeZX 4P write fix)"
                    % (out_label, data_name)
                )
                _file_data.Close()
                continue
            dataPlot = dataPlot.Clone("data_" + out_label)
            dataPlot.SetDirectory(0)

            histPlot = []
            for i in range(1, fileSize):
                mc_path = "Hist_MC_" + variableN + "_" + ProcessNames[i] + ".root"
                if not os.path.isfile(mc_path):
                    print(
                        "Skipping MC %s for %s: missing %s"
                        % (ProcessNames[i], out_label, mc_path)
                    )
                    continue
                _file = ROOT.TFile(mc_path)
                if not _file or _file.IsZombie():
                    print(
                        "Skipping MC %s for %s: cannot open %s"
                        % (ProcessNames[i], out_label, mc_path)
                    )
                    continue
                hi = _file.Get(data_name)
                if not hi:
                    print(
                        "Skipping MC %s for %s: missing %s"
                        % (ProcessNames[i], out_label, data_name)
                    )
                    _file.Close()
                    continue
                hi = hi.Clone("mc_" + ProcessNames[i] + "_" + out_label)
                hi.SetDirectory(0)
                histPlot.append(hi)
                histPlot[-1].SetFillColor(TColor.GetColor(Fillcolors[i - 1]))
                histPlot[-1].SetLineColor(TColor.GetColor(Linecolors[i - 1]))
                histPlot[-1].SetFillStyle(1001)
                histPlot[-1].SetLineStyle(0)
                histPlot[-1].SetMarkerStyle(20)
                histPlot[-1].GetXaxis().SetLabelFont(42)
                histPlot[-1].GetXaxis().SetLabelOffset(0.007)
                histPlot[-1].GetXaxis().SetLabelSize(0.05)
                histPlot[-1].GetXaxis().SetTitleSize(0.06)
                histPlot[-1].GetXaxis().SetTitleOffset(0.9)
                histPlot[-1].GetXaxis().SetTitleFont(42)
                histPlot[-1].GetYaxis().SetLabelFont(42)
                histPlot[-1].GetYaxis().SetLabelOffset(0.007)
                histPlot[-1].GetYaxis().SetLabelSize(0.05)
                histPlot[-1].GetYaxis().SetTitleSize(0.06)
                histPlot[-1].GetYaxis().SetTitleOffset(1.25)
                histPlot[-1].GetYaxis().SetTitleFont(42)
                histPlot[-1].GetZaxis().SetLabelFont(42)
                histPlot[-1].GetZaxis().SetLabelOffset(0.007)
                histPlot[-1].GetZaxis().SetLabelSize(0.05)
                histPlot[-1].GetZaxis().SetTitleSize(0.06)
                histPlot[-1].GetZaxis().SetTitleFont(42)
                print(
                    "MC 4P "
                    + out_label
                    + " "
                    + ProcessNames[i]
                    + " := "
                    + str(histPlot[-1].Integral())
                )
                A.Add(histPlot[-1])
                _file.Close()

            zxPlot = sr_tot[sr_suffix].Clone("zx_" + out_label)
            zxPlot.SetDirectory(0)
            print("Z+X estimate " + out_label + " := " + str(zxPlot.Integral()))
            zxPlot.SetFillColor(TColor.GetColor("#59A14F"))
            zxPlot.SetLineColor(TColor.GetColor("#2E6B2E"))
            zxPlot.SetLineWidth(1)
            zxPlot.SetFillStyle(1001)
            zxPlot.SetLineStyle(0)
            zxPlot.SetMarkerStyle(20)
            zxPlot.SetMarkerColor(TColor.GetColor("#2E6B2E"))
            A.Add(zxPlot)

            totalBkg = zxPlot.Clone("total_bkg_" + out_label)
            totalBkg.SetDirectory(0)
            for hist in histPlot:
                totalBkg.Add(hist)

            ROOT.gStyle.SetOptFit(1)
            ROOT.gStyle.SetOptStat(0)
            ROOT.gStyle.SetOptTitle(0)
            c1.SetFillColor(0)
            c1.SetBorderMode(0)
            c1.SetBorderSize(2)
            c1.SetTickx(1)
            c1.SetTicky(1)
            c1.SetLeftMargin(0.15)
            c1.SetRightMargin(0.05)
            c1.SetTopMargin(0.05)
            c1.SetBottomMargin(0.13)
            c1.SetFrameFillStyle(0)
            c1.SetFrameBorderMode(0)
            c1.SetFrameFillStyle(0)
            c1.SetFrameBorderMode(0)

            dataPlot.SetFillColor(TColor.GetColor("#000000"))
            dataPlot.SetLineColor(TColor.GetColor("#000000"))
            dataPlot.SetFillStyle(0)
            dataPlot.SetMarkerStyle(20)
            dataPlot.GetXaxis().SetRange(1, dataPlot.GetNbinsX())
            dataPlot.GetXaxis().SetLabelFont(42)
            dataPlot.GetXaxis().SetLabelSize(0.05)
            dataPlot.GetXaxis().SetTitleSize(0.06)
            dataPlot.GetXaxis().SetTitleOffset(0.9)
            dataPlot.GetXaxis().SetTitleFont(42)
            _style_m4l_axis(dataPlot, m4l_cfg)
            dataPlot.GetYaxis().SetLabelFont(42)
            dataPlot.GetYaxis().SetLabelSize(0.05)
            dataPlot.GetYaxis().SetTitleSize(0.06)
            dataPlot.GetYaxis().SetTitleOffset(1.25)
            dataPlot.GetYaxis().SetTitleFont(42)
            dataPlot.GetZaxis().SetLabelFont(42)
            dataPlot.GetZaxis().SetLabelSize(0.035)
            dataPlot.GetZaxis().SetTitleSize(0.035)
            dataPlot.GetZaxis().SetTitleFont(42)

            ymax = max(
                dataPlot.GetMaximum(),
                totalBkg.GetMaximum(),
                1.0,
            )
            dataPlot.GetYaxis().SetRangeUser(0, 1.3 * ymax)
            dataPlot.Draw("e1 goff")

            A.Draw("hist same")
            dataPlot.Draw("same e1 goff")

            leg.Draw("goff")

            tex1 = ROOT.TLatex(0.95, 0.96, "109.82 fb^{-1} (13.6 TeV)")
            tex1.SetNDC()
            tex1.SetTextAlign(31)
            tex1.SetTextFont(42)
            tex1.SetTextSize(0.03)
            tex1.SetLineWidth(2)
            tex1.Draw("goff")
            tex2 = ROOT.TLatex(0.15, 0.96, "CMS")
            tex2.SetNDC()
            tex2.SetTextFont(61)
            tex2.SetTextSize(0.0375)
            tex2.SetLineWidth(2)
            tex2.Draw("goff")
            tex3 = ROOT.TLatex(0.23, 0.96, "Preliminary")
            tex3.SetNDC()
            tex3.SetTextFont(52)
            tex3.SetTextSize(0.0285)
            tex3.SetLineWidth(2)
            tex3.Draw("goff")

            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.pdf")
            c1.SaveAs(SR_PLOT_DIR + "/m" + out_label + "_auto.C")
            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.root")

            _file_data.Close()


def SRTotalPlot():
    m4l_cfg = m4l_plot_config()
    plot_specs = [
        ("SR_total", ""),
        ("SR_total_4e", "_4e"),
        ("SR_total_4mu", "_4mu"),
        ("SR_total_2e2mu", "_2e2mu"),
        ("SR_total_2mu2e", "_2mu2e"),
    ]

    for fr_tag in FR_SR_TAGS:
        if not _fr_sr_estimate_available(fr_tag):
            print("SRTotalPlot: skipping %s, missing estimateZX* files" % fr_tag)
            continue
        data_base, zz_base = _fr_sr_estimate_basenames(fr_tag)
        fr_suffix = _fr_sr_file_suffix(fr_tag)
        total_label = FR_SR_TOTAL_LABELS[fr_tag]
        method_label = FR_SR_METHOD_LABELS[fr_tag]
        sr_tot = _build_sr_tot_histograms(data_base, zz_base)

        for plot_label, sr_suffix in plot_specs:
            out_label = plot_label + fr_suffix
            print("SRTotalPlot: " + out_label)
            zxPlot = sr_tot[sr_suffix].Clone("zx_total_" + out_label)
            zxPlot.SetDirectory(0)
            print("Z+X total " + out_label + " := " + str(zxPlot.Integral()))

            c1 = ROOT.TCanvas("c1_" + out_label, "myPlots", 0, 67, 600, 600)
            leg = ROOT.TLegend(0.40, 0.55, 0.95, 0.85)
            leg.SetBorderSize(0)
            leg.SetLineColor(1)
            leg.SetLineStyle(1)
            leg.SetLineWidth(1)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.032)

            flavor_label = plot_label.replace("SR_total", "").strip("_")
            leg_title = "Total Z+X" if not flavor_label else "Total Z+X " + flavor_label
            if method_label:
                leg_title += " (" + method_label + ")"
            entry = leg.AddEntry("NULL", leg_title, "h")
            entry.SetLineColor(1)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(62)

            ROOT.gStyle.SetOptFit(1)
            ROOT.gStyle.SetOptStat(0)
            ROOT.gStyle.SetOptTitle(0)
            c1.SetFillColor(0)
            c1.SetBorderMode(0)
            c1.SetBorderSize(2)
            c1.SetTickx(1)
            c1.SetTicky(1)
            c1.SetLeftMargin(0.15)
            c1.SetRightMargin(0.05)
            c1.SetTopMargin(0.05)
            c1.SetBottomMargin(0.13)
            c1.SetFrameFillStyle(0)
            c1.SetFrameBorderMode(0)
            c1.SetFrameFillStyle(0)
            c1.SetFrameBorderMode(0)

            zxPlot.SetLineColor(ROOT.kRed - 7)
            zxPlot.SetLineWidth(2)
            zxPlot.SetFillStyle(0)
            zxPlot.SetMarkerStyle(20)
            zxPlot.SetMarkerColor(ROOT.kRed - 7)
            zxPlot.SetMarkerSize(0.9)
            zxPlot.GetXaxis().SetRange(1, zxPlot.GetNbinsX())
            zxPlot.GetXaxis().SetLabelFont(42)
            zxPlot.GetXaxis().SetLabelSize(0.05)
            zxPlot.GetXaxis().SetTitleSize(0.06)
            zxPlot.GetXaxis().SetTitleOffset(0.9)
            zxPlot.GetXaxis().SetTitleFont(42)
            _style_m4l_axis(zxPlot, m4l_cfg)
            zxPlot.GetYaxis().SetLabelFont(42)
            zxPlot.GetYaxis().SetLabelSize(0.05)
            zxPlot.GetYaxis().SetTitleSize(0.06)
            zxPlot.GetYaxis().SetTitleOffset(1.25)
            zxPlot.GetYaxis().SetTitleFont(42)
            zxPlot.GetZaxis().SetLabelFont(42)
            zxPlot.GetZaxis().SetLabelSize(0.035)
            zxPlot.GetZaxis().SetTitleSize(0.035)
            zxPlot.GetZaxis().SetTitleFont(42)
            zxPlot.GetYaxis().SetRangeUser(0, 1.3 * max(zxPlot.GetMaximum(), 1.0))
            zxPlot.Draw("e1 goff")

            entry = leg.AddEntry(zxPlot, total_label, "lep")
            entry.SetLineColor(ROOT.kRed - 7)
            entry.SetLineStyle(1)
            entry.SetLineWidth(2)
            entry.SetMarkerColor(ROOT.kRed - 7)
            entry.SetMarkerStyle(20)
            entry.SetMarkerSize(1.2)
            entry.SetTextFont(42)

            leg.Draw("goff")

            tex1 = ROOT.TLatex(0.95, 0.96, "109.82 fb^{-1} (13.6 TeV)")
            tex1.SetNDC()
            tex1.SetTextAlign(31)
            tex1.SetTextFont(42)
            tex1.SetTextSize(0.03)
            tex1.SetLineWidth(2)
            tex1.Draw("goff")
            tex2 = ROOT.TLatex(0.15, 0.96, "CMS")
            tex2.SetNDC()
            tex2.SetTextFont(61)
            tex2.SetTextSize(0.0375)
            tex2.SetLineWidth(2)
            tex2.Draw("goff")
            tex3 = ROOT.TLatex(0.23, 0.96, "Preliminary")
            tex3.SetNDC()
            tex3.SetTextFont(52)
            tex3.SetTextSize(0.0285)
            tex3.SetLineWidth(2)
            tex3.Draw("goff")

            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.pdf")
            c1.SaveAs(SR_PLOT_DIR + "/m" + out_label + "_auto.C")
            c1.SaveAs(SR_PLOT_DIR + "/" + out_label + "_auto.root")


def _nd_process_hist_path(var_name, process_name, jet_tag=None):
    if jet_tag is None:
        if process_name == "Data":
            return "Hist_Data_%s_Data.root" % var_name
        return "Hist_MC_%s_%s.root" % (var_name, process_name)
    if process_name == "Data":
        return "Hist_Data_%s_VBSFR_%s_Data.root" % (var_name, jet_tag)
    return "Hist_MC_%s_VBSFR_%s_%s.root" % (var_name, jet_tag, process_name)


def PlotNDProcess(Type, VarP, jet_tag=None):

    plot_dir = MET_PLOT_DIR if VarP == "mEt" else Z1L_PLOT_DIR
    out_suffix = "" if jet_tag is None else "_%s" % jet_tag
    if not os.path.isfile(_nd_process_hist_path(VarP, "Data", jet_tag)):
        print("PlotNDProcess: skipping, missing %s" % _nd_process_hist_path(VarP, "Data", jet_tag))
        return
    histName = ["el_EB","el_EE","mu_EB","mu_EE"]
    histNameLabel = ["electron barrel","electron endcap","muon barrel","muon endcap"]
    StackTitle = "Proba"

    


    ProcessNames = ["Data", "WZ", "ZZ", "TT", "DY50"]
    fileSize = 5
    Fillcolors = ["#cc0099", "#99ccff", "#996666", "#669966"]
    Linecolors = ["#990066", "#000099", "#5f3f3f", "#003300"]
    plotSize = 4
    variableN = "ptl3"

    for  p in range (plotSize): 
        
        leg_xl = 0.58
        leg_xr = 0.98
        leg_yb = 0.55
        leg_yt = 0.85

        if (VarP == "mZ1"):
            leg_xl = 0.25
            leg_xr = 0.75
            leg_yb = 0.55
            leg_yt = 0.85 



        h = ""

        if (Type == "n"):
            h ="numerator"
        else:
            h = "denominator"
        
        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)

        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)

        entry = leg.AddEntry(ROOT.nullptr, histNameLabel[p]+" "+h,"h")
        
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        entry=leg.AddEntry(ROOT.nullptr,"Data","LP")
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(20)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry(ROOT.nullptr,"Z + jets","F")
    
        ci = ROOT.TColor.GetColor("#669966")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#003300")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry(ROOT.nullptr,"t#bar{t} + jets","F")
        
        ci = ROOT.TColor.GetColor("#996666")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#5f3f3f")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)

        entry=leg.AddEntry(ROOT.nullptr,"WZ","F")
        
        ci = ROOT.TColor.GetColor("#cc0099")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#990066")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)

        entry=leg.AddEntry(ROOT.nullptr,"Z#gamma*,ZZ","F")
        
        ci = ROOT.TColor.GetColor("#99ccff")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
    
        
        A = ROOT.THStack(StackTitle,StackTitle)
        _file = []
        for i in range (fileSize):
            
            Title = "MC"
            if(i == 0):
                Title = "Data"
            
            _file.append(ROOT.TFile(_nd_process_hist_path(VarP, ProcessNames[i], jet_tag)))
            print (_file)
            if (i==0):
                dataPlot = _file[i].Get("Data_FR"+histName[p]+"_"+Type)
                dataPlot.SetFillColor(1)
                dataPlot.SetLineColor(1)
                dataPlot.SetFillStyle(0)
                dataPlot.SetMarkerStyle(20)
                dataPlot.GetXaxis().SetRange(3,40)
                dataPlot.GetXaxis().SetLabelFont(42)
                dataPlot.GetXaxis().SetLabelSize(0.05)
                dataPlot.GetXaxis().SetTitleSize(0.06)
                dataPlot.GetXaxis().SetTitleOffset(0.9)
                dataPlot.GetXaxis().SetTitleFont(42)
                dataPlot.GetYaxis().SetLabelFont(42)
                dataPlot.GetYaxis().SetLabelSize(0.05)
                dataPlot.GetYaxis().SetTitleSize(0.06)
                dataPlot.GetYaxis().SetTitleOffset(1.6)
                dataPlot.GetYaxis().SetTitleFont(42)
                dataPlot.GetZaxis().SetLabelFont(40)
                dataPlot.GetZaxis().SetLabelSize(0.035)
                dataPlot.GetZaxis().SetTitleSize(0.035)
                dataPlot.GetZaxis().SetTitleFont(42)
                dataPlot.GetXaxis().SetTitle(VarP+" [GeV]")
                
                if (VarP=="mZ1"):
                    dataPlot.GetYaxis().SetTitle("Events/3 GeV")
                if (VarP=="mZ1"):
                    dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
                if (VarP=="mEt"):
                    dataPlot.GetYaxis().SetTitle("Events/5 GeV")
                if (VarP=="mEt"):
                    dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
                if (VarP=="ptl3"):
                    dataPlot.GetXaxis().SetRangeUser(5,80)
                    dataPlot.SetMinimum(0)
                if (VarP=="mEt"):
                    dataPlot.GetXaxis().SetRangeUser(0,50)
                    dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
                
                if (VarP=="mZ1"):
                    dataPlot.GetXaxis().SetRangeUser(60,120)
                    dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
 
            else:
                histPlot = _file[i].Get("Data_FR"+histName[p]+"_"+Type)
                histPlot.SetFillColor(ROOT.TColor.GetColor(Fillcolors[i-1]))
                histPlot.SetLineColor(ROOT.TColor.GetColor(Linecolors[i-1]))
                histPlot.SetFillStyle(1001)
                
                histPlot.SetLineStyle(0)
                histPlot.SetMarkerStyle(20)
                histPlot.GetXaxis().SetLabelFont(42)
                histPlot.GetXaxis().SetLabelOffset(0.007)
                histPlot.GetXaxis().SetLabelSize(0.05)
                histPlot.GetXaxis().SetTitleSize(0.06)
                histPlot.GetXaxis().SetTitleOffset(0.9)
                histPlot.GetXaxis().SetTitleFont(42)
                histPlot.GetYaxis().SetLabelFont(42)
                histPlot.GetYaxis().SetLabelOffset(0.007)
                histPlot.GetYaxis().SetLabelSize(0.05)
                histPlot.GetYaxis().SetTitleSize(0.06)
                histPlot.GetYaxis().SetTitleOffset(1.25)
                histPlot.GetYaxis().SetTitleFont(42)
                histPlot.GetZaxis().SetLabelFont(42)
                histPlot.GetZaxis().SetLabelOffset(0.007)
                histPlot.GetZaxis().SetLabelSize(0.05)
                histPlot.GetZaxis().SetTitleSize(0.06)
                histPlot.GetZaxis().SetTitleFont(42)
                A.Add(histPlot)

        
        # setup environment & canvas
        c1 = ROOT.TCanvas("c1_"+VarP+"_"+histName[p]+"_"+Type+"_process"+out_suffix, "myPlots",0,75,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.18)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        _overall_max = max(dataPlot.GetMaximum(), A.GetMaximum())
        dataPlot.SetMaximum(1.5*_overall_max)
        if (p>3):
            dataPlot.SetMaximum(2*_overall_max)

        dataPlot.Draw("e1 goff")

        if (VarP == "mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
        elif (VarP == "mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)

        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        leg.Draw("goff")
        
        tex = ROOT.TLatex(0.95,0.96,"109.82 fb^{-1} (13.6 TeV)")
        tex.SetNDC()
        tex.SetTextAlign(31)
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex1 = ROOT.TLatex(0.15,0.96,"CMS")
        tex1.SetNDC()
        tex1.SetTextFont(61)
        tex1.SetTextSize(0.0375)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex2.SetNDC()
        tex2.SetTextFont(52)
        tex2.SetTextSize(0.0285)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        
        
        c1.SaveAs(plot_dir+"/"+VarP+histName[p]+"_"+Type+"_process_auto"+out_suffix+".pdf")

def _nd_id_hist_path(var_name, process_name, jet_tag=None):
    if jet_tag is None:
        return "Hist_ID_%s_%s.root" % (var_name, process_name)
    return "Hist_ID_%s_VBSFR_%s_%s.root" % (var_name, jet_tag, process_name)


def PlotNDMatched(Type, VarP, jet_tag=None):

    plot_dir = MET_PLOT_DIR if VarP == "mEt" else Z1L_PLOT_DIR
    out_suffix = "" if jet_tag is None else "_%s" % jet_tag
    data_path = _nd_process_hist_path(VarP, "Data", jet_tag)
    if not os.path.isfile(data_path):
        print("PlotNDMatched: skipping, missing %s" % data_path)
        return
    if not os.path.isfile(_nd_id_hist_path(VarP, "ZZ", jet_tag)):
        print("PlotNDMatched: skipping, missing %s" % _nd_id_hist_path(VarP, "ZZ", jet_tag))
        return
    
    histNameLabel = ["electron endcap","electron barrel","muon endcap","muon barrel"]
    
    histName = ["EE","EB","ME","MB"]
    dataName = ["el_EE","el_EB","mu_EE","mu_EB"]
    StackTitle = "Proba"

    # Unused in this routine; kept only to mirror other functions’ legacy bookkeeping.
    fileSize = 5
    # fileSize = 2  # when only Data+WZ FR files exist
    plotSize = 4
    variableN = "ptl3"

    for p in range (plotSize):
        
        h = ""
        if (Type == "n"):
            h ="numerator"
        else:
            h = "denominator"
        

        leg_xl = 0.58
        leg_xr = 0.98
        leg_yb = 0.55
        leg_yt = 0.85

        if (VarP == "mZ1"):
            leg_xl = 0.25
            leg_xr = 0.75
            leg_yb = 0.55
            leg_yt = 0.85 

        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        entry = leg.AddEntry(
            ROOT.nullptr, histNameLabel[p] + " " + h, "h"
        )
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        
        _file1 = ROOT.TFile(_nd_process_hist_path(VarP, "Data", jet_tag), "READ")

        print ("Data_FR"+dataName[p]+"_"+Type)
        
        dataPlot = _file1.Get("Data_FR"+dataName[p]+"_"+Type)
        
        dataPlot.SetFillColor(1)
        dataPlot.SetLineColor(1)
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetRange(3,40)
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.6)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(40)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        dataPlot.GetXaxis().SetTitle(VarP+" [GeV]")
        
        if (VarP=="mZ1"):
            dataPlot.GetYaxis().SetTitle("Events/3 GeV")
        if (VarP=="mZ1"):
            dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
        if (VarP=="mEt"):
            dataPlot.GetYaxis().SetTitle("Events/5 GeV")
        if (VarP=="mEt"):
            dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
        if (VarP=="ptl3"):
            dataPlot.GetXaxis().SetRangeUser(5,80)
            dataPlot.SetMinimum(0)
        if (VarP=="mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)
            dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
        
        if (VarP=="mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
            dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())


        if (VarP=="mZ1"):
            dataPlot.GetYaxis().SetTitle("Events/3 GeV")
            dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
        if (VarP=="mEt"):
            dataPlot.GetYaxis().SetTitle("Events/5 GeV")
            dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
        

        Title = "Data"
        leg.AddEntry(dataPlot,Title,"LP")
        
        _fileZZ = ROOT.TFile(_nd_id_hist_path(VarP, "ZZ", jet_tag), "READ")
        _fileWZ = ROOT.TFile(_nd_id_hist_path(VarP, "WZ", jet_tag), "READ")
        _fileTT = ROOT.TFile(_nd_id_hist_path(VarP, "TT", jet_tag), "READ")
        _fileDY = ROOT.TFile(_nd_id_hist_path(VarP, "DY50", jet_tag), "READ")
        A = ROOT.THStack(StackTitle,StackTitle)
        
        histPlotC = _fileZZ.Get("Hist_conv"+histName[p]+"_"+Type).Clone("Hist_conv"+histName[p]+"_"+Type+"_sum")
        histPlotC.SetDirectory(0)
        histPlotC.Add(_fileWZ.Get("Hist_conv"+histName[p]+"_"+Type))
        histPlotC.Add(_fileTT.Get("Hist_conv"+histName[p]+"_"+Type))
        histPlotC.Add(_fileDY.Get("Hist_conv"+histName[p]+"_"+Type))
        histPlotC.SetFillColor(30)
        histPlotC.SetLineColor(1)
        histPlotC.SetFillStyle(1001)
        
        histPlotC.SetLineStyle(0)
        histPlotC.SetMarkerStyle(20)
        histPlotC.GetXaxis().SetLabelFont(42)
        histPlotC.GetXaxis().SetLabelOffset(0.007)
        histPlotC.GetXaxis().SetLabelSize(0.05)
        histPlotC.GetXaxis().SetTitleSize(0.06)
        histPlotC.GetXaxis().SetTitleOffset(0.9)
        histPlotC.GetXaxis().SetTitleFont(42)
        histPlotC.GetYaxis().SetLabelFont(42)
        histPlotC.GetYaxis().SetLabelOffset(0.007)
        histPlotC.GetYaxis().SetLabelSize(0.05)
        histPlotC.GetYaxis().SetTitleSize(0.06)
        histPlotC.GetYaxis().SetTitleOffset(1.25)
        histPlotC.GetYaxis().SetTitleFont(42)
        histPlotC.GetZaxis().SetLabelFont(42)
        histPlotC.GetZaxis().SetLabelOffset(0.007)
        histPlotC.GetZaxis().SetLabelSize(0.05)
        histPlotC.GetZaxis().SetTitleSize(0.06)
        histPlotC.GetZaxis().SetTitleFont(42)
        
        
        leg.AddEntry(histPlotC,"Conversions","F")
        A.Add(histPlotC)
        
        histPlotP = _fileZZ.Get("Hist_prompt"+histName[p]+"_"+Type).Clone("Hist_prompt"+histName[p]+"_"+Type+"_sum")
        histPlotP.SetDirectory(0)
        histPlotP.Add(_fileWZ.Get("Hist_prompt"+histName[p]+"_"+Type))
        histPlotP.Add(_fileTT.Get("Hist_prompt"+histName[p]+"_"+Type))
        histPlotP.Add(_fileDY.Get("Hist_prompt"+histName[p]+"_"+Type))
        histPlotP.SetFillColor(38)
        histPlotP.SetLineColor(1)
        histPlotP.SetFillStyle(1001)
        histPlotP.SetLineStyle(0)
        histPlotP.SetMarkerStyle(20)
        histPlotP.GetXaxis().SetLabelFont(42)
        histPlotP.GetXaxis().SetLabelOffset(0.007)
        histPlotP.GetXaxis().SetLabelSize(0.05)
        histPlotP.GetXaxis().SetTitleSize(0.06)
        histPlotP.GetXaxis().SetTitleOffset(0.9)
        histPlotP.GetXaxis().SetTitleFont(42)
        histPlotP.GetYaxis().SetLabelFont(42)
        histPlotP.GetYaxis().SetLabelOffset(0.007)
        histPlotP.GetYaxis().SetLabelSize(0.05)
        histPlotP.GetYaxis().SetTitleSize(0.06)
        histPlotP.GetYaxis().SetTitleOffset(1.25)
        histPlotP.GetYaxis().SetTitleFont(42)
        histPlotP.GetZaxis().SetLabelFont(42)
        histPlotP.GetZaxis().SetLabelOffset(0.007)
        histPlotP.GetZaxis().SetLabelSize(0.05)
        histPlotP.GetZaxis().SetTitleSize(0.06)
        histPlotP.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotP,"Prompt","F")
        A.Add(histPlotP)
        
        histPlotF = _fileZZ.Get("Hist_fakes"+histName[p]+"_"+Type).Clone("Hist_fakes"+histName[p]+"_"+Type+"_sum")
        histPlotF.SetDirectory(0)
        histPlotF.Add(_fileWZ.Get("Hist_fakes"+histName[p]+"_"+Type))
        histPlotF.Add(_fileTT.Get("Hist_fakes"+histName[p]+"_"+Type))
        histPlotF.Add(_fileDY.Get("Hist_fakes"+histName[p]+"_"+Type))
        histPlotF.SetFillColor(46)
        histPlotF.SetLineColor(1)
        histPlotF.SetFillStyle(1001)
        histPlotF.SetLineStyle(0)
        histPlotF.SetMarkerStyle(20)
        histPlotF.GetXaxis().SetLabelFont(42)
        histPlotF.GetXaxis().SetLabelOffset(0.007)
        histPlotF.GetXaxis().SetLabelSize(0.05)
        histPlotF.GetXaxis().SetTitleSize(0.06)
        histPlotF.GetXaxis().SetTitleOffset(0.9)
        histPlotF.GetXaxis().SetTitleFont(42)
        histPlotF.GetYaxis().SetLabelFont(42)
        histPlotF.GetYaxis().SetLabelOffset(0.007)
        histPlotF.GetYaxis().SetLabelSize(0.05)
        histPlotF.GetYaxis().SetTitleSize(0.06)
        histPlotF.GetYaxis().SetTitleOffset(1.25)
        histPlotF.GetYaxis().SetTitleFont(42)
        histPlotF.GetZaxis().SetLabelFont(42)
        histPlotF.GetZaxis().SetLabelOffset(0.007)
        histPlotF.GetZaxis().SetLabelSize(0.05)
        histPlotF.GetZaxis().SetTitleSize(0.06)
        histPlotF.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotF,"LF Fakes","F")
        A.Add(histPlotF)
        
        histPlotB = _fileZZ.Get("Hist_BDfakes"+histName[p]+"_"+Type).Clone("Hist_BDfakes"+histName[p]+"_"+Type+"_sum")
        histPlotB.SetDirectory(0)
        histPlotB.Add(_fileWZ.Get("Hist_BDfakes"+histName[p]+"_"+Type))
        histPlotB.Add(_fileTT.Get("Hist_BDfakes"+histName[p]+"_"+Type))
        histPlotB.Add(_fileDY.Get("Hist_BDfakes"+histName[p]+"_"+Type))
        histPlotB.SetFillColor(49)
        histPlotB.SetLineColor(1)
        histPlotB.SetFillStyle(1001)
        histPlotB.SetLineStyle(0)
        histPlotB.SetMarkerStyle(20)
        histPlotB.GetXaxis().SetLabelFont(42)
        histPlotB.GetXaxis().SetLabelOffset(0.007)
        histPlotB.GetXaxis().SetLabelSize(0.05)
        histPlotB.GetXaxis().SetTitleSize(0.06)
        histPlotB.GetXaxis().SetTitleOffset(0.9)
        histPlotB.GetXaxis().SetTitleFont(42)
        histPlotB.GetYaxis().SetLabelFont(42)
        histPlotB.GetYaxis().SetLabelOffset(0.007)
        histPlotB.GetYaxis().SetLabelSize(0.05)
        histPlotB.GetYaxis().SetTitleSize(0.06)
        histPlotB.GetYaxis().SetTitleOffset(1.25)
        histPlotB.GetYaxis().SetTitleFont(42)
        histPlotB.GetZaxis().SetLabelFont(42)
        histPlotB.GetZaxis().SetLabelOffset(0.007)
        histPlotB.GetZaxis().SetLabelSize(0.05)
        histPlotB.GetZaxis().SetTitleSize(0.06)
        histPlotB.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotB,"HF fakes","F")
        A.Add(histPlotB)
        
        # setup environment & canvas
        c1 = ROOT.TCanvas("c1_"+VarP+"_"+histName[p]+"_"+Type+"_match"+out_suffix, "myPlots",0,75,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.18)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        dataPlot.Draw("e1")

        _overall_max = max(dataPlot.GetMaximum(), A.GetMaximum())
        dataPlot.SetMaximum(1.5*_overall_max)
        if (p>3):
            dataPlot.SetMaximum(2*_overall_max)

        if (VarP == "mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
        elif (VarP == "mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)

        dataPlot.Draw("e1 goff")

        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        leg.Draw("goff")
        
        
        tex = ROOT.TLatex(0.95,0.96,"109.82 fb^{-1} (13.6 TeV)")
        tex.SetNDC()
        tex.SetTextAlign(31)
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex1 = ROOT.TLatex(0.15,0.96,"CMS")
        tex1.SetNDC()
        tex1.SetTextFont(61)
        tex1.SetTextSize(0.0375)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex2.SetNDC()
        tex2.SetTextFont(52)
        tex2.SetTextSize(0.0285)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        
        
        c1.SaveAs(plot_dir+"/"+VarP+"_"+Type+"_"+histName[p]+"_match_auto"+out_suffix+".pdf")
    
        _fileZZ.Close()
        _fileWZ.Close()
        _fileTT.Close()
        _fileDY.Close()


def PlotFRvsPtl3():
    # (data file, WZ-removed file, output-name suffix, legend tag). The ge1j/ge2j
    # variants only exist when VBS=True; they are skipped gracefully otherwise.
    fr_variants = [
        ("Hist_Data_ptl3_Data.root", "Hist_Data_ptl3_WZremoved.root", "", ""),
        (
            "Hist_Data_ptl3_VBSFR_ge1j_Data.root",
            "Hist_Data_ptl3_VBSFR_ge1j_WZremoved.root",
            "_ge1j",
            " (#geq1 jet)",
        ),
        (
            "Hist_Data_ptl3_VBSFR_ge2j_Data.root",
            "Hist_Data_ptl3_VBSFR_ge2j_WZremoved.root",
            "_ge2j",
            " (#geq2 jets)",
        ),
    ]
    for data_path, wz_path, out_suffix, tag in fr_variants:
        _plot_fr_vs_ptl3_variant(data_path, wz_path, out_suffix, tag)


def _plot_fr_vs_ptl3_variant(data_path, wz_path, out_suffix, tag):
    if not os.path.isfile(data_path) or not os.path.isfile(wz_path):
        print("PlotFRvsPtl3: skipping variant, missing %s or %s" % (data_path, wz_path))
        return

    data_file = ROOT.TFile(data_path, "READ")
    wz_removed_file = ROOT.TFile(wz_path, "READ")

    if not data_file or data_file.IsZombie():
        print("Could not open %s" % data_path)
        return
    if not wz_removed_file or wz_removed_file.IsZombie():
        print("Could not open %s" % wz_path)
        data_file.Close()
        return

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # (label, EB hist, EE hist, y_max, output basename)
    panel_cfg = [
        ("electrons", "Data_FRel_EB", "Data_FRel_EE", 0.12, "FR_vs_ptl3_electrons" + out_suffix),
        ("muons",     "Data_FRmu_EB", "Data_FRmu_EE", 0.5, "FR_vs_ptl3_muons" + out_suffix),
    ]

    for lep_label, h_eb_name, h_ee_name, y_max, out_base in panel_cfg:
        h_eb_data = data_file.Get(h_eb_name)
        h_ee_data = data_file.Get(h_ee_name)
        h_eb_wzr = wz_removed_file.Get(h_eb_name)
        h_ee_wzr = wz_removed_file.Get(h_ee_name)

        if not (h_eb_data and h_ee_data and h_eb_wzr and h_ee_wzr):
            print("Missing FR histogram(s):", h_eb_name, h_ee_name)
            continue

        # Detach from file so they survive after Close().
        for h in (h_eb_data, h_ee_data, h_eb_wzr, h_ee_wzr):
            h.SetDirectory(0)

        eb_light = ROOT.TColor.GetColor("#E2A5A5")  # electrons/muons EB Data (light red)
        ee_light = ROOT.TColor.GetColor("#A5C9E2")  # electrons/muons EE Data (light blue)
        eb_dark  = ROOT.TColor.GetColor("#7A1F1F")  # EB WZ removed (dark red)
        ee_dark  = ROOT.TColor.GetColor("#1F2F7A")  # EE WZ removed (dark blue)

        styles = [
            (h_eb_data, eb_light, 20, 0.6),
            (h_ee_data, ee_light, 20, 0.6),
            (h_eb_wzr,  eb_dark,  20, 0.6),
            (h_ee_wzr,  ee_dark,  20, 0.6),
        ]
        for h, col, mstyle, msize in styles:
            h.SetLineColor(col)
            h.SetMarkerColor(col)
            h.SetLineWidth(2)
            h.SetMarkerStyle(mstyle)
            h.SetMarkerSize(msize)

        # Use first hist for axis cosmetics.
        h_eb_data.GetXaxis().SetTitle("ptl3 [GeV]")
        h_eb_data.GetYaxis().SetTitle("Fake Rate")
        h_eb_data.GetXaxis().SetRangeUser(0, 80)
        h_eb_data.GetYaxis().SetRangeUser(0.0, y_max)
        h_eb_data.GetXaxis().SetTitleSize(0.05)
        h_eb_data.GetYaxis().SetTitleSize(0.05)
        h_eb_data.GetXaxis().SetLabelSize(0.04)
        h_eb_data.GetYaxis().SetLabelSize(0.04)
        h_eb_data.GetXaxis().SetTitleOffset(1.05)
        h_eb_data.GetYaxis().SetTitleOffset(1.10)

        c = ROOT.TCanvas("c_FR_" + lep_label + out_suffix, "FR vs ptl3 " + lep_label + tag, 700, 600)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.04)
        c.SetTopMargin(0.07)
        c.SetBottomMargin(0.13)
        c.SetTickx(1)
        c.SetTicky(1)

        # "E1" gives points + vertical err + horizontal err = bin width.
        h_eb_data.Draw("E1")
        h_ee_data.Draw("E1 SAME")
        h_eb_wzr.Draw("E1 SAME")
        h_ee_wzr.Draw("E1 SAME")

        leg = ROOT.TLegend(0.42, 0.68, 0.93, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.032)
        leg.SetHeader("#bf{pp #rightarrow Z + l_{loose}}" + tag)
        leg.AddEntry(h_eb_data, lep_label + ", EB",            "LP")
        leg.AddEntry(h_ee_data, lep_label + ", EE",            "LP")
        leg.AddEntry(h_eb_wzr,  lep_label + ", EB WZ removed", "LP")
        leg.AddEntry(h_ee_wzr,  lep_label + ", EE WZ removed", "LP")
        leg.Draw()

        tex_lumi = ROOT.TLatex(0.96, 0.94, "109.82 fb^{-1} (13.6 TeV)")
        tex_lumi.SetNDC(); tex_lumi.SetTextAlign(31); tex_lumi.SetTextFont(42); tex_lumi.SetTextSize(0.035)
        tex_lumi.Draw()
        tex_cms = ROOT.TLatex(0.14, 0.94, "CMS")
        tex_cms.SetNDC(); tex_cms.SetTextFont(61); tex_cms.SetTextSize(0.045)
        tex_cms.Draw()
        tex_prelim = ROOT.TLatex(0.23, 0.94, "Preliminary")
        tex_prelim.SetNDC(); tex_prelim.SetTextFont(52); tex_prelim.SetTextSize(0.034)
        tex_prelim.Draw()

        c.SaveAs(os.path.join(FR_PTL3_PLOT_DIR, out_base + ".pdf"))

    data_file.Close()
    wz_removed_file.Close()


def PlotFRvsMZ1():
    data_file = ROOT.TFile("Hist_Data_mZ1_Data.root", "READ")
    if not data_file or data_file.IsZombie():
        print("Could not open Hist_Data_mZ1_Data.root")
        return

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # (label, EB hist, EE hist, y_max, output basename)
    panel_cfg = [
        ("electrons", "Data_FRel_EB", "Data_FRel_EE", 0.6, "FR_vs_mZ1_electrons"),
        ("muons",     "Data_FRmu_EB", "Data_FRmu_EE", 0.5, "FR_vs_mZ1_muons"),
    ]

    for lep_label, h_eb_name, h_ee_name, y_max, out_base in panel_cfg:
        h_eb = data_file.Get(h_eb_name)
        h_ee = data_file.Get(h_ee_name)
        if not (h_eb and h_ee):
            print("Missing FR histogram(s):", h_eb_name, h_ee_name)
            continue
        for h in (h_eb, h_ee):
            h.SetDirectory(0)

        col_eb = ROOT.TColor.GetColor("#C53030")  # red
        col_ee = ROOT.TColor.GetColor("#3858B0")  # blue
        for h, col in ((h_eb, col_eb), (h_ee, col_ee)):
            h.SetLineColor(col)
            h.SetMarkerColor(col)
            h.SetLineWidth(2)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.6)

        h_eb.GetXaxis().SetTitle("m_{Z_{1}} [GeV]")
        h_eb.GetYaxis().SetTitle("Fake Rate")
        h_eb.GetXaxis().SetRangeUser(60, 120)
        h_eb.GetYaxis().SetRangeUser(0.0, y_max)
        h_eb.GetXaxis().SetTitleSize(0.05)
        h_eb.GetYaxis().SetTitleSize(0.05)
        h_eb.GetXaxis().SetLabelSize(0.04)
        h_eb.GetYaxis().SetLabelSize(0.04)
        h_eb.GetXaxis().SetTitleOffset(1.05)
        h_eb.GetYaxis().SetTitleOffset(1.10)

        c = ROOT.TCanvas("c_FRmZ1_" + lep_label, "FR vs mZ1 " + lep_label, 700, 600)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.04)
        c.SetTopMargin(0.07)
        c.SetBottomMargin(0.13)
        c.SetTickx(1)
        c.SetTicky(1)

        h_eb.Draw("E1")
        h_ee.Draw("E1 SAME")

        # Z-mass window used in the ptl3 FR cut: |mZ1 - 91.188| < 7 GeV.
        line_lo = ROOT.TLine(91.188 - 7.0, 0.0, 91.188 - 7.0, y_max)
        line_hi = ROOT.TLine(91.188 + 7.0, 0.0, 91.188 + 7.0, y_max)
        for ln in (line_lo, line_hi):
            ln.SetLineColor(ROOT.TColor.GetColor("#6A3FA0"))
            ln.SetLineWidth(2)
            ln.Draw()

        leg = ROOT.TLegend(0.42, 0.74, 0.93, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.032)
        leg.SetHeader("#bf{pp #rightarrow Z + l_{loose}}")
        leg.AddEntry(h_eb, lep_label + ", EB", "LP")
        leg.AddEntry(h_ee, lep_label + ", EE", "LP")
        leg.Draw()

        tex_lumi = ROOT.TLatex(0.96, 0.94, "109.82 fb^{-1} (13.6 TeV)")
        tex_lumi.SetNDC(); tex_lumi.SetTextAlign(31); tex_lumi.SetTextFont(42); tex_lumi.SetTextSize(0.035)
        tex_lumi.Draw()
        tex_cms = ROOT.TLatex(0.14, 0.94, "CMS")
        tex_cms.SetNDC(); tex_cms.SetTextFont(61); tex_cms.SetTextSize(0.045)
        tex_cms.Draw()
        tex_prelim = ROOT.TLatex(0.23, 0.94, "Preliminary")
        tex_prelim.SetNDC(); tex_prelim.SetTextFont(52); tex_prelim.SetTextSize(0.034)
        tex_prelim.Draw()

        c.SaveAs(os.path.join(FR_MZ1_PLOT_DIR, out_base + ".pdf"))
        c.SaveAs(os.path.join(FR_MZ1_PLOT_DIR, out_base + ".root"))

    data_file.Close()


if __name__ == "__main__":
    for output_dir in (
        CR_PLOT_DIR,
        SR_PLOT_DIR,
        Z1L_PLOT_DIR,
        MET_PLOT_DIR,
        FR_PTL3_PLOT_DIR,
        FR_MZ1_PLOT_DIR,
    ):
        os.makedirs(output_dir, exist_ok=True)


    os.chdir(STAGE1_DIR)

    CRPlot()
    CRMatchPlot()
    SRInputPlot()
    SRTotalPlot()
    SRPlot()
    SRMatchPlot()
    PlotFRvsPtl3()
    PlotFRvsMZ1()

    for var in ("ptl3", "mZ1", "mEt"):
        PlotNDProcess("n", var)
        PlotNDProcess("d", var)
        PlotNDMatched("n", var)
        PlotNDMatched("d", var)
        if var == "ptl3":
            for jet_tag in ("ge1j", "ge2j"):
                if not os.path.isfile(_nd_process_hist_path(var, "Data", jet_tag)):
                    continue
                PlotNDProcess("n", var, jet_tag=jet_tag)
                PlotNDProcess("d", var, jet_tag=jet_tag)
                PlotNDMatched("n", var, jet_tag=jet_tag)
                PlotNDMatched("d", var, jet_tag=jet_tag)