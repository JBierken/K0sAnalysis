####################################################################
# function for plotting a histogram along with its fitted function #
####################################################################

import ROOT
import sys
import os
sys.path.append('../')
import tools.plottools as pt

def plot_fit(hist,figname,fitfunc=None,backfit=None,
             label=None,paramdict=None,xaxtitle=None,yaxtitle=None,
             extrainfo='',lumi=0): 
        # args: - hist is the histogram
        #        - figname is the name of the output figure
        #        - fitfunc and backfit are two fitted functions 
        #          (typically total and background only, both can be None)
        #         - label is the legend entry for the histogram
        #        - paramdict is a dict containing names and values of fitted parameters of fitfunc
        #         - xaxtitle and yaxtitle are axis titles
        #         - extrainfo is a string with extra info to show, use '<<' to split lines
        #        - lumi is the luminosity to display in the header (in pb-1!)
        ROOT.gROOT.SetBatch(ROOT.kTRUE);
        pt.setTDRstyle()

        # initialization of canvas
        c1 = ROOT.TCanvas("c1","c1")
        c1.SetCanvasSize(800,800)
        titlefont = 4; titlesize = 30
        labelfont = 4; labelsize = 30
        axtitlefont = 4; axtitlesize = 35
        infofont = 4; infosize = 25
        legendfont = 4; legendsize = 25
        c1.SetBottomMargin(0.1)
        c1.SetTopMargin(0.05)
        c1.SetLeftMargin(0.15)
        if label is None: label = 'histogram'

        # legend
        leg = ROOT.TLegend(0.65,0.8,0.88,0.9)
        leg.SetTextFont(10*legendfont+3)
        leg.SetTextSize(legendsize)
        leg.SetBorderSize(0)

        # draw histogram
        hist.SetStats(False)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetLineWidth(2)
        hist.Sumw2()
        leg.AddEntry(hist,label,"l")
        hist.Draw()

        # X-axis layout
        xax = hist.GetXaxis()
        xax.SetNdivisions(10,4,0,ROOT.kTRUE)
        xax.SetLabelFont(10*labelfont+3)
        xax.SetLabelSize(labelsize)
        if xaxtitle is not None:
            xax.SetTitle(xaxtitle)
            xax.SetTitleFont(10*axtitlefont+3)
            xax.SetTitleSize(axtitlesize)

        # Y-axis layout
        hist.SetMaximum(hist.GetMaximum()*1.2)
        hist.SetMinimum(0.)
        yax = hist.GetYaxis()
        yax.SetMaxDigits(3)
        yax.SetNdivisions(8,4,0,ROOT.kTRUE)
        yax.SetLabelFont(10*labelfont+3)
        yax.SetLabelSize(labelsize)
        if yaxtitle is not None:
            yax.SetTitle(yaxtitle)
            yax.SetTitleFont(10*axtitlefont+3)
            yax.SetTitleSize(axtitlesize)
            yax.SetTitleOffset(1.5)
        ROOT.gPad.SetTicks(1,1)
        hist.Draw("HIST E")
        tinfo = ROOT.TLatex()

        # draw fitted function
        if backfit is not None:
            backfit.SetLineColor(ROOT.kGreen+3)
            backfit.SetLineWidth(3)
            backfit.SetLineStyle(ROOT.kDashed)
            backfit.Draw("SAME")
            leg.AddEntry(backfit,'background fit','l')
            leg.Draw()

        # draw fitted function
        if fitfunc is not None:
            fitfunc.SetLineColor(ROOT.kRed)
            fitfunc.SetLineWidth(3)
            fitfunc.Draw("SAME")
            leg.AddEntry(fitfunc,'total fit','l')
            leg.Draw()

            # display additional info
            tinfo.SetTextFont(10*infofont+3)
            tinfo.SetTextSize(infosize)
            if fitfunc.GetNDF()==0: normchi2 = 0
            else: normchi2 = fitfunc.GetChisquare()/fitfunc.GetNDF()
            for i,key in enumerate(paramdict):
                info = key+' : {0:.4E}'.format(paramdict[key])
                tinfo.DrawLatexNDC(0.65,0.75-i*0.035,info)
            info = r"#frac{#chi^{2}}{ndof}"+' of fit: {0:.2E}'.format(normchi2)
            tinfo.DrawLatexNDC(0.65,0.75-(len(paramdict)+1)*0.035,info)
        
        # display extra info
        tinfo.SetTextFont(10*infofont+3)
        tinfo.SetTextSize(infosize)
        extrainfo = extrainfo.replace('HACK_KS','PDG K^{0}_{S} mass:<<    497.614 #pm 0.024 MeV')
        extrainfo = extrainfo.replace('HACK_LA','PDG #Lambda^{0} mass:<<    1115.683 #pm 0.006 MeV')
        extrainfo = extrainfo.replace('HACK_JP','PDG J/#Psi mass:<<    3096.900 #pm 0.006 MeV')
        infolist = extrainfo.split('<<')
        for i,info in enumerate(infolist):
                tinfo.DrawLatexNDC(0.2,0.8-i*0.035,info)

        # title
        lumitext = '' if lumi==0 else '{0:.1f} '.format(float(lumi)/1000) + 'fb^{-1} (13 TeV)'
        pt.drawLumi(c1,extratext="Preliminary", lumitext=lumitext)

        c1.Update()
        c1.SaveAs(figname)
        c1.Close()
