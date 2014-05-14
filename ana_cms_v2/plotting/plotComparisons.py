#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="outtre.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=22)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.5)

(options, args) = parser.parse_args()

############################################################
def makeComparisonPlots(hists, names, canvasname, linecolors, linestyles, isLog=False, norm=False, directory='plots'):

    max = -999.;
    for typ in names:
        if norm:
            hists[typ].Scale(1./hists[typ].Integral());
        if max < hists[typ].GetMaximum():
            max = hists[typ].GetMaximum();
        hists[typ].SetLineWidth(2)
    
    # legend
    leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    # text
    latex = ROOT.TLatex(0.20,0.89,("Anti-kT (R=%.0f)"%(options.radius)))
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex2 = ROOT.TLatex(0.20,0.84,("n_{PU} = "+str(options.nPU)))
    latex2.SetNDC()
    latex2.SetTextSize(0.03)
    latex3 = ROOT.TLatex(0.20,0.79,("raw p_{T} > 25 GeV "))
    #latex3 = ROOT.TLatex(0.20,0.79,("raw p_{T} > %.0f GeV "%(options.minpt)))
    #if (options.useCorrPt == 1):
    #    latex3 = ROOT.TLatex(0.20,0.79,("p_{T} > %.0f GeV "%(options.minpt)))
    latex3.SetNDC()
    latex3.SetTextSize(0.03)

    # style and legend
    for typ in names:
        hists[typ].SetLineColor(linecolors[typ])
        hists[typ].SetLineStyle(linestyles[typ])
        hists[typ].SetLineWidth(2)
        leg.AddEntry(hists[typ], typ, "l");

    can = ROOT.TCanvas("c_"+canvasname,"c_"+canvasname,700,700);

                   
    for typ in names:
        if typ == names[0]:
            hists[typ].SetMaximum( 1.3*max );
            hists[typ].SetMinimum( 0 );    
            hists[typ].Draw("hist");
        else:    
            hists[typ].Draw("sames");
    leg.Draw();
    latex.Draw();
    latex2.Draw();
    latex3.Draw();
    if isLog:
        ROOT.gPad.SetLogy();
        hists[names[0]].SetMinimum( 1 );
        
    #can.SaveAs(directory+"/"+canvasname+".eps");
    can.SaveAs(directory+"/"+canvasname+".png");
    can.SaveAs(directory+"/"+canvasname+".pdf");

    #for hist in hists:

    #tmax    hist.Scale(1./hist.Integral());
    #raw_input('ok?')


def makePileupPlots(h,hpu,hgood,typ,name,directory):
    
    maxy = -1
    for key,hist in h.iteritems():
        tmp = hist.GetMaximum()
        #print key, " max = ", tmp
        if (tmp > maxy):
            maxy = tmp

    #col = h[typ].GetLineColor()
    col = ROOT.kGray+2
    hpu[typ].SetLineColor(col)
    hpu[typ].SetLineStyle(2)
    hgood[typ].SetLineColor(col)
    hgood[typ].SetLineStyle(1)

    # legend                                                                                                                                                                                                                       
    leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hgood[typ],"Real jets","L")
    leg.AddEntry(hpu[typ],"Pile-up","L")
    leg.AddEntry(h[typ],"All","L")

    # text                                                                                                  
    latex4 = ROOT.TLatex(0.20,0.89,typ);
    latex4.SetNDC()
    latex4.SetTextSize(0.03)
    latex = ROOT.TLatex(0.20,0.84,("Anti-kT (R=0.5)"));
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex2 = ROOT.TLatex(0.20,0.79,("n_{PU} = "+str(options.nPU)));
    latex2.SetNDC()
    latex2.SetTextSize(0.03)
    latex3 = ROOT.TLatex(0.20,0.74,"raw p_{T} > 25 GeV ");
    #latex3 = ROOT.TLatex(0.20,0.74,("raw p_{T} > %.0f GeV "%(options.minpt)));
    #if (options.useCorrPt == 1):
    #    latex3 = ROOT.TLatex(0.20,0.74,("p_{T} > %.0f GeV "%(options.minpt)));
    latex3.SetNDC()
    latex3.SetTextSize(0.03)

    can = ROOT.TCanvas("c_"+name+"_"+typ,"c_"+name+"_"+typ,700,700);

    #print "  max = ", maxy 
    h[typ].GetYaxis().SetRangeUser(0,maxy*1.05)
    h[typ].Draw()
    hpu[typ].Draw("same")
    hgood[typ].Draw("same")
    leg.Draw()
    latex.Draw();
    latex2.Draw();
    latex3.Draw();
    latex4.Draw();

    can.SaveAs(directory+"/"+name+"_"+typ+".png");
    can.SaveAs(directory+"/"+name+"_"+typ+".pdf");
    #raw_input('ok?') 

if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        sys.exit()

    f = ROOT.TFile.Open(filename);
    
    #types = ['GEN','PUPPI','PFlow','PFlowCHS','PFlow-CMSSW'];
    types = ['GEN','PUPPI','PFlow','PFlowCHS'];

    folder = {} 
    folder['GEN']='gen'
    folder['PUPPI']='puppi'
    folder['PFlow']='pf'
    folder['PFlowCHS']='pfchs'


    linecolors = {}
    linecolors['GEN'] = ROOT.kBlack
    linecolors['PUPPI'] = ROOT.kGreen+1
    linecolors['PFlow'] = ROOT.kBlue
    linecolors['PFlowCHS'] = ROOT.kMagenta
    #linecolors['PFlow-CMSSW'] = ROOT.kOrange

    linestyles = {}
    linestyles['GEN'] = 1
    linestyles['PUPPI'] = 1
    linestyles['PFlow'] = 1
    linestyles['PFlowCHS'] = 1
    #linestyles['PFlow-CMSSW'] = 1
    
    nTypes = len(types);
    
    # -- get histograms
    
    hnjets = {}

    hrawpt          = {}
    hrawpt_pu       = {}
    hrawpt_good     = {}
    hrawpt_response = {}
    
    hpt             = {}
    hpt_pu          = {}
    hpt_good        = {}
    hpt_response    = {}

    hpt_leadjet          = {}
    hpt_pu_leadjet       = {}
    hpt_good_leadjet     = {}
    hpt_response_leadjet = {}

    hrawpt_leadjet          = {}
    hrawpt_pu_leadjet       = {}
    hrawpt_good_leadjet     = {}
    hrawpt_response_leadjet = {}

    heta      = {}
    heta_pu   = {}
    heta_good = {}

    heta_leadjet      = {}
    heta_pu_leadjet   = {}
    heta_good_leadjet = {}

    hmass           = {}
    hmass_response  = {}
    hmass_leadjet            = {}
    hmass_response_leadjet   = {}
    
    hnparticles  = {}

    for typ in types:
        hnjets[typ]=f.Get(folder[typ]+'/hnjets')

        hrawpt[typ]=f.Get(folder[typ]+'/hrawpt')
        hrawpt_pu[typ]=f.Get(folder[typ]+'/hrawpt_pu')
        hrawpt_good[typ]=f.Get(folder[typ]+'/hrawpt_good')
        hrawpt_response[typ]=f.Get(folder[typ]+'/hrawpt_response')

        hpt[typ]=f.Get(folder[typ]+'/hpt')
        hpt_pu[typ]=f.Get(folder[typ]+'/hpt_pu')
        hpt_good[typ]=f.Get(folder[typ]+'/hpt_good')
        hpt_response[typ]=f.Get(folder[typ]+'/hpt_response')

        hrawpt_leadjet[typ]=f.Get(folder[typ]+'/hrawpt_leadjet')
        hrawpt_pu_leadjet[typ]=f.Get(folder[typ]+'/hrawpt_pu_leadjet')
        hrawpt_good_leadjet[typ]=f.Get(folder[typ]+'/hrawpt_good_leadjet')
        hrawpt_response_leadjet[typ]=f.Get(folder[typ]+'/hrawpt_response_leadjet')

        hpt_leadjet[typ]=f.Get(folder[typ]+'/hpt_leadjet')
        hpt_pu_leadjet[typ]=f.Get(folder[typ]+'/hpt_pu_leadjet')
        hpt_good_leadjet[typ]=f.Get(folder[typ]+'/hpt_good_leadjet')
        hpt_response_leadjet[typ]=f.Get(folder[typ]+'/hpt_response_leadjet')

        heta[typ]=f.Get(folder[typ]+'/heta')
        heta_pu[typ]=f.Get(folder[typ]+'/heta_pu')
        heta_good[typ]=f.Get(folder[typ]+'/heta_good')

        heta_leadjet[typ]=f.Get(folder[typ]+'/heta_leadjet')
        heta_pu_leadjet[typ]=f.Get(folder[typ]+'/heta_pu_leadjet')
        heta_good_leadjet[typ]=f.Get(folder[typ]+'/heta_good_leadjet')

        hmass[typ]=f.Get(folder[typ]+'/hmass')
        hmass_response[typ]=f.Get(folder[typ]+'/hmass_response')

        hmass_leadjet[typ]=f.Get(folder[typ]+'/hmass_leadjet')
        hmass_response_leadjet[typ]=f.Get(folder[typ]+'/hmass_response_leadjet')

        hnparticles[typ]=f.Get(folder[typ]+'/hnparticles')
        
    ## -----------------------------
    ## Plotting time
            
    makeComparisonPlots(hnjets, types, 'njets', linecolors,linestyles,False, False,outdir)        
    makeComparisonPlots(hrawpt, types, 'rawpt', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hpt, types, 'pt', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hrawpt_leadjet, types, 'rawpt_leadjet', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hpt_leadjet, types, 'pt_leadjet', linecolors,linestyles,False, False,outdir)

    makeComparisonPlots(hmass, types, 'mass', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hmass_leadjet, types, 'mass_leadjet', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(heta, types, 'eta', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(heta_leadjet, types, 'eta_leadjet', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hnparticles, types, 'nparticles', linecolors,linestyles,False, False,outdir)

    
    for typ in types:
        if typ!='GEN':
            makePileupPlots(hpt, hpt_pu, hpt_good,typ,'pt',outdir)
            makePileupPlots(hrawpt, hrawpt_pu, hrawpt_good,typ,'rawpt',outdir)
            makePileupPlots(heta, heta_pu, heta_good,typ,'eta',outdir)
            makePileupPlots(hpt_leadjet, hpt_pu_leadjet, hpt_good_leadjet,typ,'pt',outdir)
            makePileupPlots(hrawpt_leadjet, hrawpt_pu_leadjet, hrawpt_good_leadjet,typ,'rawpt',outdir)
            makePileupPlots(heta_leadjet, heta_pu_leadjet, heta_good_leadjet,typ,'eta',outdir)


    redtypes = types
    redtypes.remove('GEN')
    print redtypes
    
    makeComparisonPlots(hrawpt_response, redtypes, 'rawpt_response', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hrawpt_response_leadjet, redtypes, 'rawpt_response_leadjet', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hpt_response, redtypes, 'pt_response', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hpt_response_leadjet, redtypes, 'pt_response_leadjet', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hmass_response, redtypes, 'mass_response', linecolors,linestyles,False, False,outdir)
    makeComparisonPlots(hmass_response_leadjet, redtypes, 'mass_response_leadjet', linecolors,linestyles,False, False,outdir)
    

    raw_input('ok?')

        
        
