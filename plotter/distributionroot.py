import ROOT
from collections import OrderedDict
import numpy as np

time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)
#nbin = 50

lf = ROOT.TF1("landau", "TMath::Landau(x, 5, .5)", 0, time_max)              #need to add "" around function in python to make it string
gf = ROOT.TF1("gauss", "TMath::Gaus(x, 5, .5, true)", 0, time_max)           #both are normalized i think



lanh = ROOT.TH1F("lanh", "Landau distribution", nBins, 0, time_max)          #name in root, title in plot, number of bins, xmin, xmax
gaush = ROOT.TH1F("gaush", "Gaussian distribution", nBins, 0, time_max)               

for bin in range(1, nBins + 1):
    x = lanh.GetBinCenter(bin)                                          #gets landau bin center
    y = gaush.GetBinCenter(bin)                                         #gets gaussian bin center

    lanh.SetBinContent(bin, lf.Eval(x))                                 #evaluates the landau function at the bin center
    gaush.SetBinContent(bin, gf.Eval(y))                                #evaluates the gaussian function at the bin center

# Code to check if distributions are correct
ofile = ROOT.TFile("distributions.root", "RECREATE")
lanh.Write()
gaush.Write()
ofile.Close()