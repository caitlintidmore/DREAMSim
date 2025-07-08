import ROOT
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import OrderedDict
import os

#parameters to make 1024 bins
time_max = 32.0
time_per_bin = 0.03125  
nBins = int(time_max / time_per_bin)

# #parameters to make 750 bins
# time_max = 30.0
# time_per_bin = 0.04  # use 40 ps per bin
# nBins = int(time_max / time_per_bin)

sim_data = "/home/catidmor/DREAMSim/sim/build/mc_mu0_run001_003_Test_10evt_mu+_100.0_100.0.root"
ifile = ROOT.TFile(sim_data)
tree = ifile.Get("tree")

#creating 3x4 geometry (maybe)
rodsize = .4
xsize = 3 * rodsize
ysize = 4 * rodsize

#real data
p_data = "data_pulses.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
h_pulse = pfile.Get("pulse_evt957") 


#attempting to make it with an actual pulse shape
pulses = np.zeros(nBins)
for i in range(nBins):
    pulses[i] = h_pulse.GetBinContent(i+1)

histos_truth = OrderedDict()
histos_reco = OrderedDict()

nevts = tree.GetEntries()
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()
    # histos_gaus_reco[ievt] = OrderedDict()
    # histos_n_gaus_reco[ievt] = OrderedDict()
    histos_reco[ievt] = OrderedDict()

def AddPulseFromArray(hist, pulse_array, t0):
    t0_bin = hist.FindBin(t0)
    for i in range(len(pulse_array)):
        bin_idx = t0_bin + i
        if 0 < bin_idx <= hist.GetNbinsX():
            val = hist.GetBinContent(bin_idx)
            hist.SetBinContent(bin_idx, val + pulse_array[i])
    

for ievt in range(nevts):
    tree.GetEntry(ievt)

    nPhotons = tree.nOPs
    for j in range(nPhotons):

        # make it to the end of the fiber
        if not tree.OP_pos_final_z[j] > 50.0:
            continue
        # look at photons from core of Cherenkov cone
        if not tree.OP_isCoreC[j]:
            continue

        t = tree.OP_time_final[j]
        #t = tree.OP_time_produced[j]
        #t = tree.OP_time_final[j] - tree.OP_time_produced[j]
       

        #breaking up photons hits into 3x4 grid
        x = tree.OP_pos_final_x[j]
        y = tree.OP_pos_final_y[j]


        xscale = int(x/(xsize))
        yscale = int(y/(ysize))

        pulsehit = (xscale, yscale)


        #makes histograms for pulsehits (creates new ones if it doesn't already exist)
        if pulsehit not in histos_truth[ievt]:
            histos_truth[ievt][pulsehit] = ROOT.TH1D(
                f"h_truth_evt{ievt}_x{xscale}_y{yscale}", f"h_truth_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
            histos_reco[ievt][pulsehit] = ROOT.TH1D(
                f"h_reco_evt{ievt}_x{xscale}_y{yscale}",
                f"h_reco_evt{ievt}_x{xscale}_y{yscale}",
                nBins, 0, time_max
                )



        # truth photon arriving time
        histos_truth[ievt][pulsehit].Fill(t, 1.0)
        
        #trying to reco with pulse shape and multiple noise levels
        hist = histos_reco[ievt][pulsehit]  # just one clean reco hist
        AddPulseFromArray(hist, pulses, t)


ofile = ROOT.TFile("output.root", "RECREATE")
for ievt in histos_truth:
    for pulsehit in histos_truth[ievt]:
        histos_truth[ievt][pulsehit].Write()
        histos_reco[ievt][pulsehit].Write()
ofile.Close()
