import ROOT
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import OrderedDict
import os

os.makedirs("test", exist_ok=True)
noise = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
res = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)

sim_data = "/home/catidmor/DREAMSim/sim/build/mc_mu0_run001_003_Test_10evt_mu+_100.0_100.0.root"
ifile = ROOT.TFile(sim_data)
tree = ifile.Get("tree")

#creating 3x4 geometry (maybe)
rodsize = .4
xsize = 3 * rodsize
ysize = 4 * rodsize


p_data = "data_pulses.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
h_pulse = pfile.Get("pulse_evt957")

print("Number of bins in h_pulse:", h_pulse.GetNbinsX())


#attempting to make it with an actual pulse shape
nPulseBins = h_pulse.GetNbinsX()
full_pulses = np.array([h_pulse.GetBinContent(i+1) for i in range(nPulseBins)])
 

# peak = np.argmax(full_pulses)  
# start = max(0, peak - 40)
# end = min(len(full_pulses), peak + 100)
# pulses = full_pulses[start:end]


# pulses = full_pulses[start:end]


histos_truth = OrderedDict()
histos_reco = OrderedDict()

nevts = tree.GetEntries()

# def AddPulse(h_pulse_reco, t0):
#     t0_bin_center = h_pulse_reco.FindBin(t0)
#     bin_center_time = h_pulse_reco.GetBinCenter(t0_bin_center)
#     center = pulses.size // 2

#     shift_in_bins = int(round((t0 - bin_center_time) / time_per_bin))
#     t0_bin = t0_bin_center + shift_in_bins + 5

#     for i in range(pulses.size):
#         bin_idx = t0_bin + i - center
#         if 1 <= bin_idx <= h_pulse_reco.GetNbinsX():
#             val = h_pulse_reco.GetBinContent(bin_idx)
#             val += pulses[i]
#             h_pulse_reco.SetBinContent(bin_idx, val)

#     print(f"  t0 = {t0:.2f}, t0_bin = {t0_bin}, pulse center = {center}")


def AddPulse(h_pulse_reco, t0):
    t0_bin = h_pulse_reco.FindBin(t0)
    for i in range(full_pulses.size):
        val = h_pulse_reco.GetBinContent(t0_bin + i)
        val += full_pulses[i]
        h_pulse_reco.SetBinContent(t0_bin + i, val)


nevts = tree.GetEntries()
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

        #print(f"Produced: {tree.OP_time_produced[j]:.2f}, Final: {tree.OP_time_final[j]:.2f}, Δt: {t:.2f}")

        if ievt not in histos_truth:
            histos_truth[ievt] = OrderedDict()
            histos_reco[ievt] = OrderedDict()

        #makes histograms for pulsehits (creates new ones if it doesn't already exist)
        if pulsehit not in histos_truth[ievt]:
            histos_truth[ievt][pulsehit] = ROOT.TH1D(
                f"h_truth_evt{ievt}_x{xscale}_y{yscale}", f"h_truth_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            histos_reco[ievt][pulsehit] = ROOT.TH1D(
                f"h_reco_evt{ievt}_x{xscale}_y{yscale}", f"h_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
  
        # truth photon arriving time
        histos_truth[ievt][pulsehit].Fill(t, 1.0)
        # reco-ed pulse shape
        AddPulse(histos_reco[ievt][pulsehit], t)

    # break # for testing, only process the first event
        

ofile = ROOT.TFile("output.root", "RECREATE")
for ievt in histos_truth:
    for pulsehit in histos_truth[ievt]:
        histos_truth[ievt][pulsehit].Write()
        histos_reco[ievt][pulsehit].Write()
ofile.Close()