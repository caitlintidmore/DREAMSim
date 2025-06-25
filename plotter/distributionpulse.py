import ROOT
from collections import OrderedDict
import numpy as np

p_data = "distributions.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
#g_pulse = pfile.Get("gaush")
l_pulse = pfile.Get("lanh")

pulses = np.zeros(l_pulse.GetNbinsX())
for i in range(l_pulse.GetNbinsX()):
    pulses[i] = l_pulse.GetBinContent(i+1)

# pulses = np.zeros(l_pulse.GetNbinsX())
# for i in range(g_pulse.GetNbinsX()):
#     pulses[i] = g_pulse.GetBinContent(i+1)

print("pulses: ", pulses)

sim_data = "/home/catidmor/DREAMSim/sim/build/mc_mu45_run001_003_Test_10evt_mu+_100.0_100.0.root"
ifile = ROOT.TFile(sim_data)
tree = ifile.Get("tree")

time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)

# pulses for truth photons and reco with shapes
histos_truth = OrderedDict()
histos_reco = OrderedDict()

#creating 3x4 geometry (maybe)
rodsize = .4

xsize = 3 * rodsize
ysize = 4 * rodsize


nevts = tree.GetEntries()
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()
    histos_reco[ievt] = OrderedDict()

def AddPulse(h_pulse_reco, t0):
    t0_bin = h_pulse_reco.FindBin(t0)
    for i in range(pulses.size):
        val = h_pulse_reco.GetBinContent(t0_bin + i)
        val += pulses[i]
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


        #makes histograms for pulsehits (creates new ones if it doesn't already exist)
        if pulsehit not in histos_truth[ievt]:
            histos_truth[ievt][pulsehit] = ROOT.TH1D(
                f"tevt{ievt}_x{xscale}y{yscale}", f"tevt{ievt}_x{xscale}y{yscale}", nBins, 0, time_max)

            histos_reco[ievt][pulsehit] = ROOT.TH1D(
                f"revt{ievt}_x{xscale}y{yscale}", f"revt{ievt}_x{xscale}y{yscale}", nBins, 0, time_max)
       


        # truth photon arriving time
        histos_truth[ievt][pulsehit].Fill(t, 1.0)
        # reco-ed pulse shape
        AddPulse(histos_reco[ievt][pulsehit], t)



# save to file
ofile = ROOT.TFile("output.root", "RECREATE")
for ievt in range(nevts):
    for pulsehit in histos_truth[ievt]:
        histos_truth[ievt][pulsehit].Write()
        #histos_reco[ievt][pulsehit].Scale(1.0 / histos_reco[ievt][pulsehit].Integral())
        histos_reco[ievt][pulsehit].Write()
ofile.Close()