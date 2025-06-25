import ROOT
from collections import OrderedDict
import numpy as np
import random

time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)
#nbin = 50
res = 0.1  

#Set Noise Stuff
mean = 0
sigma = 0.5

lf = ROOT.TF1("landau", f"TMath::Landau(x, 5, {res})", 0, time_max)              #need to add "" around function in python to make it string
gf = ROOT.TF1("gauss", f"TMath::Gaus(x, 5, {res}, true)", 0, time_max)           #both are normalized i think



lanh = ROOT.TH1F("lanh", "Landau distribution", nBins, 0, time_max)          #name in root, title in plot, number of bins, xmin, xmax
gaush = ROOT.TH1F("gaush", "Gaussian distribution", nBins, 0, time_max)               

for bin in range(1, nBins + 1):
    x = lanh.GetBinCenter(bin)                                          #gets landau bin center
    y = gaush.GetBinCenter(bin)                                         #gets gaussian bin center

    lanh.SetBinContent(bin, lf.Eval(x))                                 #evaluates the landau function at the bin center
    gaush.SetBinContent(bin, gf.Eval(y))                                #evaluates the gaussian function at the bin center

lpulse = np.zeros(lanh.GetNbinsX())
for i in range(lanh.GetNbinsX()):
    lpulse[i] = lanh.GetBinContent(i+1)

print("landau pulses: ", lpulse)

gpulse = np.zeros(gaush.GetNbinsX())
for i in range(gaush.GetNbinsX()):
    gpulse[i] = gaush.GetBinContent(i+1)

print("gaussian pulses: ", gpulse)

sim_data = "/home/catidmor/DREAMSim/sim/build/mc_mu0_run001_003_Test_10evt_mu+_100.0_100.0.root"
ifile = ROOT.TFile(sim_data)
tree = ifile.Get("tree")


# pulses for truth photons and reco with shapes
histos_truth = OrderedDict()

histos_lan_reco = OrderedDict()
histos_gaus_reco = OrderedDict()
histos_n_gaus_reco = OrderedDict()
histos_n_lan_reco = OrderedDict()

#creating 3x4 geometry (maybe) move up later
rodsize = .4

xsize = 3 * rodsize
ysize = 4 * rodsize

nevts = tree.GetEntries()
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()
    histos_lan_reco[ievt] = OrderedDict()
    histos_gaus_reco[ievt] = OrderedDict()
    histos_n_gaus_reco[ievt] = OrderedDict()
    histos_n_lan_reco[ievt] = OrderedDict()


def AddLandauPulse(l_pulse_reco, t0):
    t0_bin = l_pulse_reco.FindBin(t0)
    for i in range(lpulse.size):
        val = l_pulse_reco.GetBinContent(t0_bin + i)
        val += lpulse[i]
        l_pulse_reco.SetBinContent(t0_bin + i, val)

def AddGaussianPulse(g_pulse_reco, t0):
    t0_bin = g_pulse_reco.FindBin(t0)
    for i in range(gpulse.size):
        val = g_pulse_reco.GetBinContent(t0_bin + i)
        val += gpulse[i]
        g_pulse_reco.SetBinContent(t0_bin + i, val)

def AddNoise(n_pulse_reco, g_pulse_reco, t0):
    t0_bin = n_pulse_reco.FindBin(t0)
    for i in range(n_pulse_reco.GetNbinsX()):
        noise = random.gauss(mean, sigma)
        val = g_pulse_reco.GetBinContent(t0_bin + i)
        val += noise
        n_pulse_reco.SetBinContent(t0_bin + i, val)

def AddNoiseLandau(n_pulse_reco, l_pulse_reco, t0):
    t0_bin = n_pulse_reco.FindBin(t0)
    for i in range(n_pulse_reco.GetNbinsX()):
        noise = random.gauss(mean, sigma)
        val = l_pulse_reco.GetBinContent(t0_bin + i)
        val += noise
        n_pulse_reco.SetBinContent(t0_bin + i, val)

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
                f"h_truth_evt{ievt}_x{xscale}_y{yscale}", f"h_truth_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)

            histos_lan_reco[ievt][pulsehit] = ROOT.TH1D(
                f"landau_reco_evt{ievt}_x{xscale}_y{yscale}", f"landau_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
            histos_gaus_reco[ievt][pulsehit] = ROOT.TH1D(
                f"gaussian_reco_evt{ievt}_x{xscale}_y{yscale}", f"gaussian_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
            histos_n_gaus_reco[ievt][pulsehit] = ROOT.TH1D(
                f"gaus_noise_reco_evt{ievt}_x{xscale}_y{yscale}", f"gaus_noise_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
            histos_n_lan_reco[ievt][pulsehit] = ROOT.TH1D(
                f"lan_noise_reco_evt{ievt}_x{xscale}_y{yscale}", f"lan_noise_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
       

        # truth photon arriving time
        histos_truth[ievt][pulsehit].Fill(t, 1.0)
        
        # reco-ed pulse shapes for landau and gaussian distributions
        AddLandauPulse(histos_lan_reco[ievt][pulsehit], t)
        AddGaussianPulse(histos_gaus_reco[ievt][pulsehit], t)
        AddNoise(histos_n_gaus_reco[ievt][pulsehit], histos_gaus_reco[ievt][pulsehit], t)
        AddNoiseLandau(histos_n_lan_reco[ievt][pulsehit], histos_lan_reco[ievt][pulsehit], t)

    #break # for testing, only process the first event
        

# save to file
ofile = ROOT.TFile("output.root", "RECREATE")
for ievt in range(nevts):
    for pulsehit in histos_truth[ievt]:
        histos_truth[ievt][pulsehit].Write()
        #histos_reco[ievt][pulsehit].Scale(1.0 / histos_reco[ievt][pulsehit].Integral())
        histos_lan_reco[ievt][pulsehit].Write()
        histos_gaus_reco[ievt][pulsehit].Write()
        histos_n_gaus_reco[ievt][pulsehit].Write()
        histos_n_lan_reco[ievt][pulsehit].Write()
ofile.Close()


# Code to check if distributions are correct
# ofile = ROOT.TFile("distriution.root", "RECREATE")
# lanh.Write()
# gaush.Write()
# ofile.Close()