import ROOT
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import OrderedDict
import os

os.makedirs("test", exist_ok=True)
noise = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
res = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

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


# for i in noise:
#         graph = np.random.normal(0, i, nBins)

        # plt.plot(graph)
        # plt.xlim(0, time_max)
        # plt.ylim(-2, 2)

        # plt.savefig(f"test/noise_res{i}.png")
        # plt.close()


# for i in res:
#         plt.figure()

#         gf = gf = ROOT.TF1("gauss", f"TMath::Gaus(x, 0, {res}, true)", -time_max / 2, time_max / 2)
#         gaush = ROOT.TH1F(f"gaush_res{i}", "Gaussian distribution", nBins, 0, time_max)

#         for bin in range(1, nBins + 1):
#                 x = gaush.GetBinCenter(bin)
#                 value = gf.Eval(x)
#                 gaush.SetBinContent(bin, value)

#         gpulse = np.zeros(gaush.GetNbinsX())
#         for j in range(gaush.GetNbinsX()):
#                 gpulse[j] = gaush.GetBinContent(j+1)

#         taxis = np.array([gaush.GetBinCenter(j + 1) for j in range(nBins)])

        # plt.plot(taxis, gpulse)
        # plt.xlim(0, time_max)
        # plt.ylim(-2, 6)
        
        # plt.savefig(f"test/gaus_res{i}.png")
        # plt.close()

#Going to add noise to gaussian and store in dictionary
#Should have all different noise levels for each resolution
#Hopefully </3

# noise_gaussians = OrderedDict()

# for i in res:
#     noise_gaussians[i] = OrderedDict()
#     gf = ROOT.TF1(f"gauss_res{i}", f"TMath::Gaus(x, 5, {i}, true)", 0, time_max)
#     gaush = ROOT.TH1F(f"gaush_res{i}", "Gaussian distribution", nBins, 0, time_max)

#     for bin in range(1, nBins + 1):
#         x = gaush.GetBinCenter(bin)
#         gaush.SetBinContent(bin, gf.Eval(x))

#     taxis = np.array([gaush.GetBinCenter(j + 1) for j in range(nBins)])
#     gpulse = np.array([gaush.GetBinContent(j+1) for j in range(nBins)])

#     for j in noise:
#         n_graph = np.random.normal(0, j, nBins)
#         ngaus = gpulse + n_graph  # vector addition
#         noise_gaussians[i][j] = ngaus

#         # plt.plot(taxis, ngaus)
#         # plt.xlim(0, time_max)
#         # plt.ylim(-2, 6)
#         # plt.title(f"Gaussian σ={i} + Noise σ={j}")
#         # plt.savefig(f"test/gaus_noise_res{i}_noise{j}.png")
#         # plt.close()

# for i in noise_gaussians:
#       for j in noise_gaussians[i]:
#         pulse = noise_gaussians[i][j]
#         print(f"\nres={i}, noise={j}, values=\n{pulse}")

noisy_pulse = OrderedDict()

p_data = "data_pulses.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
h_pulse = pfile.Get("pulse_evt957") 

#attempting to make it with an actual pulse shape
pulses = np.zeros(nBins)
for i in range(nBins):
    pulses[i] = h_pulse.GetBinContent(i+1)
 
# for i in noise:
#     n_graph = np.random.normal(0, i, nBins)
#     npulse = pulses + n_graph  # vector addition
#     noisy_pulse[i] = npulse


histos_truth = OrderedDict()
# histos_gaus_reco = OrderedDict()
# histos_n_gaus_reco = OrderedDict()



nevts = tree.GetEntries()
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()
    # histos_gaus_reco[ievt] = OrderedDict()
    # histos_n_gaus_reco[ievt] = OrderedDict()
    noisy_pulse[ievt] = OrderedDict()


# def AddGaussianPulse(g_pulse_reco, t0):
#     t0_bin = g_pulse_reco.FindBin(t0)
#     for i in range(gpulse.size):
#         val = g_pulse_reco.GetBinContent(t0_bin + i)
#         val += gpulse[i]
#         g_pulse_reco.SetBinContent(t0_bin + i, val)

# def AddNoise(n_pulse_reco, g_pulse_reco, t0):
#     t0_bin = n_pulse_reco.FindBin(t0)
#     for i in range(n_pulse_reco.GetNbinsX()):
#         noise = random.gauss(mean, sigma)
#         val = g_pulse_reco.GetBinContent(t0_bin + i)
#         val += noise
#         n_pulse_reco.SetBinContent(t0_bin + i, val)

def AddPulseFromArray(hist, pulse_array, t0):
    t0_bin = hist.FindBin(t0)

    for i in range(len(pulse_array)):
        bin_idx = t0_bin + i
        if 0 < bin_idx <= hist.GetNbinsX():
            val = hist.GetBinContent(bin_idx)
            hist.SetBinContent(bin_idx, val + pulse_array[i])


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
            
            noisy_pulse[ievt][pulsehit] = OrderedDict()
            for noise_val in noise:
                noisy_pulse[ievt][pulsehit][noise_val] = ROOT.TH1D(
                f"h_noisy_pulse_evt{ievt}_x{xscale}_y{yscale}_noise{noise_val}",
                f"h_noisy_pulse_evt{ievt}_x{xscale}_y{yscale}_noise{noise_val}",
                nBins, 0, time_max
            )
            # histos_gaus_reco[ievt][pulsehit] = ROOT.TH1D(
            #     f"gaussian_reco_evt{ievt}_x{xscale}_y{yscale}", f"gaussian_reco_evt{ievt}_x{xscale}_y{yscale}", nBins, 0, time_max)
            
            # if pulsehit not in histos_n_gaus_reco[ievt]:
            #     histos_n_gaus_reco[ievt][pulsehit] = OrderedDict()

            # for i in res:
            #     histos_n_gaus_reco[ievt][pulsehit][i] = OrderedDict()

            #     for j in noise:
            #          histos_n_gaus_reco[ievt][pulsehit][i][j] = ROOT.TH1D(
            #             f"reco_evt{ievt}_x{xscale}_y{yscale}_res{i}_noise{j}", "",
            #             nBins, 0, time_max)


       

        # truth photon arriving time
        histos_truth[ievt][pulsehit].Fill(t, 1.0)
        
        #trying to reco with pulse shape and multiple noise levels
        for noise_val in noise:
            pulse_array = noisy_pulse[noise_val]
            hist = noisy_pulse[ievt][pulsehit][noise_val]
            AddPulseFromArray(hist, pulse_array, t)




        # reco-ed pulse shapes for landau and gaussian distributions
        # for res_val in res:
        #     for noise_val in noise:
        #         pulse_array = noise_gaussians[res_val][noise_val]
        #         hist = histos_n_gaus_reco[ievt][pulsehit][res_val][noise_val]
        #         AddPulseFromArray(hist, pulse_array, t-5)
        # AddGaussianPulse(histos_gaus_reco[ievt][pulsehit], t)
        # AddNoise(histos_n_gaus_reco[ievt][pulsehit], histos_gaus_reco[ievt][pulsehit], t)

    break # for testing, only process the first event
        

# save to file
# ofile = ROOT.TFile("output.root", "RECREATE")
# for ievt in range(nevts):
#     for pulsehit in histos_truth[ievt]:
#         histos_truth[ievt][pulsehit].Write()
#         #histos_reco[ievt][pulsehit].Scale(1.0 / histos_reco[ievt][pulsehit].Integral())
#         histos_gaus_reco[ievt][pulsehit].Write()
#         histos_n_gaus_reco[ievt][pulsehit].Write()
# ofile.Close()

ofile = ROOT.TFile("output.root", "RECREATE")
for ievt in histos_truth:
    for pulsehit in histos_truth[ievt]:
        histos_truth[ievt][pulsehit].Write()
        for noise_val in noisy_pulse[ievt][pulsehit]:
            noisy_pulse[ievt][pulsehit][noise_val].Write()        
        # for res_val in res:
        #     for noise_val in noise:
        #         histos_n_gaus_reco[ievt][pulsehit][res_val][noise_val].Write()
ofile.Close()