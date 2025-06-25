import sys
from collections import OrderedDict
import ROOT
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("plots", exist_ok=True)

sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

args = {
    'dology': False,
    'donormalize': False,
    'mycolors': [2, 3, 4, 6, 7, 8, 9, 28, 30, 38, 46],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True
}

time_max = 30.0
time_per_bin = 0.04
nBins = int(time_max / time_per_bin)
nevts = 10
xmin = 10
xmax = 15

res = 0.4
noiseres = 0.5

# gf = ROOT.TF1("gauss", f"TMath::Gaus(x, 5, {res}, true)", 0, time_max)
# gaush = ROOT.TH1F("gaush", f"Gaussian distribution with {res} ns Resolution", nBins, 0, time_max)

# for bin in range(1, nBins + 1):
#     y = gaush.GetBinCenter(bin)
#     gaush.SetBinContent(bin, gf.Eval(y))

# gpulse = np.zeros(gaush.GetNbinsX())
# for i in range(gaush.GetNbinsX()):
#     gpulse[i] = gaush.GetBinContent(i+1)

# plt.plot(gpulse)
# plt.xlabel("Time [ns]")
# plt.ylabel("Arbitrary Units")
# plt.title("Gaussian Pulse Shape")
# plt.savefig("plots/gpulse.png")
# plt.close()

fname_pulses = "output.root"
f_pulses = ROOT.TFile(fname_pulses)

# Plot truth photon time
histos_truth = OrderedDict()
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()
    for hit in f_pulses.GetListOfKeys():
        if "truth" in hit.GetName() and f"evt{ievt}" in hit.GetName():
            histos_truth[ievt][hit.GetName()] = f_pulses.Get(hit.GetName())

    plt.hist(list(histos_truth[ievt].values()), histtype='bar')
    plt.xlabel("Time [ns]")
    plt.ylabel("Voltage [mV]")
    plt.title(f"Truth Photon Time - Event {ievt}")
    plt.legend(list(histos_truth[ievt].keys()), loc='upper right')
    plt.savefig(f"plots/TruthPhotonTime_{ievt}.png")
    plt.close()

# Group and plot reco histograms by (evt, res, noise)
# grouped = {}
# for key in f_pulses.GetListOfKeys():
#     name = key.GetName().split(";")[0]
#     if "reco_evt" not in name:
#         continue

#     parts = name.split("_")
#     try:
#         evt = int([p for p in parts if p.startswith("evt")][0][3:])
#         res_val = float([p[3:] for p in parts if p.startswith("res")][0])
#         noise_val = float([p[5:] for p in parts if p.startswith("noise")][0])
#     except (IndexError, ValueError):
#         continue

#     grouped.setdefault((evt, res_val, noise_val), []).append(f_pulses.Get(name))

# for (evt, res_val, noise_val), hist_list in grouped.items():
#     plt.figure()
#     for hist in hist_list:
#         x = np.array([hist.GetBinCenter(i + 1) for i in range(hist.GetNbinsX())])
#         y = np.array([hist.GetBinContent(i + 1) for i in range(hist.GetNbinsX())])
#         plt.plot(x, y, label=hist.GetName())

#     plt.xlabel("Time [ns]")
#     plt.ylabel("Voltage [mV]")
#     plt.xlim(xmin, xmax)
#     plt.title(f"Reco Pulses - Event {evt}, Res {res_val}, Noise {noise_val}")
#     plt.legend(fontsize=6)
#     plt.savefig(f"plots/reco_evt{evt}_res{res_val}_noise{noise_val}.png")
#     plt.close()

noisy_pulses = OrderedDict()
for ievt in range(nevts):
    noisy_pulses[ievt] = OrderedDict()
    for hit in f_pulses.GetListOfKeys():
        if "h_noisy" in hit.GetName() and f"evt{ievt}" in hit.GetName():
            noisy_pulses[ievt][hit.GetName()] = f_pulses.Get(hit.GetName())

    plt.hist(list(noisy_pulses[ievt].values()), histtype='bar')
    plt.xlabel("Time [ns]")
    plt.ylabel("Voltage [mV]")
    plt.title(f"Noisy Pulses - Event {ievt}")
    plt.legend(list(noisy_pulses[ievt].keys()), loc='upper right')
    plt.savefig(f"plots/NoisyPulses_{ievt}.png")
    plt.close()


