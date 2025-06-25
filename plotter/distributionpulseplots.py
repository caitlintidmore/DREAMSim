import sys                                  #talk to os
from collections import OrderedDict         #for neatness
import ROOT
sys.path.append("./CMSPLOTS")               #puts cmsplots in path to access it anywhere
from myFunction import DrawHistos
import matplotlib.pyplot as plt
import numpy as np

args = {
    'dology': False,
    'donormalize': False,
    'mycolors': [1],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True,
    'mycolors': [2, 3, 4, 6, 7, 8, 9, 28, 30, 38, 46]
}
time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)
#nbin = 5
nevts = 10
xmin = 10
xmax = 15
ymin = 0
ymax = 6
res = 0.4  # resolution of the landau and gaussian distributions
noiseres = 0.5

# Drawing landau and gaussian pulse shapes
#resolution is width of distribution REMEMBER!!
#lf = ROOT.TF1("landau", f"TMath::Landau(x, 5, {res})", 0, time_max)                   #need to add "" around function in python to make it string
gf = ROOT.TF1("gauss", f"TMath::Gaus(x, 5, {res}, true)", 0, time_max)                #both are normalized i think



#lanh = ROOT.TH1F("lanh", f"Landau distribution with {res} ns Resolution", nBins, 0, time_max)          #name in root, title in plot, number of bins, xmin, xmax
gaush = ROOT.TH1F("gaush", f"Gaussian distribution with {res} ns Resolution", nBins, 0, time_max)               

for bin in range(1, nBins + 1):
    #x = lanh.GetBinCenter(bin)                                               #gets landau bin center
    y = gaush.GetBinCenter(bin)                                              #gets gaussian bin center

    #lanh.SetBinContent(bin, lf.Eval(x))                                      #evaluates the landau function at the bin center
    gaush.SetBinContent(bin, gf.Eval(y))                                     #evaluates the gaussian function at the bin center

#lpulse = np.zeros(lanh.GetNbinsX())
#for i in range(lanh.GetNbinsX()):
    #lpulse[i] = lanh.GetBinContent(i+1)

#print("landau pulses: ", lpulse)

gpulse = np.zeros(gaush.GetNbinsX())
for i in range(gaush.GetNbinsX()):
    gpulse[i] = gaush.GetBinContent(i+1)

# plt.plot(lpulse)
# plt.xlabel("Time [ns]")
# plt.ylabel("Arbitrary Units")
# plt.title("Landau Pulse Shape")
# plt.savefig("plots/lpulse.png")
# plt.close()

plt.plot(gpulse)
plt.xlabel("Time [ns]")
plt.ylabel("Arbitrary Units")
plt.title("Gaussian Pulse Shape")
plt.savefig("plots/gpulse.png")

plt.close()

######################################################################             

# pulses for truth photons
fname_pulses = "output.root"
f_pulses = ROOT.TFile(fname_pulses)

args['dology'] = False

# truth photon time
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


   # DrawHistos(list(histos_truth[ievt].values()), list(histos_truth[ievt].keys()), 0, 20, "Measured Time [ns]",
    #            0, 25, "Voltage [mV]", f"Truth_{ievt}", **args)

######################################################################

#Pulses for reco photons

#LOG SCALE THING
args['dology'] = True

# reco photon time
#histos_lan_reco = OrderedDict()
histos_gaus_reco = OrderedDict()


# for ievt in range(nevts):
#     histos_lan_reco[ievt] = OrderedDict()
#     for hit in f_pulses.GetListOfKeys():
#         print(hit.GetName())
#         name = hit.GetName().split(";")[0]                                                  #Removes ;1 to hopefully fix my bug
#         if "landau_reco" in name and f"evt{ievt}" in name:
#             histos_lan_reco[ievt][name] = f_pulses.Get(hit.GetName())
        

#     for name, hist in histos_lan_reco[ievt].items():
#         n_bins = hist.GetNbinsX()
#         x = np.array([hist.GetBinCenter(i+1) for i in range(n_bins)])
#         y = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])
        
#         #print(f"Event {ievt} — matched pulses: {list(histos_lan_reco[ievt].keys())}")

#         print(f"Adding plot for {name}")

#         plt.plot(x, y, label=name)

#     plt.xlabel("Time [ns]")
#     plt.ylabel("Voltage [mV]")
#     plt.xlim(xmin, xmax)
#     #plt.ylim(ymin, ymax)
#     plt.title(f"Landau Reco Photon Time - Event {ievt}")

#     print("Legend entries:", list(histos_lan_reco[ievt].keys()))

#     #plt.legend(list(histos_lan_reco[ievt].keys()), loc='upper right')
#     plt.legend(f"Distribution Resolution: {res}", loc='upper left')
#     plt.savefig(f"plots/LandauReco_{ievt}.png")
#     plt.close()


for ievt in range(nevts):
    histos_gaus_reco[ievt] = OrderedDict()
    for hit in f_pulses.GetListOfKeys():
        print(hit.GetName())
        name = hit.GetName().split(";")[0]                                                  #Removes ;1 to hopefully fix my bug
        if "reco" in name and f"evt{ievt}" in name:
            histos_gaus_reco[ievt][name] = f_pulses.Get(hit.GetName())
        

    for name, hist in histos_gaus_reco[ievt].items():
        n_bins = hist.GetNbinsX()
        x = np.array([hist.GetBinCenter(i+1) for i in range(n_bins)])
        y = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])
        
        #print(f"Event {ievt} — matched pulses: {list(histos_gaus_reco[ievt].keys())}")

        print(f"Adding plot for {name}")

        plt.plot(x, y, label=name)


    # ng_rods = OrderedDict()
    # rod_count = 0

    # for name in histos_gaus_reco[ievt].keys():
    #     name = name.split("_")
    #     xpos = [p for p in name if p.startswith("x")][0]
    #     ypos = [p for p in name if p.startswith("y")][0]
    #     rod = f"{xpos}_{ypos}"

    #     ng_rods[rod_count] = rod
    #     rod_count += 1

    plt.xlabel("Time [ns]")
    plt.ylabel("Voltage [mV]")
    plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    plt.title(f"Gaussian Reco Photon Time - Event {ievt}")

    print("Legend entries:", list(histos_gaus_reco[ievt].keys()))

    plt.plot([], [], label=f"Distribution Resolution: {res}")
    plt.legend(loc='upper left')

    plt.legend(list(histos_gaus_reco[ievt].keys()), loc='upper right')
    plt.savefig(f"plots/GaussianReco_{ievt}.png")
    plt.close()

    ##################################################################
#Now add noise to the reco pulses
histos_n_gaus_reco = OrderedDict()

for ievt in range(nevts):
    histos_n_gaus_reco[ievt] = OrderedDict()
    for hit in f_pulses.GetListOfKeys():
        #print(hit.GetName())
        name = hit.GetName().split(";")[0]                                                  #Removes ;1 to hopefully fix my bug
        if "gaus_noise_reco" in name and f"evt{ievt}" in name:
            histos_n_gaus_reco[ievt][name] = f_pulses.Get(hit.GetName())
        

    for name, hist in histos_n_gaus_reco[ievt].items():
        n_bins = hist.GetNbinsX()
        x = np.array([hist.GetBinCenter(i+1) for i in range(n_bins)])
        y = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])
        
        #print(f"Event {ievt} — matched pulses: {list(histos_lan_reco[ievt].keys())}")

        print(f"Adding plot for {name}")

        plt.plot(x, y, label=name)

    plt.xlabel("Time [ns]")
    plt.ylabel("Voltage [mV]")
    plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    plt.title(f"Gaussian Reco with Noise- Event {ievt}")

    #print("Legend entries:", list(histos_n_gaus_reco[ievt].keys()))

    plt.text(0.05, 0.95, f"Noise Resolution: {noiseres}", fontsize=10, transform=plt.gca().transAxes, verticalalignment='top')
    plt.legend(list(histos_n_gaus_reco[ievt].keys()), loc='upper right')
   
    plt.savefig(f"plots/NoiseGaussianReco_{ievt}.png")
    plt.close()


histos_n_lan_reco = OrderedDict()

# for ievt in range(nevts):
#     histos_n_lan_reco[ievt] = OrderedDict()
#     for hit in f_pulses.GetListOfKeys():
#         #print(hit.GetName())
#         name = hit.GetName().split(";")[0]                                                  #Removes ;1 to hopefully fix my bug
#         if "lan_noise_reco" in name and f"evt{ievt}" in name:
#             histos_n_lan_reco[ievt][name] = f_pulses.Get(hit.GetName())
        

#     for name, hist in histos_n_lan_reco[ievt].items():
#         n_bins = hist.GetNbinsX()
#         x = np.array([hist.GetBinCenter(i+1) for i in range(n_bins)])
#         y = np.array([hist.GetBinContent(i+1) for i in range(n_bins)])
        
#         #print(f"Event {ievt} — matched pulses: {list(histos_lan_reco[ievt].keys())}")

#         print(f"Adding plot for {name}")

#         plt.plot(x, y, label=name)

#     plt.xlabel("Time [ns]")
#     plt.ylabel("Voltage [mV]")
#     plt.xlim(xmin, xmax)
#     #plt.ylim(ymin, ymax)
#     plt.title(f"Landau Reco with Noise- Event {ievt}")

    #print("Legend entries:", list(histos_n_gaus_reco[ievt].keys()))

    plt.text(0.05, 0.95, f"Noise Resolution: {noiseres}", fontsize=10, transform=plt.gca().transAxes, verticalalignment='top')
    plt.legend(list(histos_n_gaus_reco[ievt].keys()), loc='upper right')
   
    plt.savefig(f"plots/NoiseLandauReco_{ievt}.png")
    plt.close()