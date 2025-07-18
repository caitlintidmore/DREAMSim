import ROOT
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.signal import fftconvolve, find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os
import re

# Create output directory
os.makedirs("byhandsdeconv", exist_ok=True)

# Parameters
time_start = 0
time_max = 30.0
time_per_bin = (30.0/1024)  # 0.029296875 ps
nBins = int(time_max / time_per_bin)



# Load data
f = ROOT.TFile("drsoutput.root")



# Process each reconstructed signal
for key in f.GetListOfKeys():
    name = key.GetName()
    if not name.startswith("waveform"):
        continue
    match = re.match(r"waveform_(\d+)", name)
    if not match:
        continue
    
    hist = f.Get(name)
  
    data = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
    # smoothed = gaussian_filter1d(data, sigma=3)

    #next two lines look at first 100 bins, estimate noise level, and 
    #finds peaks 3x above noise, makes everything else 0
    
    # noise_level = np.std(smoothed[:100])
    # threshold = np.where(smoothed > noise_level * 3, smoothed, 0)

    noise_lvl = np.std(data[:100])
    threshld = np.where(data > noise_lvl * 3, data, 0)


    #peak finding for raw and smoothed data
    peaks, _ = find_peaks(threshld)              #_ = dummy variable for peak properties (i.e height,width,etc)
    # smoothed_peaks, _ = find_peaks(threshold)

    #deltas for regular data
    h_delta = ROOT.TH1F(f"h_delta_evt{name}", "Deconvolved Peaks", nBins, time_start, time_max)
    for p in peaks:
        h_delta.SetBinContent(int(p + 1), float(threshld[p]))
    delta = np.array([h_delta.GetBinContent(i + 1) for i in range(nBins)])  


    for i, val in enumerate(delta):
        if val != 0:
            print(f"Name: {name}, Bin: {i}, Value: {val}")


    # #deltas for smoothed data
    # h_delta_smoothed = ROOT.TH1F(f"h_delta_smoothed_evt{name}", "Deconvolved Smoothed Peaks", nBins, time_start, time_max)
    # for p in smoothed_peaks:
    #     h_delta_smoothed.SetBinContent(int(p + 1), float(threshold[p]))
    # delta_smoothed = np.array([h_delta_smoothed.GetBinContent(i + 1) for i in range(nBins)]) 



    plt.figure()
    plt.plot(data, label="Observed Signal", alpha=0.5)
    # plt.plot(smoothed, label="Smoothed Signal", linewidth=1)
    # plt.plot(threshold, label="Smoothed Thresholded Signal", linewidth=1)
    # plt.plot(threshld, label="Raw Thresholded Signal", linewidth=1)
    plt.plot(delta, label="Photon Hits", color="red", linewidth=1)
    # plt.plot(delta_smoothed, label="Deconvolved Smoothed Peaks", color="orange", linewidth=1)
    plt.xlabel("Time Slice")
    plt.savefig(f"byhanddeconv/{key}.png")
    plt.close()




    #print(f"{name}'s max: {np.max(smoothed)}, noise level: {noise_level}")

    #if np.max(smoothed) < 3*noise_level:
        #print(f"Skipping {name} due to low signal")













    # # time_axis = np.linspace(time_start + time_per_bin / 2, time_max - time_per_bin / 2, nBins)





    # plt.figure()
    # #plt.plot(time_axis, signal, label="Observed", alpha=0.5)
    # plt.plot(time_axis, pulse_array * np.max(signal), '--', label="Response (scaled)")
    # #plt.plot(time_axis, deconv, label="Deconvolved", linewidth=1)
    # #plt.plot(time_axis, delta, label="Reconstructed Truth Hits", color="red")
    # plt.xlabel("Time Slice")
    # plt.ylabel("Amplitude [arb.]")
    # plt.xlim(0, time_max/2)
    # # plt.ylim(0, 50)
    # plt.title(f"Evt {ievt}")
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig(f"deconv/deconv_evt{ievt}.png")
    # plt.close()
