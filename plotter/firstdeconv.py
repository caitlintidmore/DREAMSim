#via chatGPT

import ROOT
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.signal import fftconvolve
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.interpolate import interp1d


import os
import re


f = ROOT.TFile("data_pulses.root")
pulse = f.Get("pulse_evt957")
bins = pulse.GetNbinsX()
print(f"Number of bins in data_pulses.root: {bins}")



# # Create output directory lol
# os.makedirs("deconv", exist_ok=True)

# f = ROOT.TFile("output.root")


# # Parameters
# time_max = 30.0                                 # Full window: 30 ns
# time_per_bin = 0.04                             # 40 ps per bin
# nBins = int(time_max / time_per_bin)
# time_start = 0                                  # Start at 0 ns
# time_end = time_max                             # End at +30 ns
# # nevts = 1                                     # Number of events to process


# def get_gaussian_pulse(res):
#     center = time_max / 2  # center = 15.0 ns
#     gf = ROOT.TF1("gauss", f"TMath::Gaus(x, {center}, {res}, true)", 0, time_max)
#     h = ROOT.TH1F("template", "template", nBins, 0, time_max)
#     for i in range(1, nBins + 1):
#         h.SetBinContent(i, gf.Eval(h.GetBinCenter(i)))
#     return np.array([h.GetBinContent(i + 1) for i in range(nBins)])


# # Richardson-Lucy deconvolution
# def richardson_lucy(observed, psf, iterations=50, clip=True):
#     psf = psf / np.sum(psf)
#     psf_mirror = psf[::-1]                                       #reverses the list (in this case psf); makes a mirror of the psf, needed for deconvolution
    
#     estimate = np.full_like(observed, 0.5)                       #creates new array with the same shape/data type as observed, filled with 0.5
#                                                                  #initial estimate for the deconvolution, can be any value, here we use 0.5
    
#     for _ in range(iterations):
#         conv = fftconvolve(estimate, psf, mode='same')
#         conv = np.maximum(conv, 1e-8)
#         ratio = observed / conv
#         estimate *= fftconvolve(ratio, psf_mirror, mode='same')
#         if clip:
#             estimate = np.clip(estimate, 0, 1e6)

#     return estimate

# truth_dict = OrderedDict()

# for key in f.GetListOfKeys():
#     tname = key.GetName()
#     if not tname.startswith("h_truth"):
#         continue

#     tmatch = re.match(r"h_truth_evt(\d+)_x(\d+)_y(\d+)", tname)

#     if not tmatch:
#         continue

#     tievt, tx, ty = tmatch.groups()
#     tx = int(tx)
#     ty = int(ty)
#     tievt = int(tievt)

#     thist = f.Get(tname)
#     if not thist:
#         continue

#     truth_dict[(tievt, tx, ty)] = np.array([thist.GetBinContent(i + 1) for i in range(nBins)])


# res = 4.0 # Resolution of the Gaussian pulse shape
# for key in f.GetListOfKeys():
#     name = key.GetName()

#     if not name.startswith("h_reco"):
#         continue

#     match = re.match(r"h_reco_evt(\d+)_x(\d+)_y(\d+)", name)

#     if not match:
#         continue

#     ievt, x, y = match.groups()
#     x = int(x)
#     y = int(y)
#     ievt = int(ievt)

#     # if ievt >= nevts:
#     #     continue

#     hist = f.Get(name)
#     if not hist:
#         continue

#     signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
#     smoothed = gaussian_filter1d(signal, sigma=5.0)


#     # Extract signal and time axis
#     # signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
#     time_axis = np.linspace(time_start + time_per_bin / 2, time_end - time_per_bin / 2, nBins)
#     # shifted_time_axis = time_axis + truth_shift


#     # Build response and apply deconvolution
#     response = get_gaussian_pulse(res)
#     # response = pulses
 


#     # deconv = richardson_lucy(signal, response)
#     deconv = richardson_lucy(smoothed, response)

#     peaks, properties = find_peaks(deconv, height=10)

#     #Rough estimation of delta function
#     #gf = ROOT.TF1("gauss", f"TMath::Gaus(x, {time_max}/2 , 0.02, true)", 0, time_max)

#     h_delta = ROOT.TH1F(f"h_delta_evt{ievt}_x{x}_y{y}", "Deconvolved Peaks", nBins, time_start, time_end)

#     for p in peaks:
#         t = time_axis[p]  # Convert peak index to time
#         h_delta.Fill(t, deconv[p])

#     delta = np.array([h_delta.GetBinContent(i + 1) for i in range(nBins)])
    
#     # idelta = delta.astype(int)

#     # print(f"Event: {ievt} delta array values: {delta}")
#     # print(f"Event: {ievt} delta array integer values: {idelta}")
    
#     truthhits = truth_dict.get((ievt, x, y), np.zeros(nBins))
#     scalefactor = 50

#     nonzero_indices = np.nonzero(truthhits)[0]
#     if len(nonzero_indices) > 0:
#         first_hit_time = time_axis[nonzero_indices[0]]
#         shift = (time_max / 2) - first_hit_time  
#     else:
#         shift = 0.0

#     shifted_time_axis = time_axis - shift  



#     # Plot
#     plt.figure()
#     plt.plot(shifted_time_axis, smoothed, label="Observed Smoothed", alpha=0.5)
#     # plt.plot(shifted_time_axis, signal, label="Observed", alpha=0.5)
#     # plt.plot(shifted_time_axis, response * np.max(signal), '--', label="Response (scaled)")
#     plt.plot(shifted_time_axis, deconv, label="Deconvolved", linewidth=1)
#     plt.plot(time_axis, truthhits * scalefactor, label="Truth Hits")
#     plt.plot(shifted_time_axis, delta, label="Reconstructed Truth Hits")
#     plt.xlabel("Time Slice")
#     plt.ylabel("Amplitude [arb.]")
#     plt.xlim(0, time_max)
#     plt.title(f"Evt {ievt} | x{x} y{y}")
#     plt.legend()
#     plt.grid(True)
#     plt.tight_layout()
#     plt.savefig(f"deconv/deconv_evt{ievt}_iX{x}_iY{y}.png")
#     plt.close()












# #dirac delta jazz
# center = 0
# x = np.arange(-nBins//2, nBins//2)

# def get_dirac_delta_pulse(i):
#     if i == center:
#         return 1.0
#     else:
#         return 0.0
    
# pulse = [get_dirac_delta_pulse(i) for i in x]


# plt.plot(x, pulse)
# plt.title("Dirac Delta Pulse")
# plt.xlim(-10, 10)
# plt.ylim(-0.1, 1.1)
# plt.savefig("deconv/dirac_delta_pulse.png")
# plt.close()


# Load input file
# f = ROOT.TFile("output.root")

# # Loop over all matching histograms
# for key in f.GetListOfKeys():
#     name = key.GetName()
#     if not name.startswith("reco_evt"):
#         continue

#     match = re.match(r"reco_evt(\d+)_x(\d+)_y(\d+)_res([0-9.]+)_noise([0-9.]+)", name)      #matches pattern with the string (name)
#     if not match:
#         continue                                                                            #if not name properly, continue

#     ievt, x, y, res, noise = match.groups()
#     ievt = int(ievt)
#     res = float(res)
#     noise = float(noise)

#     if ievt >= nevts:                                                                       #if ievt is greater than the number of events skip
#         continue

#     hist = f.Get(name)                                                                      #get the histogram from the file  
#     if not hist:                                                                            #if the histogram is not found, continue    
#         continue

#     # Extract signal and time axis
#     signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
#     time_axis = np.linspace(time_start + time_per_bin / 2, time_end - time_per_bin / 2, nBins)

#     # Build response and apply deconvolution
#     response = get_gaussian_pulse(res)
#     deconv = richardson_lucy(signal, response)

#     # Plot
#     plt.figure()
#     plt.plot(time_axis, signal, label="Observed", alpha=0.5)
#     # plt.plot(time_axis, response * np.max(signal), '--', label="Response (scaled)")
#     # plt.plot(time_axis, deconv, label="Deconvolved", linewidth=2)
#     plt.xlabel("Time [ns]")
#     plt.ylabel("Amplitude [arb.]")
#     plt.title(f"Evt {ievt} | x{x} y{y} | res={res} noise={noise}")
#     plt.legend()
#     plt.grid(True)
#     plt.tight_layout()
#     plt.savefig(f"deconv/deconv_evt{ievt}_x{x}_y{y}_res{res}_noise{noise}.png")
#     plt.close()