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
os.makedirs("deconv", exist_ok=True)

# Parameters
time_start = 0
time_max = 32.0
time_per_bin = 0.03125
nBins = int(time_max / time_per_bin)

# Load and process pulse
g = ROOT.TFile("data_pulses.root")
pulse = g.Get("pulse_evt957")
pulse_array = np.array([pulse.GetBinContent(i + 1) for i in range(pulse.GetNbinsX())])

# if pulse.GetNbinsX() != nBins:
#     x_old = np.linspace(0, time_max, pulse.GetNbinsX())
#     x_new = np.linspace(0, time_max, nBins)
#     interp = interp1d(x_old, pulse_array, kind="cubic", fill_value=0, bounds_error=False)
#     pulse_array = interp(x_new)

pulse_array = np.clip(pulse_array, 1e-6, None)  # Avoid zeros
pulse_array /= np.sum(pulse_array)
pulse_array = np.roll(pulse_array, int(nBins / 2) - int(np.argmax(pulse_array)))

# Load data
f = ROOT.TFile("output.root")

def richardson_lucy(observed, psf, iterations=1000, clip=True):
    psf = psf / np.sum(psf)
    psf_mirror = psf[::-1]
    estimate = np.full_like(observed, 0.5)
    for _ in range(iterations):
        conv = fftconvolve(estimate, psf, mode='same')
        conv = np.maximum(conv, 1e-8)
        ratio = observed / conv
        estimate *= fftconvolve(ratio, psf_mirror, mode='same')
        if clip:
            estimate = np.clip(estimate, 0, 1e8)
    return estimate

# Load truth hits
truth_dict = OrderedDict()
for key in f.GetListOfKeys():
    name = key.GetName()
    if not name.startswith("h_truth"):
        continue
    match = re.match(r"h_truth_evt(\d+)_x(\d+)_y(\d+)", name)
    if not match:
        continue
    ievt, x, y = map(int, match.groups())
    hist = f.Get(name)
    truth_dict[(ievt, x, y)] = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])




# Process each reconstructed signal
for key in f.GetListOfKeys():
    name = key.GetName()
    if not name.startswith("h_reco"):
        continue
    match = re.match(r"h_reco_evt(\d+)_x(\d+)_y(\d+)", name)
    if not match:
        continue
    ievt, x, y = map(int, match.groups())
    hist = f.Get(name)
    signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])

    time_axis = np.linspace(time_start + time_per_bin / 2, time_max - time_per_bin / 2, nBins)

    smoothed = gaussian_filter1d(signal, sigma=2)
    noise_level = np.std(smoothed[:100])
    thresholded = np.where(smoothed > noise_level * 3, smoothed, 0)

    deconv = richardson_lucy(thresholded, pulse_array)

    peaks, _ = find_peaks(deconv, height=10)
    h_delta = ROOT.TH1F(f"h_delta_evt{ievt}_x{x}_y{y}", "Deconvolved Peaks", nBins, time_start, time_max)
    for p in peaks:
        h_delta.SetBinContent(int(p + 1), float(deconv[p]))
    delta = np.array([h_delta.GetBinContent(i + 1) for i in range(nBins)])

    truthhits = truth_dict.get((ievt, x, y), np.zeros(nBins))
    scalefactor = 50

    nonzero_indices = np.nonzero(truthhits)[0]
    if len(nonzero_indices) > 0:
        first_hit_time = time_axis[nonzero_indices[0]]
        pulse_peak_idx = np.argmax(pulse_array)
        pulse_peak_time = time_axis[pulse_peak_idx]
        shift = pulse_peak_time - first_hit_time
    else:
        shift = 0.0


    shifted_time_axis = time_axis - shift 
 
    plt.figure()
    plt.plot(shifted_time_axis, signal, label="Observed", alpha=0.5)
    plt.plot(shifted_time_axis, pulse_array * np.max(signal), '--', label="Response (scaled)")
    plt.plot(shifted_time_axis, deconv, label="Deconvolved", linewidth=1)
    plt.plot(shifted_time_axis, delta, label="Reconstructed Truth Hits", color="red")
    plt.plot(shifted_time_axis, truthhits * scalefactor, label="Truth Hits")
    plt.xlabel("Time Slice")
    plt.ylabel("Amplitude [arb.]")
    plt.xlim(0, time_max/2)
    # plt.ylim(0, 50)
    plt.title(f"Evt {ievt} | x{x} y{y}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"deconv/deconv_evt{ievt}_iX{x}_iY{y}.png")
    plt.close()
