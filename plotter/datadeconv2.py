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
os.makedirs("RLdeconv", exist_ok=True)

# Parameters
time_start = 0
time_max = 30.0
time_per_bin = (30.0/1024)  # 0.029296875 ps
nBins = int(time_max / time_per_bin)

# Load and process pulse
g = ROOT.TFile("data_pulses.root")
pulse = g.Get("pulse_evt957")
pulse_array = np.array([pulse.GetBinContent(i + 1) for i in range(pulse.GetNbinsX())])

pulse_array = np.clip(pulse_array, 1e-6, None)  # Avoid zeros
pulse_array /= np.sum(pulse_array)
pulse_array = np.roll(pulse_array, int(nBins / 2) - int(np.argmax(pulse_array)))

def richardson_lucy(observed, psf, iterations=5, clip=True):
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
    
    ievt = map(int, match.groups())
    hist = f.Get(name)
  
    data = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])





    time_axis = np.linspace(time_start + time_per_bin / 2, time_max - time_per_bin / 2, nBins)

    smoothed = gaussian_filter1d(data, sigma=2)
    noise_level = np.std(smoothed[:100])
    thresholded = np.where(smoothed > noise_level * 3, smoothed, 0)

    deconv = richardson_lucy(thresholded, pulse_array)

    peaks, _ = find_peaks(deconv, height=10)
    h_delta = ROOT.TH1F(f"evt{ievt}", "Photon Hits", nBins, time_start, time_max)
    for p in peaks:
        #h_delta.SetBinContent(int(p + 1), float(deconv[p]))

        amp = float(data[p])
        h_delta.SetBinContent(int(p + 1), amp)

    delta = np.array([h_delta.GetBinContent(i + 1) for i in range(nBins)])

    time_axis = np.linspace(time_start + time_per_bin / 2, time_max - time_per_bin / 2, nBins)

    plt.figure()
    plt.plot(time_axis, data, label="Observed", alpha=0.5)
    plt.plot(time_axis, pulse_array * np.max(data), '--', label="Response (scaled)")
    plt.plot(time_axis, deconv, label="Deconvolved", linewidth=1)
    plt.plot(time_axis, delta, label="Photon Hits", color="red")
    plt.xlabel("Time Slice")
    plt.ylabel("Amplitude [arb.]")
    plt.xlim(0, time_max/2)
    # plt.ylim(0, 50)
    plt.title(f"Evt {ievt}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"RLdeconv/{ievt}.png")
    plt.close()
