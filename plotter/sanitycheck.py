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

pulse_array = np.clip(pulse_array, 1e-6, None)  # Avoid zeros
pulse_array /= np.sum(pulse_array)
pulse_array = np.roll(pulse_array, int(nBins / 2) - int(np.argmax(pulse_array)))


def richardson_lucy(observed, psf, iterations=2, clip=True):
    psf = psf / np.sum(psf)
    psf_mirror = psf[::-1]
    estimate = np.full_like(observed, 0.5)
    for _ in range(iterations):
        conv = fftconvolve(estimate, psf, mode='same')
        conv = np.maximum(conv, 1e-8)
        ratio = observed / conv
        estimate *= fftconvolve(ratio, psf_mirror, mode='same')
        if clip:
            estimate = np.clip(estimate, 0, 1e6)
    return estimate

# Make a fake signal
truth = np.zeros(nBins)
truth[nBins // 2] = 1.0  # Delta function in the middle

# Convolve with your pulse shape
reco = fftconvolve(truth, pulse_array, mode='same')

# Deconvolve using the exact same pulse shape
deconv = richardson_lucy(reco, pulse_array, iterations=100)

# Compare
plt.plot(pulse_array, label="Pulse Shape")
# plt.plot(truth, label="Truth")
# plt.plot(reco, label="Convolved")
# plt.plot(deconv, label="Deconvolved")
plt.legend()
plt.savefig("deconv/test_deconv.png")
