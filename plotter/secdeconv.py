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

# Create output directory lol
os.makedirs("deconv", exist_ok=True)

# Parameters - 1024 bins
time_start = 0
time_max = 32.0                                 # Full window: 32 ns
time_per_bin = 0.03125                          # 31.25 ps per bin
nBins = int(time_max / time_per_bin)


g = ROOT.TFile("data_pulses.root")
pulse = g.Get("pulse_evt957")
pulse_array = np.array([pulse.GetBinContent(i+1) for i in range(pulse.GetNbinsX())])

# If not already 1024 bins, interpolate to match nBins
if pulse.GetNbinsX() != nBins:
    x_old = np.linspace(0, time_max, pulse.GetNbinsX())
    x_new = np.linspace(0, time_max, nBins)
    interp = interp1d(x_old, pulse_array, kind="cubic", fill_value=0, bounds_error=False)
    pulse_array = interp(x_new)

# Optional: normalize
pulse_array /= np.sum(pulse_array)

# Optional: center pulse peak to middle of array
peak_index = np.argmax(pulse_array)
shift = int(nBins / 2) - peak_index
pulse_array = np.roll(pulse_array, shift)



f = ROOT.TFile("output.root")




# def get_landau_shape(mpv=15.0, width=0.1, time_max=32.0, nBins=1024):
#     landau_fn = ROOT.TF1("landau", f"TMath::Landau(x, {mpv}, {width})", 0, time_max)
#     h_landau = ROOT.TH1F("h_landau", "Landau Pulse", nBins, 0, time_max)
    
#     for i in range(1, nBins + 1):
#         x = h_landau.GetBinCenter(i)
#         h_landau.SetBinContent(i, landau_fn.Eval(x))

#     #return for nonnormalized shape
#     return np.array([h_landau.GetBinContent(i) for i in range(1, nBins + 1)])
    
#     # #Normalize (optional, depending on your deconvolution assumptions)
#     # values = np.array([h_landau.GetBinContent(i) for i in range(1, nBins + 1)])
#     # return values / np.sum(values)  # or omit normalization if you want raw shape



# Richardson-Lucy deconvolution
def richardson_lucy(observed, psf, iterations=30, clip=True):
    psf = psf / np.sum(psf)
    psf_mirror = psf[::-1]                                       #reverses the list (in this case psf); makes a mirror of the psf, needed for deconvolution
    
    estimate = np.full_like(observed, 0.5)                       #creates new array with the same shape/data type as observed, filled with 0.5
    
    for _ in range(iterations):
        conv = fftconvolve(estimate, psf, mode='same')
        conv = np.maximum(conv, 1e-8)
        ratio = observed / conv
        estimate *= fftconvolve(ratio, psf_mirror, mode='same')
        if clip:
            estimate = np.clip(estimate, 0, 1e6)

    return estimate


truth_dict = OrderedDict()

for key in f.GetListOfKeys():
    tname = key.GetName()
    if not tname.startswith("h_truth"):
        continue

    tmatch = re.match(r"h_truth_evt(\d+)_x(\d+)_y(\d+)", tname)

    if not tmatch:
        continue

    tievt, tx, ty = tmatch.groups()
    tx = int(tx)
    ty = int(ty)
    tievt = int(tievt)

    thist = f.Get(tname)
    if not thist:
        continue

    truth_dict[(tievt, tx, ty)] = np.array([thist.GetBinContent(i + 1) for i in range(nBins)])


for key in f.GetListOfKeys():
    name = key.GetName()

    if not name.startswith("h_reco"):
        continue

    match = re.match(r"h_reco_evt(\d+)_x(\d+)_y(\d+)", name)

    if not match:
        continue

    ievt, x, y = match.groups()
    x = int(x)
    y = int(y)
    ievt = int(ievt)

    # if ievt >= nevts:
    #     continue

    hist = f.Get(name)
    if not hist:
        continue

    signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
    # smoothed = gaussian_filter1d(signal, sigma=5.0)


    # Extract signal and time axis
    # signal = np.array([hist.GetBinContent(i + 1) for i in range(nBins)])
    time_axis = np.linspace(time_start + time_per_bin / 2, time_max - time_per_bin / 2, nBins)
    # shifted_time_axis = time_axis + truth_shift


    # Build response and apply deconvolution
    # response = pulses
    # response = pulses[:nBins]  # just to be safe if lengths differ
    # response = get_landau_shape(width=0.4)
    response = pulse_array


    # Step 1: light smoothing to suppress tiny spikes
    signal_smoothed = gaussian_filter1d(signal, sigma=2)

    # Step 2: threshold out sub-threshold noise
    noise_level = np.std(signal_smoothed[:100])  # adjust region as needed
    thresholded_signal = np.where(signal_smoothed > noise_level * 3, signal_smoothed, 0)

    # Step 3: Deconvolve
    deconv = richardson_lucy(thresholded_signal, response)


    # deconv = richardson_lucy(signal, response)
    # deconv = richardson_lucy(smoothed, response)


    # Find peaks in the deconvolved signal
    peaks, properties = find_peaks(deconv, height=10)

    h_delta = ROOT.TH1F(f"h_delta_evt{ievt}_x{x}_y{y}", "Deconvolved Peaks", nBins, time_start, time_max)

    for p in peaks:
        t = time_axis[p]  # Convert peak index to time
        h_delta.Fill(t, deconv[p])

    delta = np.array([h_delta.GetBinContent(i + 1) for i in range(nBins)])
    

    
    truthhits = truth_dict.get((ievt, x, y), np.zeros(nBins))
    scalefactor = 50

    nonzero_indices = np.nonzero(truthhits)[0]
    if len(nonzero_indices) > 0:
        first_hit_time = time_axis[nonzero_indices[0]]
        shift = (time_max / 2) - first_hit_time  
    else:
        shift = 0.0

    shifted_time_axis = time_axis - shift  



    # Plot
    plt.figure()
    # plt.plot(shifted_time_axis, smoothed, label="Observed Smoothed", alpha=0.5)
    plt.plot(shifted_time_axis, signal, label="Observed", alpha=0.5)
    plt.plot(shifted_time_axis, response * np.max(signal), '--', label="Response (scaled)")
    # plt.plot(shifted_time_axis, response, '--', label="Response (scaled)")
    plt.plot(shifted_time_axis, deconv, label="Deconvolved", linewidth=1)
    # plt.plot(time_axis, truthhits * scalefactor, label="Truth Hits")
    plt.plot(shifted_time_axis, delta, label="Reconstructed Truth Hits")
    plt.xlabel("Time Slice")
    plt.ylabel("Amplitude [arb.]")
    plt.xlim(0, time_max)
    plt.title(f"Evt {ievt} | x{x} y{y}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"deconv/deconv_evt{ievt}_iX{x}_iY{y}.png")
    plt.close()