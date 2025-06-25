import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import re

# Output directory
os.makedirs("hists", exist_ok=True)

# Parameters
time_max = 30.0              # Duration in ns
time_per_bin = 0.04          # 40 ps per bin
nBins = int(time_max / time_per_bin)  # Should give 750
assert nBins == 750

# # Get ROOT file and histogram
# f = ROOT.TFile("output.root")
# h_pulse = f.Get("h_pulse")  # Replace with your actual histogram name if different

# # Extract original binning from ROOT histogram
# nbins_orig = h_pulse.GetNbinsX()
# x_vals = np.array([h_pulse.GetBinCenter(i+1) for i in range(nbins_orig)])
# y_vals = np.array([h_pulse.GetBinContent(i+1) for i in range(nbins_orig)])

# # Interpolate pulse to uniform binning over 750 bins in 0–30 ns
# interp_func = interp1d(x_vals, y_vals, kind='cubic', bounds_error=False, fill_value=0.0)
# x_interp = np.linspace(0, time_max, nBins)
# y_interp = interp_func(x_interp)

# # Optional: truncate last few bins to eliminate artificial dip at the end
# # Example: zero out values after 29.5 ns
# cutoff_time = 29.5
# cutoff_bin = int(cutoff_time / time_per_bin)
# y_interp[cutoff_bin:] = 0

# # Plot for verification
# plt.plot(x_interp, y_interp, label="Interpolated Pulse")
# plt.xlabel("Time (ns)")
# plt.ylabel("Amplitude")
# plt.title("Interpolated Pulse Shape (750 bins)")
# plt.legend()
# plt.grid(True)
# plt.savefig("hists/interpolated_pulse.png")
# plt.close()

# # Save as ROOT histogram (if needed for deconvolution or later analysis)
# h_interp = ROOT.TH1F("h_interp", "Interpolated Pulse", nBins, 0, time_max)
# for i in range(nBins):
#     h_interp.SetBinContent(i+1, y_interp[i])

# # Optional: save histogram to new ROOT file
# outf = ROOT.TFile("hists/pulse_interp_output.root", "RECREATE")
# h_interp.Write()
# outf.Close()













#original code before chatGPT messed it up
import sys
from collections import OrderedDict
import ROOT
import matplotlib.pyplot as plt
import numpy as np
import os
import re

os.makedirs("hists", exist_ok=True)
onechannel = OrderedDict()
nevt = np.array([957, 121763, 167417, 170989, 184764, 202407, 214496])

single_pulse_shapes = OrderedDict()
f = ROOT.TFile("filtered_events_board1_pulse_shapes.root")
for key in f.GetListOfKeys():
    name = key.GetName()
    print(name)
    for evt in nevt:
        if f"pulse_shape_Sci_Evt{evt}" in name and "iX1_iY1" in name:
            single_pulse_shapes[name] = f.Get(name)
            
            
ofile = ROOT.TFile("data_pulses.root", "RECREATE")
for name, hist in single_pulse_shapes.items():
    match = re.match(r"pulse_shape_Sci_Evt(\d+)_iX(\d+)_iY(\d+)", name)
    if match:
        evt, x, y = match.groups()
        new_name = f"pulse_evt{evt}"
        hist.SetName(new_name)
    hist.Write()
ofile.Close()
            

#         hist = f.Get(name)
#         if not hist:
#             continue
        
#         max_pulse = hist.GetMaximum()
#         if max_pulse < 10:
#             continue

        
#         print(hist.GetName())

#         histarray = np.zeros(hist.GetNbinsX())
#         for i in range(hist.GetNbinsX()):
#             histarray[i] = hist.GetBinContent(i + 1)

#         plt.plot(histarray)
#         plt.xlabel("Time [ns]")
#         plt.ylabel("Voltage [mV]")
#         plt.title(f"Pulse Shape - {name}")
#         plt.savefig(f"hists/{name}.png")
#         plt.close()

# #pulse_shape_Sci_Evt214496_iX1_iY1