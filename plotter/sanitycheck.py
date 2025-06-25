import ROOT
import re
import matplotlib.pyplot as plt
import numpy as np

# Open your ROOT file again
f = ROOT.TFile("output.root")
keys = f.GetListOfKeys()

# Dict to hold bin usage counts
bin_map = {}

# Loop through all histogram keys
for key in keys:
    name = key.GetName()
    if name.startswith("h_reco_evt"):  # or "h_truth_evt"
        match = re.search(r"_x(-?\d+)_y(-?\d+)", name)
        if match:
            xbin, ybin = int(match.group(1)), int(match.group(2))
            bin_map[(xbin, ybin)] = bin_map.get((xbin, ybin), 0) + 1

# Build 2D array from bin_map
xvals = [x for x, _ in bin_map]
yvals = [y for _, y in bin_map]
xrange = max(xvals) - min(xvals) + 1
yrange = max(yvals) - min(yvals) + 1

# Shift indices so they start at 0
xshift = min(xvals)
yshift = min(yvals)

heatmap = np.zeros((yrange, xrange))
for (x, y), count in bin_map.items():
    heatmap[y - yshift, x - xshift] = count  # row = y, col = x

# Plot
plt.figure(figsize=(8, 6))
plt.imshow(heatmap, origin='lower', cmap='viridis', interpolation='nearest')
plt.colorbar(label='Number of Events or Pulses')
plt.xlabel("x-bin")
plt.ylabel("y-bin")
plt.title("Filled Spatial Bins (x, y)")

plt.savefig("bin_check.png")
print("Saved plot to bin_check.png")