#with spike filtering
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os

# Create output directory
os.makedirs("drsdata", exist_ok=True)
threshold = 200

# Load ROOT file and tree
f = ROOT.TFile.Open("run1040_250705180247_converted.root")
tree = f.Get("EventTree")

# Define the branch name you want to access
branch_name = "DRS_Board1_Group0_Channel1"

# Loop over all entries and store branch content as NumPy arrays
arrays = []
for i, event in enumerate(tree):

    data = getattr(event, branch_name)

    if max(data) < threshold:
        continue

    np_data = np.array(data, dtype=np.float32)
    baseline = np.median(np_data)

    base_data = np_data - baseline  # Subtract baseline

    if max(base_data) < threshold:
        continue

    arrays.append(base_data)

    # if len(arrays) >= 50:
    #     break

# for i, base_data in enumerate(arrays):    #enumerate give the index and the value at index
#     peak = np.max(base_data)
#     plt.figure()
#     plt.plot(base_data, label=f"Entry {i} (peak = {peak:.2f})")
#     plt.xlabel("Channel index")
#     plt.ylabel("Energy")
#     plt.xlim(0, 1000)
#     plt.title(f"{branch_name} - Event {i}")
#     plt.legend()
#     plt.savefig(f"drsdata/evt{i}.png")
#     plt.close()

# save to file
ofile = ROOT.TFile("drsoutput.root", "RECREATE")

for i, base_data in enumerate(arrays):
    hist = ROOT.TH1F(f"waveform_{i}", f"Filtered waveform {i}", len(base_data), 0, len(base_data))
    for j, val in enumerate(base_data):
        hist.SetBinContent(j + 1, val)
    hist.Write()

ofile.Close()
