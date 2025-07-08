#with spike filtering
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os

# Create output directory
os.makedirs("drsdata", exist_ok=True)

# Load ROOT file and tree
f = ROOT.TFile.Open("run1033_250704164003.root")
tree = f.Get("EventTree")

# Define the branch name you want to access
branch_name = "DRS_Board1_Group0_Channel0"

# Loop over all entries and store branch content as NumPy arrays
arrays = []
for i, event in enumerate(tree):
    data = getattr(event, branch_name)
    np_data = np.array(data, dtype=np.float32)
    arrays.append(np_data)

    # Spike detection — compare max to median
    baseline = np.median(np_data)
    peak = np.max(np_data)
    ratio = peak / baseline if baseline > 0 else 0

    # Skip plotting if no significant spike (less than 25% above baseline)
    if ratio < 1.25:
        continue

    # Optional: plot the first few with spikes
    if i < 100:
        plt.figure()
        plt.plot(np_data, label=f"Entry {i} (ratio = {ratio:.2f})")
        plt.xlabel("Channel index")
        plt.ylabel("Energy")
        plt.xlim(0, 1000)
        plt.title(f"{branch_name} - Event {i}")
        plt.legend()
        plt.savefig(f"drsdata/evt{i}.png")
        plt.close()










# #without spike filtering
# import ROOT
# import numpy as np
# import matplotlib.pyplot as plt
# import os

# # Create output directory
# os.makedirs("drsdata", exist_ok=True)

# # Load ROOT file and tree
# f = ROOT.TFile.Open("run1033_250704164003.root")
# tree = f.Get("EventTree")

# # Define the branch name you want to access
# branch_name = "DRS_Board1_Group0_Channel0"

# # Loop over all entries and store branch content as NumPy arrays
# arrays = []
# for i, event in enumerate(tree):
#     # Get the branch content (usually a C++ array or vector)
#     data = getattr(event, branch_name)

#     # Convert to NumPy array
#     np_data = np.array(data, dtype=np.float32)  # Adjust dtype if needed
#     arrays.append(np_data)

#     # Optional: plot the first few entries
#     if i < 100:
#         plt.figure()
#         plt.plot(np_data, label=f"Entry {i}")
#         plt.xlabel("Channel index")
#         plt.ylabel("Energy")
#         plt.xlim(0, 1000)
#         plt.title(f"{branch_name} - Event {i}")
#         plt.legend()
#         plt.savefig(f"drsdata/evt{i}.png")
#         plt.close()
