import ROOT
import re

# res_values = [0.1, 0.2, 0.3, 0.4, 0.5]
# noise_values = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
# xmax = 3
# ymax = 4
# ievt = 0  # assuming just the first event

# f = ROOT.TFile("output.root")
# actual_names = [k.GetName() for k in f.GetListOfKeys()]

# missing = []

# for x in range(xmax):
#     for y in range(ymax):
#         for res in res_values:
#             for noise in noise_values:
#                 hname = f"reco_evt{ievt}_x{x}_y{y}_res{res}_noise{noise}"
#                 if hname not in actual_names:
#                     missing.append(hname)

# print(f"Missing {len(missing)} histograms:")
# for name in missing:
#     print(name)


f = ROOT.TFile("output.root")
print("Keys in file:")
for key in f.GetListOfKeys():
    print(key.GetName())