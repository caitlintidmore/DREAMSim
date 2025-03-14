import ROOT
from collections import OrderedDict

p_data = "/Users/yfeng/Desktop/TTU/CaloX/Simulations/DREAMSim/DREAMSim/plotter/data/Sensl_FastOut_AveragePulse_1p8GHzBandwidth.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
h_pulse = pfile.Get("hist2")

sim_data = "/Users/yfeng/Desktop/mc_testjob_run001_003_Test_20evt_pi+_100.0_100.0.root"
ifile = ROOT.TFile(sim_data)
tree = ifile.Get("tree")

nFibers = 4
time_max = 30.0
time_per_bin = 0.04  # use 40 ps per bin
nBins = int(time_max / time_per_bin)

# pulses for truth photons and reco with shapes
histos_truth = OrderedDict()
histos_reco = OrderedDict()
for i in range(nFibers):
    histos_truth[i] = ROOT.TH1F(
        f"h_truth_C_{i}", f"h_truth_C_{i}", nBins, 0, time_max)
    histos_reco[i] = ROOT.TH1F(
        f"h_reco_C_{i}", f"h_reco_C_{i}", nBins, 0, time_max)


nevts = tree.GetEntries()
for i in range(nevts):
    tree.GetEntry(i)

    if i != 5:
        continue

    nPhotons = tree.nOPs
    for j in range(nPhotons):

        # make it to the end of the fiber
        if not tree.OP_pos_final_z[j] > 50.0:
            continue
        # look at photons from core of Cherenkov cone
        if not tree.OP_isCoreC[j]:
            continue

        if not tree.OP_productionFiber[j] == 0:
            continue

        t = tree.OP_time_final[j]
        itime_bin = int(t / time_per_bin)

        histos_truth[0].Fill(t, 1.0)

        for ipulse_bin in range()
