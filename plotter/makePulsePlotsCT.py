import sys                                  #talk to os
from collections import OrderedDict         #for neatness
import ROOT
sys.path.append("./CMSPLOTS")               #puts cmsplots in path to access it anywhere
from myFunction import DrawHistos  # noqa

ROOT.gROOT.SetBatch(True)                   # Don't draw on screen, just save to file, makes code faster

# Draw single pulse
p_data = "data/Sensl_FastOut_AveragePulse_1p8GHzBandwidth.root"
pfile = ROOT.TFile(p_data)
# pulse shape file per photon (?)
h_pulse = pfile.Get("hist2")

args = {
    'dology': True,
    'donormalize': False,
    'mycolors': [1],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True,
    'mycolors': [2, 3, 4, 6, 7, 8, 9, 28, 30, 38, 46]
}
#Drawing pulse histogram from known shape 
DrawHistos([h_pulse], ["Single Pulse"], 0, 20, "Measured Time [ns]",
           1, 2.5e5, "Voltage [mV]", "SinglePulse_log", **args)
args['dology'] = False
DrawHistos([h_pulse], ["Single Pulse"], 0, 20, "Measured Time [ns]",
           -2e4, 2.5e5, "Voltage [mV]", "SinglePulse", **args)             

fname_pulses = "output.root"
f_pulses = ROOT.TFile(fname_pulses)

args['dology'] = False

# truth photon time
histos_truth = OrderedDict()
nevts = 1
for ievt in range(nevts):
    histos_truth[ievt] = OrderedDict()

    for hit in f_pulses.GetListOfKeys():
        if "truth" in hit.GetName() and f"evt{ievt}" in hit.GetName():
            histos_truth[ievt][hit.GetName()] = f_pulses.Get(hit.GetName())

    DrawHistos(list(histos_truth[ievt].values()), list(histos_truth[ievt].keys()), 0, 30, "Measured Time [ns]",
               0, 25, "Voltage [mV]", f"Truth_{ievt}", **args)

#LOG SCALE THING
args['dology'] = True

# reco photon time
histos_reco = OrderedDict()
for ievt in range(nevts):
    histos_reco[ievt] = OrderedDict()

    for hit in f_pulses.GetListOfKeys():
        if "reco" in hit.GetName() and f"evt{ievt}" in hit.GetName():
            histos_reco[ievt][hit.GetName()] = f_pulses.Get(hit.GetName())
            histos_reco[ievt][hit.GetName()].SetFillColor(0)

# Draw reco histograms
    DrawHistos(list(histos_reco[ievt].values()), list(histos_reco[ievt].keys()), 0, 30, "Measured Time [ns]",
               0.001, 1, "Voltage [mV]", f"Reco_{ievt}", **args)
