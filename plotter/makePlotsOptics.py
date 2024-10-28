import ROOT
from collections import OrderedDict
import sys
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

print("Starting")

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(4)

elefile = "inputs/optics/electrons.txt"
pionfile = "inputs/optics/pions.txt"

loops = [('ele', elefile), ('pion', pionfile)]
loops = [('ele', elefile)]

rdfs = OrderedDict()
chains = OrderedDict()

for part, filename in loops:
    chains[part] = ROOT.TChain("tree")
    nfiles = 0
    with open(filename) as f:
        print(f"Reading {filename}")
        elefiles = f.readlines()
        for line in elefiles:
            line = line.strip()
            if line.startswith("#"):
                continue
            nfiles += 1
            print(f"{part} " + line)
            chains[part].Add(line)

            if nfiles > 100:
                break
    rdfs[part] = ROOT.RDataFrame(chains[part])

nEvts = OrderedDict()
for part, _ in loops:
    nEvts[part] = rdfs[part].Count().GetValue()
    print(f"Number of events for {part}: ", nEvts[part])

    rdfs[part] = rdfs[part].Define("eWeight", " 1.0 / " + str(nEvts[part]))
    rdfs[part] = rdfs[part].Define("OP_time_delta", "OP_time_final - OP_time_produced") # in ns
    rdfs[part] = rdfs[part].Define("OP_pos_delta_z", "OP_pos_final_z - OP_pos_produced_z")  # in cm
    rdfs[part] = rdfs[part].Define("OP_time_per_meter", "OP_time_delta / OP_pos_delta_z * 100.0") # in ns/m

histos = OrderedDict()
figures = ["nOPs", "OP_time_produced", "OP_time_final", "OP_time_delta", "OP_pos_delta_z", "OP_time_per_meter",
           "OP_pos_produced_x_vs_y", "OP_pos_final_x_vs_y", "OP_mom_produced_x_vs_y", "OP_mom_final_x_vs_y"]
# for event displays
evtlist = [1, 3, 5, 10, 15]
evtlist = []
for i in evtlist:
    figures.append(f"event_{i}_OP_pos_produced_x_vs_y")
    figures.append(f"event_{i}_OP_pos_final_x_vs_y")
    figures.append(f"event_{i}_OP_mom_produced_x_vs_y")
    figures.append(f"event_{i}_OP_mom_final_x_vs_y")

for fig in figures:
    histos[fig] = OrderedDict()

x_range = 0.15
px_range = 4e-9
t_range = 20

for part, rdf in rdfs.items():
    suffix = "_" + part

    histos['nOPs'][part] = rdf.Histo1D(
        ("nOPs" + suffix, "nOPs", 100, 0, 10000), "nOPs")
    histos['OP_time_produced'][part] = rdf.Histo1D(
        ("OP_time_produced" + suffix, "OP_time_produced", 100, 0, t_range), "OP_time_produced", "eWeight")
    histos['OP_time_final'][part] = rdf.Histo1D(
        ("OP_time_final" + suffix, "OP_time_final", 100, 0, t_range), "OP_time_final", "eWeight")
    histos['OP_time_delta'][part] = rdf.Histo1D(
        ("OP_time_delta" + suffix, "OP_time_delta", 100, 0, t_range), "OP_time_delta", "eWeight")
    histos["OP_pos_delta_z"][part] = rdf.Histo1D(
        ("OP_pos_delta_z" + suffix, "OP_pos_delta_z", 100, 0, 100.0), "OP_pos_delta_z", "eWeight")
    histos["OP_time_per_meter"][part] = rdf.Histo1D(
        ("OP_time_per_meter" + suffix, "OP_time_per_meter", 100, 0, 12.0), "OP_time_per_meter", "eWeight")

    histos["OP_pos_produced_x_vs_y"][part] = rdf.Histo2D(
        ("OP_pos_produced_x_vs_y" + suffix, "OP_pos_produced_x_vs_y", 50, -x_range, x_range, 50, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "eWeight")
    histos["OP_pos_final_x_vs_y"][part] = rdf.Histo2D(
        ("OP_pos_final_x_vs_y" + suffix, "OP_pos_final_x_vs_y", 50, -x_range, x_range, 50, -x_range, x_range), "OP_pos_final_x", "OP_pos_final_y", "eWeight")
    histos["OP_mom_produced_x_vs_y"][part] = rdf.Histo2D(("OP_mom_produced_x_vs_y" + suffix, "OP_mom_produced_x_vs_y",
                                                         50, -px_range, px_range, 50, -px_range, px_range), "OP_mom_produced_x", "OP_mom_produced_y", "eWeight")
    histos["OP_mom_final_x_vs_y"][part] = rdf.Histo2D(("OP_mom_final_x_vs_y" + suffix, "OP_mom_final_x_vs_y",
                                                      50, -px_range, px_range, 50, -px_range, px_range), "OP_mom_final_x", "OP_mom_final_y", "eWeight")

    # some event displays
    for i in evtlist:
        rdf_event = rdf.Filter(f"rdfentry_ == {i}")
        histos[f"event_{i}_OP_pos_produced_x_vs_y"][part] = rdf_event.Histo2D((f"event_{i}_produced_x_vs_y" + suffix,
                                                                               f"event_{i}_OP_pos_produced_x_vs_y", 50, -x_range, x_range, 50, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y")
        histos[f"event_{i}_OP_pos_final_x_vs_y"][part] = rdf_event.Histo2D((f"event_{i}_final_x_vs_y" + suffix,
                                                                            f"event_{i}_OP_pos_final_x_vs_y", 50, -x_range, x_range, 50, -x_range, x_range), "OP_pos_final_x", "OP_pos_final_y")
        histos[f"event_{i}_OP_mom_produced_x_vs_y"][part] = rdf_event.Histo2D((f"event_{i}_produced_px_vs_py" + suffix,
                                                                               f"event_{i}_OP_mom_produced_x_vs_y", 50, -px_range, px_range, 50, -px_range, px_range), "OP_mom_produced_x", "OP_mom_produced_y")
        histos[f"event_{i}_OP_mom_final_x_vs_y"][part] = rdf_event.Histo2D((f"event_{i}_final_px_vs_py" + suffix,
                                                                            f"event_{i}_OP_final_x_vs_y", 50, -px_range, px_range, 50, -px_range, px_range), "OP_mom_final_x", "OP_mom_final_y")

colormaps = {
    'ele': 2,
    'pion': 3
}


def GetColors(ene_fracs):
    colors = []
    for str in ene_fracs.keys():
        part, _ = str.split("_")
        color = colormaps[part]
        colors.append(color)
    return colors


args = {
    'dology': True,
    'mycolors': [colormaps[part] for part in rdfs.keys()],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True,
    'donormalize': False
}

print("Drawing")

DrawHistos(list(histos['nOPs'].values()), list(histos['nOPs'].keys(
)), 0, 10000, "Number of OPs", 1e-1, 1e4, "Fraction of events", "nOPs", **args)
DrawHistos(list(histos['OP_time_produced'].values()), list(histos['OP_time_produced'].keys(
)), 0, t_range, "Time [ns]", 1e-1, 1e7, "Fraction of events", "OP_time_produced", **args)
DrawHistos(list(histos['OP_time_final'].values()), list(histos['OP_time_final'].keys(
)), 0, t_range, "Time [ns]", 1e-1, 1e7, "Fraction of events", "OP_time_final", **args)
DrawHistos(list(histos['OP_time_delta'].values()), list(histos['OP_time_delta'].keys(
)), 0, t_range, "Time [ns]", 1e-1, 1e7, "Fraction of events", "OP_time_delta", **args)
DrawHistos(list(histos['OP_pos_delta_z'].values()), list(histos['OP_pos_delta_z'].keys(
)), 0, 100, "z [cm]", 1e-1, 1e7, "Fraction of events", "OP_pos_delta_z", **args)
DrawHistos(list(histos['OP_time_per_meter'].values()), list(histos['OP_time_per_meter'].keys(
)), 0, 12, "Time [ns/m]", 1e-1, 1e7, "Fraction of events", "OP_time_per_meter", **args)


# 2D plots
args['dology'] = False
args['drawoptions'] = "colz"
args['dologz'] = True
args['zmax'] = 1e3
args['zmin'] = 1e0
args['doth2'] = True
args['addOverflow'] = False
args['addUnderflow'] = False
for part in rdfs.keys():
    DrawHistos([histos['OP_pos_produced_x_vs_y'][part]], [], -x_range, x_range,
               "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_produced_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_pos_final_x_vs_y'][part]], [], -x_range, x_range,
               "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_final_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_mom_produced_x_vs_y'][part]], [], -px_range, px_range,
               "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_produced_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_mom_final_x_vs_y'][part]], [], -px_range, px_range,
               "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_final_x_vs_y_{part}", **args)

    # event displays
    for i in evtlist:
        DrawHistos([histos[f"event_{i}_OP_pos_produced_x_vs_y"][part]], [], -x_range, x_range,
                   "x [cm]", -x_range, x_range, "y [cm]", f"event_{i}_OP_pos_produced_x_vs_y_{part}", **args)

        DrawHistos([histos[f"event_{i}_OP_pos_final_x_vs_y"][part]], [], -x_range, x_range,
                   "x [cm]", -x_range, x_range, "y [cm]", f"event_{i}_OP_pos_final_x_vs_y_{part}", **args)
        DrawHistos([histos[f"event_{i}_OP_mom_produced_x_vs_y"][part]], [], -px_range, px_range,
                   "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"event_{i}_OP_mom_produced_x_vs_y_{part}", **args)
        DrawHistos([histos[f"event_{i}_OP_mom_final_x_vs_y"][part]], [], -px_range, px_range,
                   "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"event_{i}_OP_final_x_vs_y_{part}", **args)

print("Done")
