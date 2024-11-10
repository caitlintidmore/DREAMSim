import ROOT
from collections import OrderedDict, defaultdict
import sys
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

doOP = True

print("Starting")

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(8)

elefile = "inputs/optics/electrons.txt"
pionfile = "inputs/optics/pions.txt"
opfile = "inputs/optics/ops.txt"

loops = [('ele', elefile), ('pion', pionfile)]
#loops = [('ele', elefile)]
if doOP:
    loops = [('op', opfile)]

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
    nEvts[part] = float(nEvts[part])
    print(f"Number of events for {part}: ", nEvts[part])
    
    if doOP:
        # for OPs, it is always 1 OP per event
        # no need to normalize
        nEvts[part] = 1.0

    rdfs[part] = rdfs[part].Define("OP_passEnd", "OP_pos_final_z > 49.0")
    rdfs[part] = rdfs[part].Define("eWeight", " OP_passEnd / " + str(nEvts[part]))
    rdfs[part] = rdfs[part].Define("eWeightTotal", "1.0 / " + str(nEvts[part]))
    rdfs[part] = rdfs[part].Define("OP_time_delta", "OP_time_final - OP_time_produced") # in ns
    rdfs[part] = rdfs[part].Define("OP_pos_delta_z", "OP_pos_final_z - OP_pos_produced_z")  # in cm
    rdfs[part] = rdfs[part].Define("OP_time_per_meter", "OP_time_delta / OP_pos_delta_z * 100.0") # in ns/m
    rdfs[part] = rdfs[part].Define("OP_mom_produced", "sqrt(OP_mom_produced_x*OP_mom_produced_x + OP_mom_produced_y*OP_mom_produced_y + OP_mom_produced_z*OP_mom_produced_z)")
    rdfs[part] = rdfs[part].Define("OP_sinTheta_produced", "OP_mom_produced_z / OP_mom_produced")
    rdfs[part] = rdfs[part].Define("OP_sinThetaInv_produced", "1.0 / OP_sinTheta_produced")
    rdfs[part] = rdfs[part].Define("OP_time_per_meter_sinTheta_produced", "OP_time_per_meter * OP_sinTheta_produced")

histos = defaultdict(OrderedDict)

# for event displays
evtlist = [1, 3, 5, 10, 15]
evtlist = []

x_range = 0.045
nx_bins = 100
px_range = 1e-8
if doOP:
    px_range = 1e-6
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
    histos["OP_sinTheta_produced"][part] = rdf.Histo1D(
        ("OP_sinTheta_produced" + suffix, "OP_sinTheta_produced", 100, 0, 1), "OP_sinTheta_produced", "eWeight")
    histos["OP_sinTheta_produced_total"][part] = rdf.Histo1D(
        ("OP_sinTheta_produced_total" + suffix, "OP_sinTheta_produced_total", 100, 0, 1), "OP_sinTheta_produced", "eWeightTotal")
    histos["OP_time_per_meter_sinTheta_produced"][part] = rdf.Histo1D(
        ("OP_time_per_meter_sinTheta_produced" + suffix, "OP_time_per_meter_sinTheta_produced", 100, 4.5, 6.5), "OP_time_per_meter_sinTheta_produced", "eWeight")

    histos["OP_pos_produced_x_vs_y"][part] = rdf.Histo2D(
        ("OP_pos_produced_x_vs_y" + suffix, "OP_pos_produced_x_vs_y", nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "eWeight")
    histos["OP_pos_produced_x_vs_y_total"][part] = rdf.Histo2D(
        ("OP_pos_produced_x_vs_y_total" + suffix, "OP_pos_produced_x_vs_y_total", nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "eWeightTotal")
    histos["OP_pos_final_x_vs_y"][part] = rdf.Histo2D(
        ("OP_pos_final_x_vs_y" + suffix, "OP_pos_final_x_vs_y", nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_final_x", "OP_pos_final_y", "eWeight")
    histos["OP_mom_produced_x_vs_y"][part] = rdf.Histo2D(("OP_mom_produced_x_vs_y" + suffix, "OP_mom_produced_x_vs_y",
                                                         nx_bins, -px_range, px_range, nx_bins, -px_range, px_range), "OP_mom_produced_x", "OP_mom_produced_y", "eWeight")
    histos["OP_mom_final_x_vs_y"][part] = rdf.Histo2D(("OP_mom_final_x_vs_y" + suffix, "OP_mom_final_x_vs_y",
                                                      nx_bins, -px_range, px_range, nx_bins, -px_range, px_range), "OP_mom_final_x", "OP_mom_final_y", "eWeight")
    histos["OP_mom_produced_x_vs_y_total"][part] = rdf.Histo2D(("OP_mom_produced_x_vs_y_total" + suffix, "OP_mom_produced_x_vs_y_total",
                                                                nx_bins, -px_range, px_range, nx_bins, -px_range, px_range), "OP_mom_produced_x", "OP_mom_produced_y", "eWeightTotal")
    
    histos["OP_time_per_meter_vs_sinTheta_produced"][part] = rdf.Histo2D(("OP_time_per_meter_vs_sinTheta_produced" + suffix, "OP_time_per_meter_vs_sinTheta_produced",
                                                         100, 0, 12, 100, 1, 3.0), "OP_time_per_meter", "OP_sinThetaInv_produced", "eWeight")
    
    histos["OP_profilePx_produced_x_vs_y"][part] = rdf.Profile2D(("OP_profilePx_produced_x_vs_y" + suffix, "OP_profilePx_produced_x_vs_y",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "OP_mom_produced_x", "eWeight")
    histos["OP_profilePy_produced_x_vs_y"][part] = rdf.Profile2D(("OP_profilePy_produced_x_vs_y" + suffix, "OP_profilePy_produced_x_vs_y",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "OP_mom_produced_y", "eWeight")
    histos["OP_profilePx_final_x_vs_y"][part] = rdf.Profile2D(("OP_profilePx_final_x_vs_y" + suffix, "OP_profilePx_final_x_vs_y",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_final_x", "OP_pos_final_y", "OP_mom_final_x", "eWeight")
    histos["OP_profilePy_final_x_vs_y"][part] = rdf.Profile2D(("OP_profilePy_final_x_vs_y" + suffix, "OP_profilePy_final_x_vs_y",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_final_x", "OP_pos_final_y", "OP_mom_final_y", "eWeight")
    
    histos["OP_profilePx_produced_x_vs_y_total"][part] = rdf.Profile2D(("OP_profilePx_produced_x_vs_y_total" + suffix, "OP_profilePx_produced_x_vs_y_total",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "OP_mom_produced_x", "eWeightTotal")
    histos["OP_profilePy_produced_x_vs_y_total"][part] = rdf.Profile2D(("OP_profilePy_produced_x_vs_y_total" + suffix, "OP_profilePy_produced_x_vs_y_total",nx_bins, -x_range, x_range, nx_bins, -x_range, x_range), "OP_pos_produced_x", "OP_pos_produced_y", "OP_mom_produced_y", "eWeightTotal")
    
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
    'pion': 3,
    "op": 4
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
DrawHistos(list(histos['OP_sinTheta_produced'].values()), list(histos['OP_sinTheta_produced'].keys(
)), 0, 1, "sin(#theta)", 1e-1, 1e7, "Fraction of events", "OP_sinTheta_produced", **args)
DrawHistos(list(histos['OP_sinTheta_produced_total'].values()), list(histos['OP_sinTheta_produced_total'].keys(
)), 0, 1, "sin(#theta)", 1e-1, 1e7, "Fraction of events", "OP_sinTheta_produced_total", **args)
DrawHistos(list(histos['OP_time_per_meter_sinTheta_produced'].values()), list(histos['OP_time_per_meter_sinTheta_produced'].keys(
)), 4.5, 6.5, "Time [ns/m]", 1e-1, 1e7, "Fraction of events", "OP_time_per_meter_sinTheta_produced", **args)

def GetRatio(histos, hname, parts=None):
    if hname not in histos:
        return None
    if hname + "_total" not in histos:
        return None
    
    if parts is None:
        parts = histos[hname].keys()
    
    h_ratios = []
    for part in parts:
        h = histos[hname][part].GetValue()
        hden = histos[hname + "_total"][part].GetValue()
    
        h_ratio = h.Clone(h.GetName() + "_ratio")
        h_ratio.Divide(hden)
        
        h_ratios.append(h_ratio)
    return h_ratios

h_ratios_sinTheta = GetRatio(histos, "OP_sinTheta_produced")
DrawHistos(h_ratios_sinTheta, [], 0, 1, "sin(#theta)", 1e-3, 10, "Fraction of events", "OP_sinTheta_produced_ratio", **args)

def makeArrowPlots(hprof2d_x, hprof2d_y, min_entries=1, min_value= 1e-10, scale = 1e7):
    # assumes hprof2d_x and hprof2d_y have same binning
    arrows = []
    nbinsX = hprof2d_x.GetNbinsX()
    nbinsY = hprof2d_x.GetNbinsY()
    for i in range(1, nbinsX+1):
        for j in range(1, nbinsY+1):
            bin_ij = hprof2d_x.GetBin(i, j)
            if hprof2d_x.GetBinEffectiveEntries(bin_ij) <= min_entries:
                continue
            avg_px = hprof2d_x.GetBinContent(i, j)
            avg_py = hprof2d_y.GetBinContent(i, j)
            avg_pxy = (avg_px**2 + avg_py**2)**0.5
            if avg_pxy <= min_value:
                continue
            
            x = hprof2d_x.GetXaxis().GetBinCenter(i)
            y = hprof2d_x.GetYaxis().GetBinCenter(j)
            dx = avg_px * scale
            dy = avg_py * scale
            arrow = ROOT.TArrow(x, y, x+dx, y+dy, 0.01, "|>")
            arrow.SetLineColor(ROOT.kPink)
            arrow.SetFillColor(ROOT.kPink)
            arrows.append(arrow)
    return arrows

# 2D plots
args['dology'] = False
args['drawoptions'] = "colz"
args['dologz'] = True
args['zmax'] = 1e2
args['zmin'] = 1e-4
if doOP:
    args['zmax'] = 1e3
    args['zmin'] = 1.0
args['doth2'] = True
args['addOverflow'] = False
args['addUnderflow'] = False

scale = 1e7
if doOP:
    scale = 2e4
    
for part in rdfs.keys():
    DrawHistos([histos['OP_pos_produced_x_vs_y'][part]], [], -x_range, x_range,
               "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_produced_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_pos_produced_x_vs_y_total'][part]], [], -x_range, x_range,
                "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_produced_x_vs_y_total_{part}", **args)
    hratio_x_vs_y = GetRatio(histos, "OP_pos_produced_x_vs_y", parts=[part])[0]
    DrawHistos([hratio_x_vs_y], [], -x_range, x_range, "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_produced_x_vs_y_ratio", **{**args, 'zmax': 1.0, 'zmin': 0.0})
    DrawHistos([histos['OP_pos_final_x_vs_y'][part]], [], -x_range, x_range,
               "x [cm]", -x_range, x_range, "y [cm]", f"OP_pos_final_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_mom_produced_x_vs_y'][part]], [], -px_range, px_range,
               "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_produced_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_mom_final_x_vs_y'][part]], [], -px_range, px_range,
               "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_final_x_vs_y_{part}", **args)
    DrawHistos([histos['OP_mom_produced_x_vs_y_total'][part]], [], -px_range, px_range,
                "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_produced_x_vs_y_total_{part}", **args)
    hratio_px_vs_py = GetRatio(histos, "OP_mom_produced_x_vs_y", parts=[part])[0]
    DrawHistos([hratio_px_vs_py], [], -px_range, px_range, "px [GeV/c]", -px_range, px_range, "py [GeV/c]", f"OP_mom_produced_x_vs_y_ratio", **{**args, 'zmax': 1.0, 'zmin': 0.0})
    
    DrawHistos([histos['OP_time_per_meter_vs_sinTheta_produced'][part]], [], 0, 12, "Time [ns/m]", 1, 2.5, "1.0/sin(#theta)", f"OP_time_per_meter_vs_sinTheta_produced_{part}", **args)
    
    # profile plots
    arrows = makeArrowPlots(histos['OP_profilePx_produced_x_vs_y'][part].GetValue(), histos['OP_profilePy_produced_x_vs_y'][part].GetValue(), scale=scale)
    DrawHistos([histos['OP_pos_produced_x_vs_y'][part]], [], -x_range, x_range,"x [cm]", -x_range, x_range, "y [cm]", f"OP_produced_x_vs_y_{part}_withPXY", **args, extraToDraw = arrows)
    
    arrows = makeArrowPlots(histos['OP_profilePx_final_x_vs_y'][part].GetValue(), histos['OP_profilePy_final_x_vs_y'][part].GetValue(), scale=scale)
    DrawHistos([histos['OP_pos_final_x_vs_y'][part]], [], -x_range, x_range,"x [cm]", -x_range, x_range, "y [cm]", f"OP_final_x_vs_y_{part}_withPXY", **args, extraToDraw = arrows)
    
    arrows = makeArrowPlots(histos['OP_profilePx_produced_x_vs_y_total'][part].GetValue(), histos['OP_profilePy_produced_x_vs_y_total'][part].GetValue(), scale=scale)
    DrawHistos([histos['OP_pos_produced_x_vs_y_total'][part]], [], -x_range, x_range,"x [cm]", -x_range, x_range, "y [cm]", f"OP_produced_x_vs_y_total_{part}_withPXY", **args, extraToDraw = arrows)
    
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
