import ROOT
from collections import OrderedDict
import sys
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

print("Starting")

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(6)

elefile = "inputs/electrons.txt"
pionfile = "inputs/pions.txt"

chains = OrderedDict()
chains['ele'] = ROOT.TChain("tree")
chains['pion'] = ROOT.TChain("tree")

loops = [('ele', elefile), ('pion', pionfile)]

rdfs = OrderedDict()

for part, filename in loops:
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
            
rdfs['ele'] = rdfs['ele'].Define("eTotaltruth", "eLeaktruth + eCalotruth + eWorldtruth")
rdfs['pion'] = rdfs['pion'].Define("eTotaltruth", "eLeaktruth + eCalotruth + eWorldtruth")

rdfs['ele'] = rdfs['ele'].Define("rtruth", "sqrt(xtruth*xtruth + ytruth*ytruth)")
rdfs['pion'] = rdfs['pion'].Define("rtruth", "sqrt(xtruth*xtruth + ytruth*ytruth)")

nEvts = OrderedDict()
nEvts['ele'] = rdfs['ele'].Count().GetValue()
nEvts['pion'] = rdfs['pion'].Count().GetValue()
print("Number of events for electrons: ", nEvts['ele'])
print("Number of events for pions: ", nEvts['pion'])

rdfs['ele'] = rdfs['ele'].Define("eweight", f"edeptruth/ {nEvts['ele']}")
rdfs['pion'] = rdfs['pion'].Define("eweight", f"edeptruth/ {nEvts['pion']}")

histos = OrderedDict()
figures = ['eLeaktruth', 'eCalotruth', 'eTotaltruth',
           'eRodtruth', 'eCentruth', 'eScintruth',
           'xtruth', 'ytruth', 'ztruth', 'rtruth', 
           'time', 'time_zoomed',
           'xtruth_vs_ytruth', 'xtruth_vs_ztruth', 'rtruth_vs_ztruth',
           'time_vs_ztruth', 'time_vs_rtruth']
for fig in figures:
    histos[fig] = OrderedDict()

for part, rdf in rdfs.items():
    suffix = "_" + part
    
    histos['eLeaktruth'][part] = rdf.Histo1D(("eLeaktruth" + suffix, "eLeaktruth", 50, 0, 20.0), "eLeaktruth")
    histos['eCalotruth'][part] = rdf.Histo1D(("eCalotruth" + suffix, "eCalotruth", 50, 0., 102.0), "eCalotruth")
    histos['eTotaltruth'][part] = rdf.Histo1D(("eTotaltruth" + suffix, "eTotaltruth", 50, 70., 102.0), "eTotaltruth")
    histos['eRodtruth'][part] = rdf.Histo1D(("eRodtruth" + suffix, "eRodtruth", 50, 50, 102.0), "eRodtruth")
    histos['eCentruth'][part] = rdf.Histo1D(("eCentruth" + suffix, "eCentruth", 50, 1.0, 5.0), "eCentruth")
    histos['eScintruth'][part] = rdf.Histo1D(("eScintruth" + suffix, "eScintruth", 50, 1.0, 5.0), "eScintruth")
    
    
    # energy weighted
    histos['xtruth'][part] = rdf.Histo1D(("xtruth" + suffix, "xtruth", 50, -20, 20), "xtruth", "eweight")
    histos['ytruth'][part] = rdf.Histo1D(("ytruth" + suffix, "ytruth", 50, -20, 20), "ytruth", "eweight")
    histos['ztruth'][part] = rdf.Histo1D(("ztruth" + suffix, "ztruth", 100, -100, 100), "ztruth", "eweight")
    histos['rtruth'][part] = rdf.Histo1D(("rtruth" + suffix, "rtruth", 50, 0, 30), "rtruth", "eweight")
    
    histos['time'][part] = rdf.Histo1D(("time" + suffix, "time", 50, 0, 100), "globaltimetruth", "eweight")
    histos['time_zoomed'][part] = rdf.Histo1D(("time_zoomed" + suffix, "time_zoomed", 50, 0, 10), "globaltimetruth", "eweight")
    
    histos['xtruth_vs_ytruth'][part] = rdf.Histo2D(("xtruth_vs_ytruth" + suffix, "xtruth_vs_ytruth", 50, -20, 20, 50, -20, 20), "xtruth", "ytruth", "eweight")
    histos['xtruth_vs_ztruth'][part] = rdf.Histo2D(("xtruth_vs_ztruth" + suffix, "xtruth_vs_ztruth", 50, -20, 20, 100, -100, 100), "xtruth", "ztruth", "eweight")
    histos['rtruth_vs_ztruth'][part] = rdf.Histo2D(("rtruth_vs_ztruth" + suffix, "rtruth_vs_ztruth", 50, 0, 30, 100, -100, 100), "rtruth", "ztruth", "eweight")
    
    histos['time_vs_ztruth'][part] = rdf.Histo2D(("time_vs_ztruth" + suffix, "time_vs_ztruth", 50, 0, 100, 100, -100, 100), "globaltimetruth", "ztruth", "eweight")
    histos['time_vs_rtruth'][part] = rdf.Histo2D(("time_vs_rtruth" + suffix, "time_vs_rtruth", 50, 0, 100, 50, 0, 30), "globaltimetruth", "rtruth", "eweight")
    
    
colormaps ={
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
    'donormalize': True,
    'mycolors': [colormaps[part] for part in rdfs.keys()],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True
}

print("Drawing")

DrawHistos(list(histos['eLeaktruth'].values()), list(histos['eLeaktruth'].keys()), 0, 20, "Leakage Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eLeaktruth", **args)
DrawHistos(list(histos['eCalotruth'].values()), list(histos['eCalotruth'].keys()), 0., 102, "Calo Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCalotruth", **args)
DrawHistos(list(histos['eTotaltruth'].values()), list(histos['eTotaltruth'].keys()), 70., 102, "Total Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eTotaltruth", **args)
DrawHistos(list(histos['eRodtruth'].values()), list(histos['eRodtruth'].keys()), 50, 102, "Rod Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eRodtruth", **args)
DrawHistos(list(histos['eCentruth'].values()), list(histos['eCentruth'].keys()), 1, 5, "CFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCentruth", **args)
DrawHistos(list(histos['eScintruth'].values()), list(histos['eScintruth'].keys()), 1, 5, "SFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eScintruth", **args)

args['donormalize'] = False
DrawHistos(list(histos['xtruth'].values()), list(histos['xtruth'].keys()), -20, 20, "x [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "xtruth", **args)
DrawHistos(list(histos['ytruth'].values()), list(histos['ytruth'].keys()), -20, 20, "y [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "ytruth", **args)
DrawHistos(list(histos['ztruth'].values()), list(histos['ztruth'].keys()), -100, 100, "z [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "ztruth", **args)
DrawHistos(list(histos['rtruth'].values()), list(histos['rtruth'].keys()), 0, 30, "r [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "rtruth", **args)

DrawHistos(list(histos['time'].values()), list(histos['time'].keys()), 0, 100, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time", **args)
DrawHistos(list(histos['time_zoomed'].values()), list(histos['time_zoomed'].keys()), 0, 10, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time_zoomed", **args)

# 2D plots
args['dology'] = False
args['drawoptions'] = "colz"
args['dologz'] = True
for part in rdfs.keys():
    DrawHistos([histos['xtruth_vs_ytruth'][part]], [], -20, 20, "x [cm]", -20, 20, "y [cm]", f"xtruth_vs_ytruth_{part}", **args)
    DrawHistos([histos['xtruth_vs_ztruth'][part]], [], -20, 20, "x [cm]", -100, 100, "z [cm]", f"xtruth_vs_ztruth_{part}", **args)
    DrawHistos([histos['rtruth_vs_ztruth'][part]], [], 0, 30, "r [cm]", -100, 100, "z [cm]", f"rtruth_vs_ztruth_{part}", **args)    
    
    DrawHistos([histos['time_vs_ztruth'][part]], [], 0, 100, "Time [ns]", -100, 100, "z [cm]", f"time_vs_ztruth_{part}", **args)
    DrawHistos([histos['time_vs_rtruth'][part]], [], 0, 100, "Time [ns]", 0, 30, "r [cm]", f"time_vs_rtruth_{part}", **args)

print("Done")
