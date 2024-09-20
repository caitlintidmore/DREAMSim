import ROOT
from collections import OrderedDict
import sys
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

ROOT.gROOT.SetBatch(True)

all_files = OrderedDict()
all_files['ele'] = OrderedDict()
all_files['pion'] = OrderedDict()

all_files['pion']['10GeV'] = "/afs/cern.ch/work/y/yofeng/public/CaloX/outputs/mc_testjob_run002_003_Test_100evt_pi+_10.0_10.0.root"
all_files['pion']['100GeV'] = "/afs/cern.ch/work/y/yofeng/public/CaloX/outputs/mc_testjob_run002_003_Test_100evt_pi+_100.0_100.0.root"
all_files['ele']['10GeV'] = "/afs/cern.ch/work/y/yofeng/public/CaloX/outputs/mc_testjob_run001_003_Test_100evt_e+_10.0_10.0.root"
all_files['ele']['100GeV'] = "/afs/cern.ch/work/y/yofeng/public/CaloX/outputs/mc_testjob_run001_003_Test_100evt_e+_100.0_100.0.root"

particles = ['ele', 'pion']
energies = ['10GeV', '100GeV']

files = OrderedDict()
trees = OrderedDict()
for part in particles:
    trees[part] = OrderedDict()
    files[part] = OrderedDict()
    for ene in energies:
        files[part][ene] = ROOT.TFile.Open(all_files[part][ene])
        trees[part][ene] = files[part][ene].Get("tree").Clone()

ene_fracs = OrderedDict()
eneC_fracs = OrderedDict()
eneS_fracs = OrderedDict()

hists_z = OrderedDict()
hists_z_C = OrderedDict()
hists_z_S = OrderedDict()

hists_x = OrderedDict()
hists_x_C = OrderedDict()
hists_x_S = OrderedDict()

hists_z_weighted = OrderedDict()
hists_z_C_weighted = OrderedDict()
hists_z_S_weighted = OrderedDict()

particles = ['ele', 'pion']
energies = ['10GeV']

for part in particles:
    for ene in energies:
        str = part + "_" + ene
        
        e = float(ene.replace("GeV", ""))
        
        hname = "energyfrac_" + str
        ene_fracs[str] = ROOT.TH1F(hname, hname, 100, 0, 0.1)
        eneC_fracs[str] = ROOT.TH1F(hname + "_C", hname + "_C", 100, 0, 0.1)
        eneS_fracs[str] = ROOT.TH1F(hname + "_S", hname + "_S", 100, 0, 0.1)
        
        print("file: " + all_files[part][ene])
        print("Filling " + str)
        trees[part][ene].Draw(f"Sum$(edeptruth*(calotypetruth==2 || calotypetruth==3)) / {e} >> {ene_fracs[str].GetName()}")
        trees[part][ene].Draw(f"Sum$(edeptruth*(calotypetruth==3)) / {e} >> {eneC_fracs[str].GetName()}")
        trees[part][ene].Draw(f"Sum$(edeptruth*(calotypetruth==2)) / {e} >> {eneS_fracs[str].GetName()}")
        
        hname = "zdistance_" + str
        hists_z[str] = ROOT.TH1F(hname + "_z", hname + "_z", 100, -100., 100.0)
        hists_z_C[str] = ROOT.TH1F(hname + "_C_z", hname + "_C_z", 100, -100., 100.0)
        hists_z_S[str] = ROOT.TH1F(hname + "_S_z", hname + "_S_z", 100, -100., 100.0)
        
        trees[part][ene].Draw("ztruth*(calotypetruth==2 || calotypetruth==3)>> " + hists_z[str].GetName())
        trees[part][ene].Draw("ztruth*(calotypetruth==3)>> " + hists_z_C[str].GetName())
        trees[part][ene].Draw("ztruth*(calotypetruth==2)>> " + hists_z_S[str].GetName())
        
        hname = "xdistance_" + str
        hists_x[str] = ROOT.TH1F(hname + "_x", hname + "_x", 100, -15., 15.0)
        hists_x_C[str] = ROOT.TH1F(hname + "_C_x", hname + "_C_x", 100, -15., 15.0)
        hists_x_S[str] = ROOT.TH1F(hname + "_S_x", hname + "_S_x", 100, -15., 15.0)
        
        trees[part][ene].Draw("xtruth*(calotypetruth==2 || calotypetruth==3)>> " + hists_x[str].GetName())
        trees[part][ene].Draw("xtruth*(calotypetruth==3)>> " + hists_x_C[str].GetName())
        trees[part][ene].Draw("xtruth*(calotypetruth==2)>> " + hists_x_S[str].GetName())
        
        hname = "zdistance_weighted_" + str
        hists_z_weighted[str] = ROOT.TH1F(hname + "_z_weighted", hname + "_z_weighted", 100, -100., 100.0)
        hists_z_C_weighted[str] = ROOT.TH1F(hname + "_C_z_weighted", hname + "_C_z_weighted", 100, -100., 100.0)
        hists_z_S_weighted[str] = ROOT.TH1F(hname + "_S_z_weighted", hname + "_S_z_weighted", 100, -100., 100.0)
        
        trees[part][ene].Draw("ztruth*(calotypetruth==2 || calotypetruth==3)>> " + hists_z_weighted[str].GetName(), "edeptruth")
        trees[part][ene].Draw("ztruth*(calotypetruth==3)>> " + hists_z_C_weighted[str].GetName(), "edeptruth")
        trees[part][ene].Draw("ztruth*(calotypetruth==2)>> " + hists_z_S_weighted[str].GetName(), "edeptruth")
        
        
colormaps ={
    'ele': 2,
    'pion': 3
}
linestylemaps = {
    '10GeV': 1,
    '100GeV': 2
}

def GetColors(ene_fracs):
    colors = []
    for str in ene_fracs.keys():
        part, _ = str.split("_")
        color = colormaps[part]
        colors.append(color)
    return colors

def GetLineStyles(ene_fracs):
    lineStyles = []
    for str in ene_fracs.keys():
        _, ene = str.split("_")
        lineStyle = linestylemaps[ene]
        lineStyles.append(lineStyle)
    return lineStyles
        
args = {
    'dology': True,
    'donormalize': True,
    'mycolors': GetColors(ene_fracs),
    'linestyles': GetLineStyles(ene_fracs),
    "MCOnly": True,
}

DrawHistos(list(ene_fracs.values()), list(ene_fracs.keys()), 0, 0.08, "Fraction of energy", 1e-3, 1e2, "Fraction of events", "ene_frac", **args)

DrawHistos(list(eneC_fracs.values()), list(eneC_fracs.keys()), 0, 0.05, "Fraction of energy", 1e-3, 1e2, "Fraction of events", "eneC_frac", **args)

DrawHistos(list(eneS_fracs.values()), list(eneS_fracs.keys()), 0, 0.05, "Fraction of energy", 1e-3, 1e2, "Fraction of events", "eneS_frac", **args)

DrawHistos(list(hists_z.values()), list(hists_z.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z", **args)
DrawHistos(list(hists_z_C.values()), list(hists_z_C.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z_C", **args)
DrawHistos(list(hists_z_S.values()), list(hists_z_S.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z_S", **args)

DrawHistos(list(hists_x.values()), list(hists_x.keys()), -15, 15, "x [cm]", 1e-6, 1e4, "Fraction of events", "x", **args)
DrawHistos(list(hists_x_C.values()), list(hists_x_C.keys()), -15, 15, "x [cm]", 1e-6, 1e4, "Fraction of events", "x_C", **args)
DrawHistos(list(hists_x_S.values()), list(hists_x_S.keys()), -15, 15, "x [cm]", 1e-6, 1e4, "Fraction of events", "x_S", **args)


args['donormalize'] = False
DrawHistos(list(hists_z_weighted.values()), list(hists_z_weighted.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z_weighted", **args)
DrawHistos(list(hists_z_C_weighted.values()), list(hists_z_C_weighted.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z_C_weighted", **args)
DrawHistos(list(hists_z_S_weighted.values()), list(hists_z_S_weighted.keys()), -100, 100, "z [cm]", 1e-6, 1e4, "Fraction of events", "z_S_weighted", **args)
