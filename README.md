### HG-DREAM G4 simulation: dream 2.04

Enivronment:

(Requirement: pyroot and GEANT4 are installed.)

ToDO: Finish the installation guide.

Then modify `muonSetup.sh` in dream2.03 to set environment for pyroot and geant4. Run
```
source muonSetupMac.sh
```

Compile code:

build program in  "build" area,

```
mkdir build01
cd build01
cmake ..
make
```

Modify code and run.

after a modification of code, re-build and run.

```
make
source runBatch03_single_param.sh
```

The output files are:
- root: histograms
- csv: hits in each readout cell (2D and 3D)

Run parameters:  

All run paramteres are defined in "paramBatch03_single.mac" and
may be overloaded in "runBatch03_single_param.sh".

Structure of software:

`sim/exampleB4b.cc`: main program
`sim/src/B4DetectorConstruction.cc`:definition of the detector
`sim/src/B4bSteppingAction.cc`:access hits at each step
`sim/src/CaloTree.cc`:analysis and hit handling


Simple analysis code (jupyter notebook), hgdream-3d-01.ipynb is under `plotter` subdirectory.

