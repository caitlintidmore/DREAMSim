### HG-DREAM G4 simulation: dream 2.04

#### Environment:

(Requirement: pyroot and GEANT4 are installed.)

For machines with cvmfs mounted, can directly source the environment
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-dbg/setup.sh
```
Change the path according to the OS.

On HPCC, everything is compiled inside the singularity environment, run
```
singularity run /lustre/research/hep/yofeng/SimulationEnv/alma9forgeant4_sbox
```
The corresponding docker image can be found [here](https://hub.docker.com/repository/docker/yongbinfeng/alma9geant/general), with the build file [here](https://github.com/TTU-HEP/SimulationEnv).


#### Compile:

build program in "build" area,
```
cd sim
mkdir build
cd build
cmake ..
make -j 4
```

Structure of software:

- `sim/exampleB4b.cc`: main program
- `sim/src/B4DetectorConstruction.cc`:definition of the detector
- `sim/src/B4bSteppingAction.cc`:access hits at each step
- `sim/src/CaloTree.cc`:analysis and hit handling

#### Run the code

```
source runBatch03_single_param.sh
```

The output files are:
- root: histograms
- csv: hits in each readout cell (2D and 3D)

Run parameters: 
- All run paramteres are defined in "paramBatch03_single.mac" and
may be overloaded in "runBatch03_single_param.sh".


#### Analysis

Simple analysis code (jupyter notebook), hgdream-3d-01.ipynb is under `plotter` subdirectory.

