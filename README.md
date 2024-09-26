### HG-DREAM G4 simulation: dream 2.06

#### Environment:

For machines with cvmfs mounted, can directly source the environment
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-dbg/setup.sh
```
Change the path according to the OS.

**On HPCC**, everything (ROOT and GEANT4) is compiled inside the singularity environment. Log into the interactive node ([more information](https://www.depts.ttu.edu/hpcc/userguides/Job_User_Guide.pdf)) with e.g.

```
interactive -p nocona
```
From there run the singularity container with the following command:
```
singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_sbox/
```
The corresponding docker image can be found [here](https://hub.docker.com/repository/docker/yongbinfeng/alma9geant/general), with the build file [here](https://github.com/TTU-HEP/SimulationEnv).

**Note** if you have conda installed, exit the conda environment before running the singularity container, otherwise it might cause conflicts with different ROOT versions etc.


#### Compile:

Inside the singularity environment, build program in "build" area,
```
cd /path/to/DREAMSIM/directory
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
./exampleB4b -b paramBatch03_single.mac  \
    -jobName testjob -runNumber 001 -runSeq 003  \
    -numberOfEvents 10  -eventsInNtupe 100    \
    -gun_particle e+ -gun_energy_min 100.0 -gun_energy_max 100.0 \
    -sipmType 1
```

The output files are:
- root: histograms
- csv: hits in each readout cell (2D and 3D)

Run parameters: 
- All run parameters are defined in "paramBatch03_single.mac" and
may be overloaded in "runBatch03_single_param.sh".

#### Job submission on HPCC
the script `jobs/jobSubmission.py` handles that. Run
```
cd jobs
python jobSubmission.py
# It will produce the shell scripts to submit and run the jobs on HPCC.
bash submit_all.sh
```
to submit the jobs, which can then be monitored with `squeue -u $USER`. More information on the HPCC batch system can be found [here](https://www.depts.ttu.edu/hpcc/userguides/Job_User_Guide.pdf).

#### Analysis

The script `plotter/makePlots.py` handles this.