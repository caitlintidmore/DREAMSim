#!/bin/bash

conda activate my_root_env

cd /Users/kunori/skdir/hep/g4/geant4-v11.2.2-install/share/Geant4/geant4make
source geant4make.sh
cd -
export G4BIN="$PWD"


