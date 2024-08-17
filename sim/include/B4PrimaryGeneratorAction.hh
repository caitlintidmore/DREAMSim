//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4PrimaryGeneratorAction.hh
/// \brief Definition of the B4PrimaryGeneratorAction class

#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ParticleGun.hh"
#include "B4DetectorConstruction.hh"

#include <stdlib.h> /* getenv */

class TFile;
class Py8Jet;

class G4ParticleGun;
class G4Event;

class CaloTree;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class B4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  B4PrimaryGeneratorAction(B4DetectorConstruction *det, CaloTree *histo);
  virtual ~B4PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event *event);

  // set methods
  void SetRandomFlag(G4bool value);
  G4ParticleGun *GetParticleGun() { return fParticleGun; };

  void FillHEPparticles(std::vector<int> *mHepPID, std::vector<int> *mHepStatus,
                        std::vector<int> *mHepMother1, std::vector<int> *mHepMother2,
                        std::vector<int> *mHepDaughter1, std::vector<int> *mHepDaughter2,
                        std::vector<float> *mHepPx, std::vector<float> *mHepPy,
                        std::vector<float> *mHepPz, std::vector<float> *mHepE,
                        std::vector<float> *mHepMass);

private:
  G4ParticleGun *fParticleGun; // G4 particle gun
  B4DetectorConstruction *fDetector;
  CaloTree *hh;

  G4ParticleTable *particleTable;

  double worldZHalfLength;
  void getPy8Event(G4Event *event);
  void printPy8Event();

  // parameters from env. variables
  void getParamFromEnvVars();
  int CaloXPythiaON;                       //  0=singl particle gun, 1=Pythia8 Root file.
  double CaloXPythiaXmin, CaloXPythiaXmax; //  vertex point smearing.
  double CaloXPythiaYmin, CaloXPythiaYmax;
  double CaloXPythiaZmin, CaloXPythiaZmax;
  int CaloXPythiaSkip;  //  number of events to skip
  int CaloXPythiaPrint; //  number of events to print
  std::string CaloXPythiaFile;

  TFile *finPy8;
  Py8Jet *py8evt;
  int py8eventCounter;
  int py8eventNumber;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
