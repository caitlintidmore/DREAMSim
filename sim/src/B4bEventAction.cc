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
/// \file B4bEventAction.cc
/// \brief Implementation of the B4bEventAction class

#include "B4bEventAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "G4Step.hh"

#include "B4PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include <iomanip>

// -- for CaloX data --
#include "CaloID.h"
#include "CaloHit.h"
#include "CaloTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::B4bEventAction(B4DetectorConstruction *det, B4PrimaryGeneratorAction *prim, CaloTree *histo)
    : G4UserEventAction(), fDetector(det), primary(prim), hh(histo)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bEventAction::~B4bEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::BeginOfEventAction(const G4Event *event /*event*/)
{
   // std::cout<<"B4bEventAction::BeginOfEventAction-  starting..."<<std::endl;
   // G4Random::showEngineStatus();
   hh->BeginEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bEventAction::EndOfEventAction(const G4Event *event)
{
   // std::cout<<"B4bEventAction::EndOfEventAction-  starting..."<<std::endl;
   hh->EndEvent();

} //  end of B4bEventAction::EndOfEventAction

// -----------------------------------------------------------------------
