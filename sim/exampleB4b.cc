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
/// \file exampleB4b.cc
/// \brief Main program of the B4b example

// #include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
// #include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
// #include "G4Cerenkov.hh"
#include "Randomize.hh"

#include "G4HadronicProcess.hh"
#include "G4GammaGeneralProcess.hh"

#include "B4DetectorConstruction.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "B4bRunAction.hh"
#include "B4bEventAction.hh"
#include "B4bSteppingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "CaloTree.h"
#include "B4bPhysicsList.hh"

#
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{

  bool batchJob = false;

  string macro;

  for (G4int i = 1; i < argc; i = i + 2)
  {
    string a = argv[i];
    if (G4String(argv[i]) == "-b")
    {
      macro = argv[i + 1];
      batchJob = true;
    }
    else if (G4String(argv[i]) == "-i")
    {
      macro = argv[i + 1];
      batchJob = false;
    }
    else if (a.substr(0, 1) != "-")
    {
      std::cout << "argument error: parameter shoudl start with -. " << a << std::endl;
      return 0;
    }
  }

  CaloTree *histo = new CaloTree(macro, argc, argv);

  G4UIExecutive *ui = nullptr;
  if (!batchJob)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Generate rndom number seeds;
  long seeds[2];
  // std::chrono::system_clock::time_point now=std::chrono::system_clock::now();
  // auto duration = now.time_since_epoch();
  // auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
  // long long  micros = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
  // std::cout<<"millis "<<millis <<std::endl;
  // std::cout<<"micros "<<micros <<std::endl;
  // long long t1=micros/10000000;
  // long long t2=micros-t1*10000000;
  long long t1 = 1234;
  long long t2 = 456;
  std::cout << "t1=" << t1 << "   t2=" << t2 << std::endl;
  int kseed = histo->getParamI("runNumber") + histo->getParamI("runSeq") * 3333;
  seeds[0] = long(t2) + long(kseed);
  seeds[1] = seeds[0] + 8134;

  // seeds[0]=2345;
  // seeds[1]=7999;

  G4Random::setTheSeeds(seeds);
  G4Random::showEngineStatus();
  std::cout << "seeds[0]=" << seeds[0] << "   seeds[1]=" << seeds[1] << std::endl;

  // Construct a serial run manager
  //
  // auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  auto *runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  auto detector = new B4DetectorConstruction(histo);
  runManager->SetUserInitialization(detector);

  //   optical physics from examples/extended/optical/OpNovice2
  //  auto physicsList = new FTFP_BERT;
  auto physicsList = new QGSP_BERT;
  // auto physicsList = new PhysListEmStandard();
  //  auto physicsList = new CustomPhysicsList();

  // physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);

  runManager->SetUserInitialization(physicsList);

  // G4Cerenkov* theCerenkovProcess=new G4Cerenkov("Cerenkov");
  // theCerenkovProcess->SetTrackSecondariesFirst(true);
  // nt MaxNumberPhotons=300;
  // theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumberPhotons);
  // physicsList->RegisterPhysics(theCerenkovProcess);

  auto gen_action = new B4PrimaryGeneratorAction(detector, histo);
  runManager->SetUserAction(gen_action);

  auto run_action = new B4bRunAction(histo);
  runManager->SetUserAction(run_action);
  //
  auto event_action = new B4bEventAction(detector, gen_action, histo);
  runManager->SetUserAction(event_action);
  //
  auto stepping_action = new B4bSteppingAction(event_action, histo);
  runManager->SetUserAction(stepping_action);

  runManager->Initialize();
  // auto actionInitialization = new B4bActionInitialization(detConstruction);
  // runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  if (0)
  {
    // Get the photon (gamma) particle definition
    G4ParticleDefinition *gamma = G4Gamma::GammaDefinition();
    G4ProcessManager *pmanager = gamma->GetProcessManager();
    std::cout << "pmanager=" << pmanager << std::endl;
    if (pmanager)
    {
      G4int nprocesses = pmanager->GetProcessListLength();
      std::cout << "nprocesses=" << nprocesses << std::endl;
      for (G4int i = 0; i < nprocesses; i++)
      {
        G4VProcess *process = (*pmanager->GetProcessList())[i];
        std::cout << "process name: " << process->GetProcessName() << std::endl;

        if (process->GetProcessName() == "GammaGeneralProc")
        {
          if (dynamic_cast<G4GammaGeneralProcess *>(process))
          {
            // https://github.com/Geant4/geant4/blob/v11.2.2/source/physics_lists/constructors/electromagnetic/include/G4GammaGeneralProcess.hh#L78
            // and https://github.com/Geant4/geant4/blob/v11.2.2/source/physics_lists/constructors/electromagnetic/src/G4GammaGeneralProcess.cc#L619
            ((G4GammaGeneralProcess *)process)->AddHadProcess(nullptr);
            std::cout << "Set hadronic process of gamma to null" << std::endl;
          }
        }
      }
    }
    std::cout << "done" << std::endl;
  }

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if (batchJob)
  {
    // batch mode
    std::cout << "debug:  batch mode" << std::endl;
    G4String command = "/control/execute ";
    cout << "command: " << command << endl;
    UImanager->ApplyCommand(command + macro); // macro does not have /run/beamOn 100

    //    beams are now defined in B4PrimaryGeneratorAction...
    // command="/gun/particle "+histo->getParamS("gun_particle");
    // cout<<"command: "<<command<<endl;
    // UImanager->ApplyCommand(command);

    // command="/gun/energy "+histo->getParamS("gun_energy")+" GeV";
    // cout<<"command: "<<command<<endl;
    // UImanager->ApplyCommand(command);;

    // string evtmax="100";
    command = "/run/beamOn " + histo->getParamS("numberOfEvents");
    cout << "command: " << command << endl;
    UImanager->ApplyCommand(command);
  }
  else
  {
    // interactive mode
    std::cout << "debug:  interactive mode" << std::endl;
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    std::cout << "debug:  interactive mode, step 2" << std::endl;
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  std::cout << "Job termination..." << std::endl;
  histo->EndJob();
  std::cout << "after histo->EndJob()..." << std::endl;

  G4Random::showEngineStatus();
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete histo;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
