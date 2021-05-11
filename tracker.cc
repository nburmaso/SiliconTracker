#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "include/STActionInitialization.hh"
#include "include/STDetectorConstruction.hh"

int main(int argc, char** argv)
{
  auto* ui = new G4UIExecutive(argc, argv);

  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed = time(nullptr);
  G4Random::setTheSeed(seed);
  CLHEP::HepRandom::setTheSeed(seed);

  auto* runManager = new G4RunManager;
  runManager->SetUserInitialization(new QBBC()); // physics list
  runManager->SetUserInitialization(new STActionInitialization());
  runManager->SetUserInitialization(new STDetectorConstruction());

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/execute run.mac");
//  ui->SessionStart();

//  delete ui;
  delete visManager;
  delete runManager;
}
