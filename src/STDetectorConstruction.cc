#include "STDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

STDetectorConstruction::STDetectorConstruction()
  : G4VUserDetectorConstruction(),
    fScoringVolume(0)
{
}

STDetectorConstruction::~STDetectorConstruction()
{
}

G4VPhysicalVolume* STDetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* si_mat = nist->FindOrBuildMaterial("G4_Si");

  printf("\nLOG(INFO): detector material parameters: Z = %g, A = %g, Density = %g\n",
         si_mat->GetZ(), si_mat->GetA(), si_mat->GetDensity());

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* pb_mat = nist->FindOrBuildMaterial("G4_Pb");

  auto* solidWorld = new G4Box("World", 60 * cm, 60 * cm, 120 * cm); // half-dimensions
  auto* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,
                                                   "World", 0, false, 0, 0);

  // half-thickness in [mm]; 4.7 -> ~10% rad. length for Si
  fLayersThic = 2.35; // 5%
  // fLayersThic = 4.7;   // 10%
  // fLayersThic = 7.05;  // 15%
  // fLayersThic = 9.4;   // 20%
  auto* solidSi = new G4Box("SiLayer", 50 * cm, 50 * cm, fLayersThic * mm); // half-dimensions
  auto* logicSi = new G4LogicalVolume(solidSi, si_mat, "SiLayer");

  // lead foil for conversions
  /*
  auto *absorber = new G4Box("absorber", 40 * cm, 40 * cm, 1.5 * mm); // half-dimensions in mm
  auto *absorberLogic = new G4LogicalVolume(absorber, pb_mat, "absorber");
  G4VPhysicalVolume *absorberPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, 5 * cm),
                                                      absorberLogic, "absorber", logicWorld, false, 0, 0);
  */

  fNSiLayers = 10;
  fDistCM = 10; // [cm]
  for (uint i = 0; i < fNSiLayers; i++) {
    new G4PVPlacement(0, G4ThreeVector(0., 0., fDistCM * (i + 1) * cm), logicSi,
                      "SiLayer" + std::to_string(i), logicWorld, false, 0, 0);
  }

  // setting uniform magnetic field to the world
  auto* magField = new G4UniformMagField(G4ThreeVector(fMagField[0], fMagField[1] * tesla, fMagField[2]));
  auto* fieldMgr = new G4FieldManager();
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);
  logicWorld->SetFieldManager(fieldMgr, true);

  fScoringVolume = logicSi;
  return physWorld;
}
