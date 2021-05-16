#ifndef S1DetectorConstruction_h
#define S1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class STDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  STDetectorConstruction();
  virtual ~STDetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

 protected:
  G4LogicalVolume* fScoringVolume;
};

#endif
