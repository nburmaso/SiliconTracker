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

  void setNSiLayers(uint nLayers) { fNSiLayers = nLayers; }
  uint getNSiLayers() const { return fNSiLayers; }
  void setLayersDistance(double distCM) { fDistCM = distCM; }
  double getLayersDistance() const { return fDistCM; }
  double getLayersThic() const { return fLayersThic; }

 protected:
  G4LogicalVolume* fScoringVolume;

 private:
  uint fNSiLayers{5};
  double fDistCM{20.};
  double fLayersThic{10.};
};

#endif
