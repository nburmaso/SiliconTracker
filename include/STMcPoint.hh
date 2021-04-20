#ifndef S1McPoint_h
#define S1McPoint_h 1

#include "G4ThreeVector.hh"

class STMcPoint
{
 public:
  STMcPoint(int trackID, int fVolumeIndex, double eLoss, double t, G4ThreeVector posIn, G4ThreeVector posOut);
  virtual ~STMcPoint() {}

  int GetTrackID() { return fTrackID; }
  int GetVolumeIndex() { return fVolumeIndex; }
  double GetELoss() { return fELoss; }
  double GetTime() { return fTime; }
  G4ThreeVector GetPosIn() { return fPosIn; }
  G4ThreeVector GetPosOut() { return fPosOut; }
  void Print();

 private:
  int fTrackID;
  int fVolumeIndex;
  double fELoss;
  double fTime;
  G4ThreeVector fPosIn;
  G4ThreeVector fPosOut;
};

#endif
