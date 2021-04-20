#ifndef S1SteppingAction_h
#define S1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"

class STEventAction;

class STSteppingAction : public G4UserSteppingAction
{
 public:
  STSteppingAction(STEventAction* eventAction);
  virtual ~STSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

 private:
  STEventAction* fEventAction;
  G4LogicalVolume* fScoringVolume;
  int fTrackID;
  double fELoss;
  double fTime;
  G4ThreeVector fPosIn;
  G4ThreeVector fPosOut;
};

#endif
