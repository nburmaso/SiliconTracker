#include "STSteppingAction.hh"
#include "STDetectorConstruction.hh"
#include "STEventAction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

STSteppingAction::STSteppingAction(STEventAction* eventAction)
  : G4UserSteppingAction(),
    fEventAction(eventAction),
    fScoringVolume(nullptr),
    fTrackID(),
    fELoss(),
    fTime(),
    fPosIn(),
    fPosOut()
{
}

STSteppingAction::~STSteppingAction()
{
}

void STSteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { // need to do it once
    G4RunManager* gRunManager = G4RunManager::GetRunManager();
    const auto* detectorConstruction = (const STDetectorConstruction*)gRunManager->GetUserDetectorConstruction();
    fScoringVolume = detectorConstruction->GetScoringVolume();
    fEventAction->setNSiLayers(detectorConstruction->getNSiLayers());
    fEventAction->setLayersDist(detectorConstruction->getLayersDistance());
    fEventAction->setLayersThic(detectorConstruction->getLayersThic());
  }

  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (volume != fScoringVolume) {
    return;
  }

  G4StepPoint* pointIn = step->GetPreStepPoint();
  G4StepPoint* pointOut = step->GetPostStepPoint();
  G4Track* track = step->GetTrack();
  //  G4ParticleDefinition* part = track->GetDefinition();
  //  int pdg = part->GetPDGEncoding();
  G4StepStatus statusIn = pointIn->GetStepStatus();
  G4StepStatus statusOut = pointOut->GetStepStatus();
  G4TrackStatus statusTrack = track->GetTrackStatus();

  if (track->GetParentID() == 0) { // fixme : considering only primary particles for now
    if (statusIn == fGeomBoundary || statusIn == fUndefined) {
      fELoss = 0;
      fTime = pointIn->GetGlobalTime();
      fPosIn = pointIn->GetPosition();
    }

    fELoss += step->GetTotalEnergyDeposit();

    if (statusOut == fGeomBoundary || statusTrack != fAlive) {
      fTrackID = track->GetTrackID();
      fPosOut = pointOut->GetPosition();
      fEventAction->AddMcPoint(STMcPoint(fTrackID, -1, fELoss, fTime, fPosIn, fPosOut));
    }

    double p = std::sqrt(track->GetMomentum().getX()*track->GetMomentum().getX() +
                          track->GetMomentum().getY()*track->GetMomentum().getY() +
                          track->GetMomentum().getZ()*track->GetMomentum().getZ());
    printf("LOG(INFO): STSteppingAction::UserSteppingAction(): track momentum = %.3f\n", p);
  }
}
