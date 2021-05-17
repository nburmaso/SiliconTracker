#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "STDetectorConstruction.hh"
#include "STEventAction.hh"
#include "STRunAction.hh"
#include "STSteppingAction.hh"

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
    auto* gRunManager = G4RunManager::GetRunManager();
    gRunManager->GetCurrentEvent()->GetEventID();
    const auto* detectorConstruction = (const STDetectorConstruction*)gRunManager->GetUserDetectorConstruction();
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (volume != fScoringVolume) {
    return;
  }

  G4StepPoint* pointIn = step->GetPreStepPoint();
  G4StepPoint* pointOut = step->GetPostStepPoint();
  G4Track* track = step->GetTrack();
  G4StepStatus statusIn = pointIn->GetStepStatus();
  G4StepStatus statusOut = pointOut->GetStepStatus();
  G4TrackStatus statusTrack = track->GetTrackStatus();

  if (track->GetParentID() == 0 || track->GetParentID() == 1) { // considering only primary and secondary particles for now
    if (statusIn == fGeomBoundary || statusIn == fUndefined) {
      fELoss = 0;
      fTime = pointIn->GetGlobalTime();
      fPosIn = pointIn->GetPosition();
    }

    if (step->IsFirstStepInVolume()) {
      if (!fEventAction->isTrackStored(track->GetTrackID())) {
        auto* gRunManager = G4RunManager::GetRunManager();
        int eventID = gRunManager->GetCurrentEvent()->GetEventID();
        double kine = track->GetVertexKineticEnergy();
        G4ParticleDefinition* part = track->GetDefinition();
        double mass = part->GetPDGMass();
        double energy = kine + mass;
        printf("Mass = %f\n", part->GetPDGMass());
        double p = std::sqrt(energy * energy - mass * mass);
        int pdg = part->GetPDGEncoding();
        fEventAction->getRunAction()->mcTrack.mcEventID = eventID;
        fEventAction->getRunAction()->mcTrack.mcTrackID = track->GetTrackID();
        fEventAction->getRunAction()->mcTrack.pdgID = pdg;
        fEventAction->getRunAction()->mcTrack.motherID = track->GetParentID();
        fEventAction->getRunAction()->mcTrack.vx = track->GetVertexPosition().getX();
        fEventAction->getRunAction()->mcTrack.vy = track->GetVertexPosition().getY();
        fEventAction->getRunAction()->mcTrack.vz = track->GetVertexPosition().getZ();
        fEventAction->getRunAction()->mcTrack.px = track->GetVertexMomentumDirection().getX() * p;
        fEventAction->getRunAction()->mcTrack.py = track->GetVertexMomentumDirection().getY() * p;
        fEventAction->getRunAction()->mcTrack.pz = track->GetVertexMomentumDirection().getZ() * p;
        fEventAction->getRunAction()->tMCTracks->Fill();
        fEventAction->addTrackID(track->GetTrackID());
      }
    }

    fELoss += step->GetTotalEnergyDeposit();

    if (statusOut == fGeomBoundary || statusTrack != fAlive) {
      fTrackID = track->GetTrackID();
      fPosOut = pointOut->GetPosition();
      fEventAction->AddMcPoint(STMcPoint(fTrackID, -1, fELoss, fTime, fPosIn, fPosOut));
    }
  }
}
