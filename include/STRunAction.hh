#ifndef S1RunAction_h
#define S1RunAction_h 1

#include "G4UserRunAction.hh"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

class G4Run;
class TList;

class STRunAction : public G4UserRunAction
{
 public:
  STRunAction();
  virtual ~STRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  TList* GetList() { return fList; }

  // residuals/pulls histograms
  TH1F* hXResiduals;
  TH1F* hYResiduals;
  TH1F* hTxResiduals;
  TH1F* hTyResiduals;
  TH1F* hPResiduals;
  TH1F* hXPulls;
  TH1F* hYPulls;
  TH1F* hTxPulls;
  TH1F* hTyPulls;
  TH1F* hPPulls;
  TH1F* hPSigmas;
  TH1F* hChi2s;

  TFile* fTracks;
  TTree* tTracks;

  TFile* fMCTracks;
  TTree* tMCTracks;

  struct
  {
    int eventID;
    int trackID;
    int mcLabel;
    double z;
    double x;
    double y;
    double tx;
    double ty;
    double qp;
    double covM[15];
  } recoTrack;

  struct
  {
    int mcEventID;
    int mcTrackID;
    int pdgID;
    int motherID;
    double vx; // vertex coordinates
    double vy; // vertex coordinates
    double vz; // vertex coordinates
    double px;
    double py;
    double pz;
  } mcTrack;

 private:
  TList* fList;
};

#endif
