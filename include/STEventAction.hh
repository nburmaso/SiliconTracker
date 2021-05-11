#ifndef S1EventAction_h
#define S1EventAction_h 1

#include "G4UserEventAction.hh"
#include "STDigi.hh"
#include "STDigiMatch.hh"
#include "STHit.hh"
#include "STMcPoint.hh"
#include "STTrack.hh"

#include "TGraph.h"

#include <vector>

class TH2D;
class TH1D;
class STRunAction;

class STEventAction : public G4UserEventAction
{
 public:
  STEventAction(STRunAction* runAction);

  virtual ~STEventAction();

  virtual void BeginOfEventAction(const G4Event* event);

  virtual void EndOfEventAction(const G4Event* event);

  void AddMcPoint(const STMcPoint& point) { vMcPoints.push_back(point); }

  void setNSiLayers(uint nLayers) { fNSiLayers = nLayers; }
  void setLayersDist(double distCM) { fDistCM = distCM; }
  void setLayersThic(double thick) { fLayersThic = thick; }
  void drawQA(STTrack& preFitTrack,
              std::vector<double>& coordsZ,
              std::vector<std::vector<double>>& recoStates,
              TGraph* trajectory,
              const std::string& filename);

 private:
  void Digitizer();
  void HitProducer();
  void TrackFinderIdeal();
  void TrackFitter();
  void AddNeighbours(STHit& hit, std::vector<bool>& digiUsed, int iDigi);

  // call KF fitter
  void fitTracksKF();

  STRunAction* fRunAction;
  std::vector<STMcPoint> vMcPoints;
  std::vector<STDigi> vDigis;
  std::vector<STDigiMatch> vDigiMatches;
  std::vector<STHit> vHits;
  std::vector<STTrack> vTracks;

  double fLayerSizeX;
  double fLayerSizeY;
  double fLayerCenterX;
  double fLayerCenterY;
  TH2D* fHistDigis4;
  TH1D* fHistoP;

  // fixme: clean up?
  uint fNSiLayers;
  double fDistCM;
  double fLayersThic;

  bool fDebug{false};
  bool fPulls{false};
  bool fDrawTraj{false};
};

#endif
