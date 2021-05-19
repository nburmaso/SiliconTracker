#include "STEventAction.hh"
#include "STRunAction.hh"
#include "STTrackFitter.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "STSimulationPars.hh"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TRandom.h"
#include "TRandomGen.h"
#include "TLegend.h"

#include <vector>
#include <map>

using namespace sim_pars;

STEventAction::STEventAction(STRunAction* runAction)
  : G4UserEventAction(),
    fRunAction(runAction),
    vMcPoints(),
    vDigis(),
    vDigiMatches(),
    vHits(),
    vTracks(),
    vTrackIDs(),
    fHistoP(nullptr),
    fHistDigis4(nullptr),
    // TODO read these parameters from geometry
    fLayerSizeX(1000.),
    fLayerSizeY(1000.),
    fLayerCenterX(0.),
    fLayerCenterY(0.)
{
}

STEventAction::~STEventAction()
{
  delete fHistDigis4;
  delete fHistoP;
}

void STEventAction::BeginOfEventAction(const G4Event* /* iEvent */)
{
  vMcPoints.clear();
  vDigis.clear();
  vDigiMatches.clear();
  vHits.clear();
  vTracks.clear();
  vTrackIDs.clear();
}

void STEventAction::EndOfEventAction(const G4Event* /* iEvent */)
{
  Digitizer();
  HitProducer();
  TrackFinderIdeal();
  //  TrackFitter();

  fitTracksKF();
}

void STEventAction::Digitizer()
{
  if (fDebug) {
    printf("LOG(INFO): STEventAction::Digitizer(): n layers = %d\n", eNSiLayers);
  }

  TList* list = fRunAction->GetList();
  fHistDigis4 = (TH2D*)list->FindObject("hDigis4");
  fHistoP = (TH1D*)list->FindObject("hP");

  double xmin = fLayerCenterX - fLayerSizeX / 2.;
  double ymin = fLayerCenterY - fLayerSizeY / 2.;
  double padSizeX = fLayerSizeX / gnRowsX; // 5 mm
  double padSizeY = fLayerSizeY / gnRowsY; // 10 mm

  // print MC points
  for (auto& vMcPoint : vMcPoints) {
    vMcPoint.Print();
  }

  /*
  std::vector<bool> isVolumeFired(eNSiLayers * 2, false);
  for (uint i = 0; i < vMcPoints.size(); i++) {
    STMcPoint& point = vMcPoints[i];
    int volumeIndex = point.GetVolumeIndex();
    isVolumeFired[volumeIndex] = true;
  }
*/

  // if (isVolumeFired[0] && isVolumeFired[4]) {
  //   histo->Fill(4);
  // }

  std::map<int, int> mapDigiIndex; // channel to index in digi array
  for (unsigned int i = 0; i < vMcPoints.size(); i++) {
    STMcPoint& point = vMcPoints[i];
    G4ThreeVector posIn = point.GetPosIn();
    G4ThreeVector posOut = point.GetPosOut();
    G4ThreeVector center = (posIn + posOut) / 2.;
    int iLayer = round(center.z() / (eLayersDist * 10)) - 1; // [cm] to [mm]
    printf("LOG(INFO): STEventAction::Digitizer(): iLayer = %d\n", iLayer);
    int iRowX = floor((center.x() - xmin) / padSizeX);
    int iRowY = floor((center.y() - ymin) / padSizeY);
    int iChannel = iLayer * gnRowsX * gnRowsY + iRowX * gnRowsY + iRowY;
    if (mapDigiIndex.find(iChannel) == mapDigiIndex.end()) {
      // add digi match for channel iChannel if it doesn't exist yet
      vDigis.push_back(STDigi(iChannel));
      vDigiMatches.push_back(STDigiMatch());
      mapDigiIndex[iChannel] = vDigiMatches.size() - 1;
    }
    vDigiMatches[mapDigiIndex[iChannel]].AddMcPointIndex(i);
  }

  for (unsigned int iDigi = 0; iDigi < vDigis.size(); iDigi++) { // loop over digis
    STDigi& digi = vDigis[iDigi];
    STDigiMatch& digiMatch = vDigiMatches[iDigi];
    double sumEloss = 0;
    double sumTimeEloss = 0; // needed to calculated eloss-weighted average
    for (int i = 0; i < digiMatch.GetNIndices(); i++) {
      // loop over MC points contributing to this digi
      Int_t iMcPoint = digiMatch.GetMcPointIndex(i);
      STMcPoint& point = vMcPoints[iMcPoint];
      sumEloss += point.GetELoss();
      sumTimeEloss += point.GetTime() * point.GetELoss();
    }

    // TODO store these parameters in a separate place (RunAction?)
    double elossMax = 2; // MeV
    int adcMax = 1024;   //  10-bit adc
    int adc = sumEloss > elossMax ? adcMax - 1 : floor(sumEloss / elossMax * adcMax);

    double ftime = sumTimeEloss / sumEloss; // simple eloss weighted average. Unrealistic but...
    int t = floor(ftime * 10);              // assuming 0.1 ns clock
    digi.SetAmplitude(adc);
    digi.SetTime(t);
  }

  for (unsigned int iDigi = 0; iDigi < vDigis.size(); iDigi++) { // loop over digis
    STDigi& digi = vDigis[iDigi];
    STDigiMatch& digiMatch = vDigiMatches[iDigi];
    printf("Digi %2d:", iDigi);
    printf(" ch=%6d", digi.GetChannel());
    printf(" time=%2d", digi.GetTime());
    printf(" adc=%4d", digi.GetAmplitude());
    printf(" contributors=%d\n", digiMatch.GetNIndices());
    if (digi.GetLayer() != eNSiLayers - 1) {
      continue;
    }
    fHistDigis4->SetBinContent(digi.GetRowX() + 1, digi.GetRowY() + 1, digi.GetAmplitude());
  }
}

void STEventAction::AddNeighbours(STHit& hit, std::vector<bool>& digiUsed, int iDigi)
{
  digiUsed[iDigi] = true;
  hit.AddDigiIndex(iDigi);
  STDigi& digi = vDigis[iDigi];
  int iLayer = digi.GetLayer();
  int iRowX = digi.GetRowX();
  int iRowY = digi.GetRowY();
  for (unsigned int iDigi2 = 0; iDigi2 < vDigis.size(); iDigi2++) {
    if (digiUsed[iDigi2]) {
      continue;
    }
    STDigi& digi2 = vDigis[iDigi2];
    if (digi2.GetLayer() != iLayer) {
      continue;
    }
    if (digi2.GetRowX() < iRowX - 1 || digi2.GetRowX() > iRowX + 1) {
      continue;
    }
    if (digi2.GetRowY() < iRowY - 1 || digi2.GetRowY() > iRowY + 1) {
      continue;
    }
    AddNeighbours(hit, digiUsed, iDigi2);
  }
}

void STEventAction::HitProducer()
{
  // Find clusters
  std::vector<bool> digiUsed; // flags to mark digis already attached to a cluster
  for (unsigned int iDigi = 0; iDigi < vDigis.size(); iDigi++) {
    digiUsed.push_back(0);
  }
  // Collecting neighbors into one cluster. Store cluster contributors in the hit
  for (unsigned int iDigi = 0; iDigi < vDigis.size(); iDigi++) { // loop over digis
    if (digiUsed[iDigi]) {
      continue;
    }
    STHit hit;
    AddNeighbours(hit, digiUsed, iDigi);
    vHits.push_back(hit);
  }

  double xmin = fLayerSizeX / 2. - fLayerSizeX;
  double ymin = fLayerSizeY / 2. - fLayerSizeY;
  double padSizeX = fLayerSizeX / gnRowsX;
  double padSizeY = fLayerSizeY / gnRowsY;
  // Setting hit coordinates as weighted average over contributing digis
  // (no cluster deconvolution)
  for (auto& hit : vHits) {
    double sumAdcX = 0;
    double sumAdcY = 0;
    double sumAdcZ = 0;
    double sumAdc = 0;
    int layer = -1;
    //  loop over digits attached to cluster
    for (int i = 0; i < hit.GetNDigiIndices(); i++) {
      int iDigi = hit.GetDigiIndex(i);
      STDigi& digi = vDigis[iDigi];
      layer = digi.GetLayer();
      double x = xmin + padSizeX * (digi.GetRowX() + 0.5f); // using pad centers
      double y = ymin + padSizeY * (digi.GetRowY() + 0.5f); // using pad centers
      double z = eLayersDist * 10 * (layer + 1);
      double adc = digi.GetAmplitude();
      sumAdcX += adc * x;
      sumAdcY += adc * y;
      sumAdcZ += adc * z;
      sumAdc += adc;
    }
    hit.SetLayer(layer);
    hit.SetAmplitude(sumAdc);
    hit.SetXYZ(sumAdcX / sumAdc, sumAdcY / sumAdc, sumAdcZ / sumAdc);

    // setting uncertainties of hit coordinates
    // assuming single-pad clusters
    // variance of the uniform distribution is (b-a)^2/12, see
    // https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    // using pad_width/sqrt(12)
    // can be improved for large clusters using cluster widths
    double dx = padSizeX / std::sqrt(12);
    double dy = padSizeY / std::sqrt(12);
    double dz = 10 / std::sqrt(12); // TODO Read from geometry
    hit.SetDxDyDz(dx, dy, dz);
  }

  for (unsigned int iHit = 0; iHit < vHits.size(); iHit++) {
    STHit& hit = vHits[iHit];
    printf("Hit %3d:", iHit);
    printf(" (x,y,z)=(%5.1f,%5.1f,%6.1f)", hit.X(), hit.Y(), hit.Z());
    printf(" (dx,dy,dz)=(%.2f,%.2f,%.2f)", hit.Dx(), hit.Dy(), hit.Dz());
    printf(" amplitude=%4d", hit.GetAmplitude());
    printf(" nDigits=%d\n", hit.GetNDigiIndices());
  }
}

void STEventAction::TrackFinderIdeal()
{
  // Hit <- Digi1   <- McPoint1 <- McTrack1
  //                <- McPoint2 <- McTrack2
  //        Digi2   <- McPoint3 <- McTrack3

  // Collect hits originating from the same MC track
  std::map<int, std::map<int, bool>> mapHitsToTrackID;
  for (unsigned int iHit = 0; iHit < vHits.size(); iHit++) { //hits
    STHit& hit = vHits[iHit];
    for (int iDigi = 0; iDigi < hit.GetNDigiIndices(); iDigi++) { //digis
      int iDigiIndex = hit.GetDigiIndex(iDigi);
      STDigiMatch& match = vDigiMatches[iDigiIndex];
      for (int iMcPoint = 0; iMcPoint < match.GetNIndices(); iMcPoint++) { // mc points
        int iMcPointIndex = match.GetMcPointIndex(iMcPoint);
        STMcPoint& point = vMcPoints[iMcPointIndex];
        int trackID = point.GetTrackID();
        mapHitsToTrackID[trackID][iHit] = true;
      } // mc points
    }   // digis
  }     // hits

  // Found tracks
  for (auto& itMap : mapHitsToTrackID) {
    std::map<int, bool>& mapTrackHits = itMap.second;
    // Skip tracks with less than 3 hits
    if (mapTrackHits.size() < 3) {
      continue;
    }

    // count number of layers with hits
    std::vector<int> hitIndex(eNSiLayers, -1);
    for (auto& itTrackHits : mapTrackHits) {
      int iHit = itTrackHits.first;
      STHit& hit = vHits[iHit];
      int iLayer = hit.GetLayer();
      hitIndex[iLayer] = iHit;
    }

    int nLayers = 0;
    for (int iLayer : hitIndex) {
      if (iLayer < 0) {
        continue;
      }
      nLayers++;
    }
    // reject tracks with less than 3 layers
    if (nLayers < 3) {
      continue;
    }

    // create tracks and store them
    STTrack track;
    for (int iHit : hitIndex) {
      if (iHit < 0) {
        continue;
      }
      track.addHitIndex(iHit);
    }
    track.setMCTrackID(itMap.first);
    vTracks.push_back(track);
  }
}

void STEventAction::TrackFitter()
{
  // Using simple helix model
  // Ignoring multiple scattering effects
  // zx plane:
  // (x-x0)^2 + (z-z0)^2 = R^2
  TF1 funcCircle(
    "funcCircle",
    "[3]*sqrt([0]^2-(x-[1])^2)+[2]", 0, 1200); // x(z)
  TF1 funcLinear(
    "funcLinear",
    "[0]+[1]*x", 0, 1200); // y(z)
  for (uint iTrack = 0; iTrack < vTracks.size(); iTrack++) {
    STTrack& track = vTracks[iTrack];
    int nHits = track.getNHitIndices();
    TGraphErrors graphZX(nHits); // circle
    TGraphErrors graphZY(nHits); // straight line

    for (int iHit = 0; iHit < nHits; iHit++) {
      int hitIndex = track.getHitIndex(iHit);
      STHit& hit = vHits[hitIndex];
      graphZX.SetPoint(iHit, hit.Z(), hit.X());
      graphZY.SetPoint(iHit, hit.Z(), hit.Y());
      graphZX.SetPointError(iHit, hit.Dz(), hit.Dx());
      graphZY.SetPointError(iHit, hit.Dz(), hit.Dy());
    }
    // using circle with first 3 hits to deduce initial x0,y0 and R values
    // http://algolist.ru/maths/geom/equation/circle.php
    double z1 = graphZX.GetX()[0];
    double x1 = graphZX.GetY()[0];
    double z2 = graphZX.GetX()[1];
    double x2 = graphZX.GetY()[1];
    double z3 = graphZX.GetX()[2];
    double x3 = graphZX.GetY()[2];
    double ma = (x2 - x1) / (z2 - z1);
    double mb = (x3 - x2) / (z3 - z2);
    double z0 = (ma * mb * (x1 - x3) + mb * (z1 + z2) - ma * (z2 + z3)) / 2. / (mb - ma);
    double x0 = (x1 + x2) / 2. - (z0 - (z1 + z2) / 2.) / ma;
    double r0 = sqrt((z1 - z0) * (z1 - z0) + (x1 - x0) * (x1 - x0));
    funcCircle.SetParameters(r0, z0, x0, 1.);
    funcCircle.FixParameter(3, 1.);
    if (fabs(funcCircle.Eval(z1) - x1) > 1e-3)
      funcCircle.FixParameter(3, -1);
    if (fabs(funcCircle.Eval(z1) - x1) > 1e-3)
      printf(
        "ERROR\n");
    graphZX.Fit(&funcCircle, "Q");
    double r = funcCircle.GetParameter(0);
    double zc = funcCircle.GetParameter(1);
    double xc = funcCircle.GetParameter(2);
    int q = std::round(funcCircle.GetParameter(3));
    // track parameters at z=0
    double z = 0;
    double x = funcCircle.Eval(z);    // x coordinate at 0.
    double kx = -(z - zc) / (x - xc); // kx = px/pz
    double B = 1.;                    // 1 Tesla magnetic field
    // using pt^2 = px^2 + pz^2
    // R(m)  = pt (GeV)/0.3/B(T)
    // R(mm) = pt (MeV)/0.3/B(T)
    // pt(MeV) = R(mm)*0.3*B(T)
    double pt = 0.3 * B * r;
    double pz = pt / sqrt(1 + kx * kx);
    double px = kx * pz;

    // determine y and ky=py/pz from ZY projection
    graphZY.Fit(&funcLinear, "Q");
    double ky = funcLinear.GetParameter(1);
    double y = funcLinear.Eval(z);
    double py = ky * pz;

    printf("Track %i:", iTrack);
    printf(" at z=%.2f", z);
    printf(" (x,y,px,py,pz)=(%.2f,%.2f,%7.2f,%7.2f,%7.2f)", x, y, px, py, pz);
    printf(" q=%i\n", q);
    printf("\n");
    track.setState(z, x, y, px, py, pz);

    fHistoP->Fill(track.getP());

    if (fDebug) {
      TCanvas c1;
      graphZX.Draw("ap");
      c1.Print("zx.png");

      TCanvas c2;
      graphZY.Draw("ap");
      c2.Print("zy.png");
    }
  }
}

int iCanvasFit = 0;
int iCanvasRefit = 0;
int iEvent = 1;
void STEventAction::fitTracksKF()
{
  gRandom = new TRandomMT64();
  gRandom->SetSeed(std::time(nullptr));

  // todo: set appropriate covariance matrix
  TMatrixT<double> measCovMatrix(2, 2);
  // sigma^2 x
  measCovMatrix(0, 0) = 0.04;
  // sigma^2 y
  measCovMatrix(1, 1) = 0.04;

  // initial covariance matrix
  TMatrixT<double> initCovMatrix(5, 5);
  for (uint i = 0; i < 5; i++) {
    initCovMatrix(i, i) = 100000.;
  }

  auto* trackFitter = new STTrackFitter();
  trackFitter->setMagFieldHomogeneous(eHomoMagField[0], eHomoMagField[1], eHomoMagField[2]);
  trackFitter->setCalcLoss(true);
  trackFitter->setDebugLevel(0);

  // double mass = 0.5109989461; // electron mass, [mev]
  double mass = 105.6583745;  // muon mass, [mev]
  trackFitter->setMass(mass); // [mev]; // todo : set proper hypothesis
  trackFitter->setCharge(1.);

  // todo : get material parameters from geant4 geometry
  // parameters for Si from pdg
  // https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/silicon_Si.html
  trackFitter->setMaterialPars(14.,     // [charge]
                               28.0855, // [at. units]
                               2.3290,  // [g/cm^3]
                               173.,    // [ev]
                               9.370);  // [cm]

  std::vector<double> layers(eNSiLayers);
  double layerZ = eLayersDist * 10.;
  // central coordinates
  for (auto& layer : layers) {
    layer = layerZ;
    layerZ += eLayersDist * 10.;
  }
  trackFitter->setDetectorGeom(layers);
  trackFitter->setLayerTH(eLayersThic);

  std::vector<STMcPoint> trackMcPoints;
  int nTracks = vTracks.size();
  for (int iTrack = 0; iTrack < nTracks; iTrack++) {
    auto& track = vTracks[iTrack];

    // collect relevant mc points
    for (auto mcPoint : vMcPoints) {
      if (mcPoint.GetTrackID() == track.getMCTrackID()) {
        trackMcPoints.push_back(mcPoint);
      }
    }

    // construct measurements out of mc points
    std::vector<std::vector<double>> measVec(trackMcPoints.size(), std::vector<double>(2, 0));
    std::vector<double> coordsZ(trackMcPoints.size(), 0);

    // construct measurements out of mc points
    for (uint iPoint = 0; iPoint < trackMcPoints.size(); iPoint++) {
      auto& mcPoint = trackMcPoints[iPoint];
      double measX = mcPoint.GetPosOut().getX() + gRandom->Gaus(0., 0.2);
      double measY = mcPoint.GetPosOut().getY() + gRandom->Gaus(0., 0.2);
      measVec[iPoint][0] = measX;
      measVec[iPoint][1] = measY;
      double measZ = mcPoint.GetPosOut().getZ();
      coordsZ[iPoint] = measZ;
    }

    // debug
    for (uint i = 0; i < measVec.size(); i++) {
      printf("LOG(INFO): STEventAction::fitTracksKF(): z = %.5f, m[0] = %.5f, m[1] = %.5f\n",
             coordsZ[i], measVec[i][0], measVec[i][1]);
    }

    track.setStateFromKF(0, std::vector<double>{0, 0, 0, 0, 1e-6});
    track.setCovMatrix(initCovMatrix);
    STTrack fittedTrack; // todo: propagate fitted tracks to vTrack vector
    std::vector<std::vector<double>> fitStates;
    trackFitter->setRefit(false);
    trackFitter->setZStepLength(0.01);
    trackFitter->fitTrack(track, coordsZ, measVec, measCovMatrix, fitStates, fittedTrack);

    // draw trajectory if needed
    fDrawTraj = false;
    if (fDrawTraj) {
      drawQA(track, coordsZ, fitStates, trackFitter->fTrajectoryZXrk4, Form("kf_qa_fit_%d", iCanvasFit));
      iCanvasFit++;
    }
    trackFitter->clearTrajectory();

    for (uint iState = 0; iState < fitStates.size(); iState++) {
      printf("LOG(INFO): STEventAction::fitTracksKF(): fitted state: z = %.5f, [x] = %.5f, [y] = %.5f, [tx] = %.5f, [ty] = %.5f, [q/p] = %.5f\n",
             coordsZ[iState], fitStates[iState][0], fitStates[iState][1], fitStates[iState][2], fitStates[iState][3], fitStates[iState][4]);
    }
    printf("\nLOG(INFO): STEventAction::fitTracksKF(): track chi2: %.3f\n", fittedTrack.getChi2());

    // pull histograms at the last position
    fPulls = false;
    if (fPulls) {
      TMatrixT<double> cov(5, 5);
      fittedTrack.getCovMatrix(cov);
      auto& mcPoint = trackMcPoints.back();
      double mcX = mcPoint.GetPosOut().getX();
      double mcY = mcPoint.GetPosOut().getY();
      double mcTx = (mcPoint.GetPosOut().getX() - mcPoint.GetPosIn().getX()) / (mcPoint.GetPosOut().getZ() - mcPoint.GetPosIn().getZ());
      double mcTy = (mcPoint.GetPosOut().getY() - mcPoint.GetPosIn().getY()) / (mcPoint.GetPosOut().getZ() - mcPoint.GetPosIn().getZ());
      double recoX = fitStates.back()[0];
      double recoY = fitStates.back()[1];
      double recoTx = fitStates.back()[2];
      double recoTy = fitStates.back()[3];
      fRunAction->hXResiduals->Fill(recoX - mcX);
      fRunAction->hYResiduals->Fill(recoY - mcY);
      fRunAction->hXPulls->Fill((recoX - mcX) / std::sqrt(cov(0, 0)));
      fRunAction->hYPulls->Fill((recoY - mcY) / std::sqrt(cov(1, 1)));
      double mcQP = 1. / eParticleMom; // [1/mev]
      double recoQP = fitStates.back()[4];
      fRunAction->hPResiduals->Fill(recoQP - mcQP);
      fRunAction->hPPulls->Fill((recoQP - mcQP) / std::sqrt(cov(4, 4)));
      fRunAction->hPSigmas->Fill(std::sqrt(cov(4, 4)));
      fRunAction->hChi2s->Fill(fittedTrack.getChi2());
      fRunAction->hTxResiduals->Fill(recoTx - mcTx);
      fRunAction->hTyResiduals->Fill(recoTy - mcTy);
      fRunAction->hTxPulls->Fill((recoTx - mcTx) / std::sqrt(cov(2, 2)));
      fRunAction->hTyPulls->Fill((recoTy - mcTy) / std::sqrt(cov(3, 3)));
    }

    double sumMomLoss = trackFitter->getSumMomLoss();
    printf("\nLOG(INFO): sum momentum loss = %.3f MeV\n\n", sumMomLoss);

    fitStates.clear();

    // ---------------------------------------------------------------
    // refitting track
    trackFitter->extrapolateTrackToZ(fittedTrack, coordsZ.back() + 50.);
    for (int iPoint = trackMcPoints.size() - 1; iPoint >= 0; iPoint--) {
      auto& mcPoint = trackMcPoints[iPoint];
      int iMeas = trackMcPoints.size() - 1 - iPoint;
      double measX = mcPoint.GetPosIn().getX() + gRandom->Gaus(0., 0.2);
      double measY = mcPoint.GetPosIn().getY() + gRandom->Gaus(0., 0.2);
      measVec[iMeas][0] = measX;
      measVec[iMeas][1] = measY;
      double measZ = mcPoint.GetPosIn().getZ();
      coordsZ[iMeas] = measZ;
    }

    // extrapolating to (0, 0, 0)
    measVec.emplace_back(2, 0);
    coordsZ.emplace_back(0);

    std::vector<std::vector<double>> refitStates;
    STTrack refittedTrack;

    fittedTrack.setCovMatrix(initCovMatrix);
    trackFitter->setZStepLength(-0.01);
    trackFitter->setRefit(true);
    trackFitter->fitTrack(fittedTrack, coordsZ, measVec, measCovMatrix, refitStates, refittedTrack);

    fDrawTraj = false;
    if (fDrawTraj) {
      drawQA(track, coordsZ, refitStates, trackFitter->fTrajectoryZXrk4, Form("kf_qa_fit_%d", iCanvasFit));
      iCanvasFit++;
    }
    trackFitter->clearTrajectory();

    for (uint iState = 0; iState < refitStates.size(); iState++) {
      printf("LOG(INFO): STEventAction::fitTracksKF(): refit state: z = %.5f, [x] = %.5f, [y] = %.5f, [tx] = %.5f, [ty] = %.5f, [q/p] = %.5f\n",
             coordsZ[iState], refitStates[iState][0], refitStates[iState][1], refitStates[iState][2], refitStates[iState][3], refitStates[iState][4]);
    }
    printf("\nLOG(INFO): STEventAction::fitTracksKF(): track chi2: %.3f\n", std::abs(refittedTrack.getChi2()));

    // pull histograms at the first position
    fPulls = false;
    if (fPulls) {
      TMatrixT<double> cov(5, 5);
      refittedTrack.getCovMatrix(cov);
      double mcX = 0;
      double mcY = 0;
      double mcTx = 0;
      double mcTy = 0;
      double recoX = refitStates.back()[0];
      double recoY = refitStates.back()[1];
      double recoTx = refitStates.back()[2];
      double recoTy = refitStates.back()[3];
      fRunAction->hXResiduals->Fill(recoX - mcX);
      fRunAction->hYResiduals->Fill(recoY - mcY);
      fRunAction->hXPulls->Fill((recoX - mcX) / std::sqrt(cov(0, 0)));
      fRunAction->hYPulls->Fill((recoY - mcY) / std::sqrt(cov(1, 1)));
      double mcQP = 1. / eParticleMom; // [1/mev]
      double recoQP = refitStates.back()[4];
      fRunAction->hPResiduals->Fill(recoQP - mcQP);
      fRunAction->hPPulls->Fill((recoQP - mcQP) / std::sqrt(cov(4, 4)));
      fRunAction->hPSigmas->Fill(std::sqrt(cov(4, 4)));
      fRunAction->hChi2s->Fill(refittedTrack.getChi2());
      fRunAction->hTxResiduals->Fill(recoTx - mcTx);
      fRunAction->hTyResiduals->Fill(recoTy - mcTy);
      fRunAction->hTxPulls->Fill((recoTx - mcTx) / std::sqrt(cov(2, 2)));
      fRunAction->hTyPulls->Fill((recoTy - mcTy) / std::sqrt(cov(3, 3)));
    }

    // ---------------------------------------------------------------
    // write tracks into tree
    fRunAction->recoTrack.eventID = iEvent;
    fRunAction->recoTrack.trackID = iTrack + 1;
    fRunAction->recoTrack.mcLabel = track.getMCTrackID();
    fRunAction->recoTrack.z = coordsZ.back();
    fRunAction->recoTrack.x = refitStates.back()[0];
    fRunAction->recoTrack.y = refitStates.back()[1];
    fRunAction->recoTrack.tx = refitStates.back()[2];
    fRunAction->recoTrack.ty = refitStates.back()[3];
    fRunAction->recoTrack.qp = refitStates.back()[4];
    refittedTrack.getCovMatrix(fRunAction->recoTrack.covM);
    fRunAction->tTracks->Fill();

    trackMcPoints.clear();
    refitStates.clear();
    measVec.clear();
    coordsZ.clear();

    sumMomLoss = trackFitter->getSumMomLoss();
    printf("\nLOG(INFO): sum momentum loss = %.3f MeV\n\n", sumMomLoss);
  }

  iEvent++;

  delete gRandom;
  delete trackFitter;
}

/**
 * draw QA figure with hits and reconstructed states
 */
void STEventAction::drawQA(STTrack& preFitTrack,
                           std::vector<double>& coordsZ,
                           std::vector<std::vector<double>>& recoStates,
                           TGraph* trajectory,
                           const std::string& filename)
{
  auto* graphZX = new TGraph();
  graphZX->SetMarkerStyle(21);
  graphZX->SetMarkerSize(1);
  graphZX->SetLineWidth(0);
  graphZX->SetMarkerColor(kRed);

  int nHitsTest = preFitTrack.getNHitIndices();
  for (uint iHit = 0; iHit < nHitsTest; iHit++) {
    int hitIndex = preFitTrack.getHitIndex(iHit);
    STHit& hit = vHits[hitIndex];
    graphZX->SetPoint(iHit, hit.Z(), hit.X());
  }

  auto* graphZX_KF = new TGraph();
  graphZX_KF->SetMarkerStyle(21);
  graphZX_KF->SetMarkerSize(2);
  graphZX_KF->SetLineWidth(0);
  graphZX_KF->SetMarkerColor(kGreen);

  for (uint iState = 0; iState < recoStates.size(); iState++) {
    graphZX_KF->SetPoint(iState, coordsZ[iState], recoStates[iState][0]);
  }

  auto* magmg = new TMultiGraph();
  graphZX_KF->GetYaxis()->SetLimits(-1000, 500);
  graphZX_KF->GetXaxis()->SetLimits(-100, 1100);
  magmg->Add(trajectory);
  magmg->Add(graphZX_KF);
  magmg->Add(graphZX);
  magmg->SetTitle((filename + ";z [mm]; x [mm];").c_str());

  auto* trajectories = new TFile("trajectories.root", "update");
  trajectories->cd();
  magmg->Write(filename.c_str());
  trajectories->Close();
  delete trajectories;
  delete graphZX_KF;
  delete graphZX;
}
