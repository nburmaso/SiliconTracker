//
// Created by nburmaso on 3/30/21.
//

#ifndef SILICONTRACKER_STTRACKFITTER_HH
#define SILICONTRACKER_STTRACKFITTER_HH

#include "STTrack.hh"

#include <vector>
#include "TGraph.h"
#include "TMatrixT.h"
#include "TRandomGen.h"

class STTrackFitter
{
 public:
  STTrackFitter();

  virtual ~STTrackFitter();

  // right side of the differential equation
  void extrapEquation(std::vector<double>& f, double tx, double ty, double qp);

  // RK4 implementation for a particle in a magnetic field
  void integrateMagRK4(std::vector<double>& stateVectorKF,
                       std::vector<double>& nextStateVectorKF);

  // Bethe-Bloch
  double eLossBetheBloch(double beta2, double gamma, double gamma2) const;

  // Bethe-Bloch correction for noise matrix
  void noiseBetheBloch(TMatrixT<double>& noiseMat, double p, double beta2, double gamma, double gamma2) const;

  // energy loss implementation
  double energyLoss(double energy, double p) const;

  // momentum loss implementation
  double momentumLoss(double p, double step);

  // check if inside a material
  void checkMaterial(double z);

  // propagate particle with RK4
  void extrapolateTrackToZ(STTrack& s1Track, double nextZ);

  // propagate covariance matrix using jacobian of system of ODEs
  void propagateCovMatrix(const std::vector<double>& stateKF,
                          TMatrixT<double>& covMatrix);

  void fitTrack(STTrack& inS1Track,
                const std::vector<double>& vZcoords,
                const std::vector<std::vector<double>>& vMeas,
                const TMatrixT<double>& covMatrixMeas,
                std::vector<std::vector<double>>& filteredStates,
                STTrack& outS1Track,
                bool isRefit = false);

  // z step length in mm
  void setZStepLength(double hz) { fDz = hz; }

  void setMagFieldHomogeneous(double bx, double by, double bz)
  {
    fMagFieldHomogeneous[0] = bx;
    fMagFieldHomogeneous[1] = by;
    fMagFieldHomogeneous[2] = bz;
  }

  /*
   * input: vector of z coordinates in [mm]
   */
  void setDetectorGeom(std::vector<double>& detZCoords)
  {
    if (!fDetectorLayers.empty()) {
      fDetectorLayers.clear();
    }
    for (auto& item : detZCoords) {
      fDetectorLayers.push_back(item);
    }
  }

  // half-thickness in [mm]
  void setLayerTH(double thick) { fLayerTH = thick; }

  // set mass in [mev]
  void setMass(double mass) { fMass = mass; }

  void setCharge(double charge) { fCharge = charge; }

  // set material parameters
  void setMaterialPars(double matZ, double matA, double matDensity, double matExEnrg, double radLength)
  {
    fMatZ = matZ;
    fMatA = matA;
    fMatDensity = matDensity;
    fMatExEnrg = matExEnrg;
    fMatRadL = radLength;
  }

  void setCalcLoss(bool calcLoss)
  {
    fCalcLoss = calcLoss;
    fELossBetheBloch = calcLoss;
  }

  // print debug output: intermediate matrices, states etc. ...
  void setDebugLevel(uint level) { fDebugLevel = level; }

  double getSumMomLoss() { return fSumMomLoss; }

  void clearTrajectory();

  // trajectory of a particle
  // only for debug purposes
  // todo: remove or move to private
  TGraph* fTrajectoryZXrk4;

 private:
  // hypothetical particle parameters
  double fMass;
  double fCharge;

  // material parameters
  double fMatZ;       // material z
  double fMatA;       // material a
  double fMatDensity; // material density
  double fMatExEnrg;  // mean excitation energy
  double fMatRadL;    // rad. length

  // detector geometry
  std::vector<double> fDetectorLayers; // z coordinates of detector layers
  double fLayerTH; // layer thickness

  // sum momentum loss : for debug purposes
  double fSumMomLoss{0};

  // internal service variable:
  // energy loss per length unit
  double fDEdx;

  // extrapolation parameters
  double fMagFieldHomogeneous[3];
  double fDz{1};

  // service flags
  short fDebugLevel{0};
  bool fCalcLoss{false};
  bool fInMaterial{false};
  bool fELossBetheBloch{false};
  bool fELossBrems{false};

  TRandom* fRandomGen;
};

#endif //SILICONTRACKER_STTRACKFITTER_HH
