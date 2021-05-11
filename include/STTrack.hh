#ifndef S1Track_h
#define S1Track_h 1

#include <vector>
#include "TVector3.h"
#include "TMatrixT.h"

class STTrack
{
 public:
  STTrack();
  virtual ~STTrack();

  void addHitIndex(int index) { fHitIndices.push_back(index); }

  int getNHitIndices() { return fHitIndices.size(); }

  int getHitIndex(int i) { return fHitIndices[i]; }

  // set common parameters
  void setState(double z, double x, double y, double px, double py, double pz);

  // common parameters: z, x, y, px, py, pz
  void getState(std::vector<double>& state);

  void setCovMatrix(const TMatrixT<double>& covMatrix);

  void getCovMatrix(TMatrixT<double>& covMatrix);

  // array size must be consistent
  void getCovMatrix(double* covMatrixLT) { std::copy(fCov, fCov + 15, covMatrixLT); }

  double getZ() { return fZ; };

  double getCharge() { return fCharge; };

  // parameters for KF: x, y, tx, ty, q/p
  void getStateKF(std::vector<double>& stateVector);

  // calculate common parameters from KF state
  void setStateFromKF(double z, std::vector<double> stateVectorKF);

  double getP() { return std::sqrt(fV[2] * fV[2] + fV[3] * fV[3] + fV[4] * fV[4]); };

  double getChi2() { return fChi2; };

  void setChi2(double chi2) { fChi2 = chi2; };

  void setMCTrackID(int mcTrackID) { fMCTrackID = mcTrackID; };

  int getMCTrackID() { return fMCTrackID; };

 private:
  double fZ;
  double fV[5]; // x,y,px,py,pz at given z
  double fCov[15]{};
  double fChi2{};
  double fCharge{1}; // todo: set real charge
  double fQP{};
  int fMCTrackID;
  std::vector<int> fHitIndices;
};

#endif
