#include "STTrack.hh"

STTrack::STTrack()
  : fZ(0),
    fV{0, 0, 0, 0, 0},
    fChi2(0),
    fHitIndices()
{
}

STTrack::~STTrack()
{
  fHitIndices.clear();
}

void STTrack::setState(double z, double x, double y, double px, double py, double pz)
{
  fZ = z;
  fV[0] = x;
  fV[1] = y;
  fV[2] = px;
  fV[3] = py;
  fV[4] = pz;
}

void STTrack::getStateKF(std::vector<double>& stateVector)
{
  stateVector[0] = fV[0];         // x
  stateVector[1] = fV[1];         // y
  stateVector[2] = fV[2] / fV[4]; // tx
  stateVector[3] = fV[3] / fV[4]; // ty
  stateVector[4] = fQP;           //fCharge / getP(); // q/p
}

void STTrack::getState(std::vector<double>& state)
{
  state[0] = fZ;    // z
  state[1] = fV[0]; // x
  state[2] = fV[1]; // y
  state[3] = fV[2]; // px
  state[4] = fV[3]; // py
  state[5] = fV[4]; // pz
}

void STTrack::setCovMatrix(const TMatrixT<double> covMatrix)
{
  // sigma^2 x
  fCov[0] = covMatrix(0, 0);

  fCov[1] = 0.5f * covMatrix(1, 0);

  // sigma^2 y
  fCov[2] = covMatrix(1, 1);

  fCov[3] = 0.5f * covMatrix(2, 0);

  fCov[4] = 0.5f * covMatrix(2, 1);

  // sigma^2 tx
  fCov[5] = covMatrix(2, 2);

  fCov[6] = 0.5f * covMatrix(3, 0);

  fCov[7] = 0.5f * covMatrix(3, 1);

  fCov[8] = 0.5f * covMatrix(3, 2);

  // sigma^2 ty
  fCov[9] = covMatrix(3, 3);

  fCov[10] = 0.5f * covMatrix(4, 0);

  fCov[11] = 0.5f * covMatrix(4, 1);

  fCov[12] = 0.5f * covMatrix(4, 2);

  fCov[13] = 0.5f * covMatrix(4, 3);

  // sigma^2 q/p
  fCov[14] = covMatrix(4, 4);
}

void STTrack::getCovMatrix(TMatrixT<double>& covMatrix)
{
  // sigma^2 x
  covMatrix(0, 0) = fCov[0];

  covMatrix(1, 0) = fCov[1] * 2;
  covMatrix(0, 1) = fCov[1] * 2;

  // sigma^2 y
  covMatrix(1, 1) = fCov[2];

  covMatrix(2, 0) = fCov[3] * 2;
  covMatrix(2, 1) = fCov[4] * 2;
  covMatrix(0, 2) = fCov[3] * 2;
  covMatrix(1, 2) = fCov[4] * 2;

  // sigma^2 tx
  covMatrix(2, 2) = fCov[5];

  covMatrix(3, 0) = fCov[6] * 2;
  covMatrix(3, 1) = fCov[7] * 2;
  covMatrix(3, 2) = fCov[8] * 2;
  covMatrix(0, 3) = fCov[6] * 2;
  covMatrix(1, 3) = fCov[7] * 2;
  covMatrix(2, 3) = fCov[8] * 2;

  // sigma^2 ty
  covMatrix(3, 3) = fCov[9];

  covMatrix(4, 0) = fCov[10] * 2;
  covMatrix(4, 1) = fCov[11] * 2;
  covMatrix(4, 2) = fCov[12] * 2;
  covMatrix(4, 3) = fCov[13] * 2;
  covMatrix(0, 4) = fCov[10] * 2;
  covMatrix(1, 4) = fCov[11] * 2;
  covMatrix(2, 4) = fCov[12] * 2;
  covMatrix(3, 4) = fCov[13] * 2;

  // sigma^2 q/p
  covMatrix(4, 4) = fCov[14];
}

void STTrack::setStateFromKF(double z, std::vector<double> stateVectorKF)
{
  double p = 1 / stateVectorKF[4] / getCharge();
  double pz = std::sqrt(p * p / (1 + stateVectorKF[2] * stateVectorKF[2] + stateVectorKF[3] * stateVectorKF[3]));
  double px = pz * stateVectorKF[2];
  double py = pz * stateVectorKF[3];
  double qp = stateVectorKF[4];
  setState(z, stateVectorKF[0], stateVectorKF[1], px, py, pz);
  fQP = qp;
}