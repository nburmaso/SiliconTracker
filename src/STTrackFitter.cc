//
// Created by nburmaso on 3/30/21.
//
#include "STTrackFitter.hh"

#include "TMatrixT.h"
#include "TRandom.h"
#include "TRandomGen.h"

#include "TGraph.h"
#include "TCanvas.h"

const double me = 0.5109989461; // electron mass [mev]

STTrackFitter::STTrackFitter() : fMagFieldHomogeneous{0, 0, 0}
{
  fRandomGen = new TRandomMT64();
  fRandomGen->SetSeed(std::time(nullptr));

  // debug
  fTrajectoryZXrk4 = new TGraph();
  fTrajectoryZXrk4->SetMarkerStyle(21);
  fTrajectoryZXrk4->SetMarkerSize(0.5);
  fTrajectoryZXrk4->SetLineWidth(0);
  fTrajectoryZXrk4->SetMarkerColor(kBlue);
}

STTrackFitter::~STTrackFitter()
{
  delete fRandomGen;
  delete fTrajectoryZXrk4;
  fDetectorLayers.clear();
}

void STTrackFitter::extrapEquation(std::vector<double>& f, double tx, double ty, double qp)
{
  double a = 2.99792458 * 1e-1; // [MeV/c] * [Tesla^-1] * [mm^-1]
  f[0] = tx;
  f[1] = ty;
  f[2] = a * qp * std::sqrt(1 + tx * tx + ty * ty) *
         (tx * ty * fMagFieldHomogeneous[0] - (1 + tx * tx) * fMagFieldHomogeneous[1] + ty * fMagFieldHomogeneous[2]);
  f[3] = a * qp * std::sqrt(1 + tx * tx + ty * ty) *
         ((1 + ty * ty) * fMagFieldHomogeneous[0] - tx * ty * fMagFieldHomogeneous[1] - tx * fMagFieldHomogeneous[2]);
  f[4] = 0;
}

/*
 * implementation of RK4 for a particle in a magnetic field
 */
void STTrackFitter::integrateMagRK4(std::vector<double>& stateVectorKF,
                                    std::vector<double>& nextStateVectorKF)
{
  std::vector<double> k1(5, 0);
  std::vector<double> k2(5, 0);
  std::vector<double> k3(5, 0);
  std::vector<double> k4(5, 0);

  // right side of the equation
  std::vector<double> vecFunc(5, 0);

  // calculating k1
  extrapEquation(vecFunc,
                 stateVectorKF[2],
                 stateVectorKF[3],
                 stateVectorKF[4]);
  for (uint i = 0; i < 5; i++) {
    k1[i] = fDz * vecFunc[i];
  }

  // calculating k2
  extrapEquation(vecFunc,
                 stateVectorKF[2] + 0.5 * k1[2],
                 stateVectorKF[3] + 0.5 * k1[3],
                 stateVectorKF[4] + 0.5 * k1[4]);
  for (uint i = 0; i < 5; i++) {
    k2[i] = fDz * vecFunc[i];
  }

  // calculating k3
  extrapEquation(vecFunc,
                 stateVectorKF[2] + 0.5 * k2[2],
                 stateVectorKF[3] + 0.5 * k2[3],
                 stateVectorKF[4] + 0.5 * k2[4]);
  for (uint i = 0; i < 5; i++) {
    k3[i] = fDz * vecFunc[i];
  }

  // calculating k4
  extrapEquation(vecFunc,
                 stateVectorKF[2] + k3[2],
                 stateVectorKF[3] + k3[3],
                 stateVectorKF[4] + k3[4]);
  for (uint i = 0; i < 5; i++) {
    k4[i] = fDz * vecFunc[i];
  }

  // calculating next state
  for (uint i = 0; i < 5; i++) {
    nextStateVectorKF[i] = stateVectorKF[i] + 1. / 6. * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
  }
}

double STTrackFitter::eLossBetheBloch(double beta2, double gamma, double gamma2) const
{
  double eloss = 0.307075 * fMatZ / fMatA * fMatDensity / beta2 * fCharge * fCharge;
  double massRatio = me / fMass;
  double argument = gamma2 * beta2 * me * 2. / ((1e-6 * fMatExEnrg) * std::sqrt(1. + 2. * gamma * massRatio + massRatio * massRatio));
  eloss *= std::log(argument) - beta2; // Bethe-Bloch [MeV/cm]

  if (eloss < 0.) {
    eloss = 0.;
  }

  return eloss;
}

// from genfit
void STTrackFitter::noiseBetheBloch(TMatrixT<double>& noiseMat, double p, double beta2, double gamma, double gamma2) const
{
  // Code ported from GEANT 3 (erland)

  // ENERGY LOSS FLUCTUATIONS; calculate sigma^2(E);
  double sigma2E = 0.;
  double zeta = 153.4e3 * fCharge * fCharge / beta2 * fMatZ / fMatA * fMatDensity * std::abs(fDz / 10.);   // eV
  double Emax = 2.e9 * me * beta2 * gamma2 / (1. + 2. * gamma * me / fMass + (me / fMass) * (me / fMass)); // eV
  double kappa = zeta / Emax;

  if (kappa > 0.01) {                           // Vavilov-Gaussian regime
    sigma2E += zeta * Emax * (1. - beta2 / 2.); // eV^2
  } else {                                      // Urban/Landau approximation
    // calculate number of collisions Nc
    double I = 16. * std::pow(fMatZ, 0.9); // eV
    double f2 = 0.;
    if (fMatZ > 2.) {
      f2 = 2. / fMatZ;
    }
    double f1 = 1. - f2;
    double e2 = 10. * fMatZ * fMatZ;            // eV
    double e1 = pow((I / pow(e2, f2)), 1 / f1); // eV

    double mbbgg2 = 2.e9 * fMass * beta2 * gamma2;                                                          // eV
    double Sigma1 = fDEdx * 1.0e9 * f1 / e1 * (log(mbbgg2 / e1) - beta2) / (log(mbbgg2 / I) - beta2) * 0.6; // 1/cm
    double Sigma2 = fDEdx * 1.0e9 * f2 / e2 * (log(mbbgg2 / e2) - beta2) / (log(mbbgg2 / I) - beta2) * 0.6; // 1/cm
    double Sigma3 = fDEdx * 1.0e9 * Emax / (I * (Emax + I) * log((Emax + I) / I)) * 0.4;                    // 1/cm

    double Nc = (Sigma1 + Sigma2 + Sigma3) * std::abs(fDz / 10);

    if (Nc > 50.) { // truncated Landau distribution
      double sigmaalpha = 15.76;
      // calculate sigmaalpha  (see GEANT3 manual W5013)
      double RLAMED = -0.422784 - beta2 - log(zeta / Emax);
      double RLAMAX = 0.60715 + 1.1934 * RLAMED + (0.67794 + 0.052382 * RLAMED) * std::exp(0.94753 + 0.74442 * RLAMED);
      // from lambda max to sigmaalpha=sigma (empirical polynomial)
      if (RLAMAX <= 1010) {
        sigmaalpha = 1.975560 + 9.898841e-02 * RLAMAX - 2.828670e-04 * RLAMAX * RLAMAX + 5.345406e-07 * std::pow(RLAMAX, 3.) - 4.942035e-10 * std::pow(RLAMAX, 4) + 1.729807e-13 * std::pow(RLAMAX, 5);
      } else {
        sigmaalpha = 1.871887e+01 + 1.296254e-02 * RLAMAX;
      }
      // alpha=54.6  corresponds to a 0.9996 maximum cut
      if (sigmaalpha > 54.6) {
        sigmaalpha = 54.6;
      }
      sigma2E += sigmaalpha * sigmaalpha * zeta * zeta; // eV^2
    } else {                                            // Urban model
      static const double alpha = 0.996;
      double Ealpha = I / (1. - (alpha * Emax / (Emax + I)));                                    // eV
      double meanE32 = I * (Emax + I) / Emax * (Ealpha - I);                                     // eV^2
      sigma2E += std::abs(fDz / 10.) * (Sigma1 * e1 * e1 + Sigma2 * e2 * e2 + Sigma3 * meanE32); // eV^2
    }
  }

  sigma2E *= 1.e-15; // [eV] -> [MeV]

  // update noise matrix, using linear error propagation from E to q/p
  noiseMat(4, 4) += fCharge * fCharge / beta2 / std::pow(p, 4.) * sigma2E;
}

/*
 * energy loss per length unit
 */
double STTrackFitter::energyLoss(double energy, double p) const
{
  double eloss = 0;
  double gamma = energy / fMass;
  double gamma2 = gamma * gamma;
  double beta2 = 1. - 1. / gamma2;

  if (fELossBetheBloch) {
    eloss += eLossBetheBloch(beta2, gamma, gamma2);
  }

  // todo : implement Bremsstrahlung
  /*
  if (fELossBrems) {
    eloss += eLossBrems(p);
  }
  */

  return eloss;
}

// from Genfit
double STTrackFitter::momentumLoss(double p, double step)
{
  double E0 = std::sqrt(p * p + fMass * fMass);
  double dEdx1 = energyLoss(E0, p);

  // RK4
  double E1 = E0 - dEdx1 * step / 10. / 2.;
  double dEdx2 = energyLoss(E1, p); // dEdx(x0 + h/2, E0 + h/2 * dEdx1)

  double E2 = E0 - dEdx2 * step / 10. / 2.;
  double dEdx3 = energyLoss(E2, p); // dEdx(x0 + h/2, E0 + h/2 * dEdx2)

  double E3 = E0 - dEdx3 * step / 10.;
  double dEdx4 = energyLoss(E3, p); // dEdx(x0 + h, E0 + h * dEdx3)

  fDEdx = (dEdx1 + 2. * dEdx2 + 2. * dEdx3 + dEdx4) / 6.;

  double dE = step / 10. * fDEdx; // positive for positive stepSign

  double momLoss = 0.;

  if (E0 - dE <= fMass) {
    // Step would stop particle (E_kin <= 0).
    return momLoss = p;
  } else {
    momLoss = p - std::sqrt((E0 - dE) * (E0 - dE) - fMass * fMass); // momLoss; positive for positive step
  }

  if (fDebugLevel > 2) {
    printf("LOG(INFO): STTrackFitter::momentumLoss: p = %.3f, E0 = %.3f, dEdx = %.6f, dE = %.6f, mass = %.3f\n",
           p, E0, fDEdx, dE, fMass);
  }

  return momLoss; // loss in [mev]
}

void STTrackFitter::checkMaterial(double z)
{
  fInMaterial = false;
  for (auto& coord : fDetectorLayers) {
    if (std::abs(coord - z) < fLayerTH) {
      fInMaterial = true;
      break;
    }
  }
}

/**
 * Extrapolating track parameters using RK4 integrator
 * and propagating covariance matrix
 * @param s1Track : track to update
 * @param nextZ : z coordinate to extrapolate to
 */
void STTrackFitter::extrapolateTrackToZ(STTrack& s1Track, double nextZ)
{
  if (fDebugLevel > 0) {
    printf("LOG(INFO): extrapolate: in\n");
  }

  // state vector: x, y, tx, ty, q/p
  std::vector<double> trackStateKF(5, 0);
  std::vector<double> nextStateKF(5, 0);
  s1Track.getStateKF(trackStateKF);

  std::vector<double> state(6, 0);
  s1Track.getState(state);
  double z = state[0];
  uint nSteps = std::abs(std::floor((nextZ - z) / fDz));

  TMatrixT<double> covMatrix(5, 5);
  s1Track.getCovMatrix(covMatrix);

  if (fDebugLevel > 0) {
    printf("LOG(INFO): integrating with RK4...\n");
  }

  for (uint step = 0; step <= nSteps; step++) {
    // check if inside a material
    checkMaterial(z);

    // track propagation
    integrateMagRK4(trackStateKF, nextStateKF);

    if (fDebugLevel > 1) {
      if (step % 1000 == 0)
        printf("LOG(INFO): STTrackFitter::extrapolateTrackToZ: pre ext.: %.3f %.3f %.3f %.3f %.3f %.3f\n",
               z, trackStateKF[0], trackStateKF[1], trackStateKF[2], trackStateKF[3], trackStateKF[4]);
    }

    // cov matrix propagation
    propagateCovMatrix(trackStateKF, covMatrix);

    // correction for energy loss
    // todo : account for Bremsstrahlung (?)
    if (fCalcLoss && fInMaterial) {
      // correction for qp
      double tx = nextStateKF[2];
      double ty = nextStateKF[3];
      double t = std::sqrt(1. + tx * tx + ty * ty);
      double dl = fDz * t; // signed dLength
      double qp = trackStateKF[4];
      double momLoss = momentumLoss(std::abs(fCharge / qp), dl);
      fSumMomLoss += momLoss;
      if (fDebugLevel > 2) {
        printf("LOG(INFO): STTrackFitter::extrapolateTrackToZ: momLoss = %.6f\n", momLoss);
      }
      nextStateKF[4] = fCharge / (std::abs(fCharge / qp) - momLoss);

      // correction for cov matrix
      // from cbm code
      double p = std::abs(fCharge / qp);
      double effDL = dl / 10. / fMatRadL;
      double s0 = (effDL > std::exp(-1. / 0.038) ) ? qp * 0.0136 * (1. + 0.038 * std::log(effDL)) : 0.;
      double a  = (1. + fMass * fMass * qp * qp) * s0 * s0 * t * t * effDL;
      double Q5 = a * (1. + tx * tx);
      double Q8 = a * tx * ty;
      double Q9 = a * (1. + ty * ty);
      // double eloss = std::sqrt(p * p + fMass * fMass);
      double corr = (1. - std::abs(qp) * effDL * fDEdx);
      double Ecor = corr > 1e-3 * std::abs(qp) ? 1. / corr : 1e3 * p;
      covMatrix(2, 2) += Q5;
      covMatrix(3, 2) += Q8;
      covMatrix(3, 3) += Q9;
      covMatrix(4, 0) *= Ecor;
      covMatrix(4, 1) *= Ecor;
      covMatrix(4, 2) *= Ecor;
      covMatrix(4, 3) *= Ecor;
      covMatrix(4, 4) *= Ecor * Ecor;
    }

    for (uint i = 0; i < 5; i++) {
      trackStateKF[i] = nextStateKF[i];
    }

    if (fDebugLevel > 1) {
      if (step % 1000 == 0) {
        printf("LOG(INFO): STTrackFitter::extrapolateTrackToZ: %.3f %.3f %.3f %.3f %.3f %.3f\n",
               z, trackStateKF[0], trackStateKF[1], trackStateKF[2], trackStateKF[3], trackStateKF[4]);
      }
    }

    fTrajectoryZXrk4->SetPoint(fTrajectoryZXrk4->GetN(), z, trackStateKF[0]);
    z += fDz;
  }

  // update track state
  s1Track.setStateFromKF(z, trackStateKF);
  s1Track.setCovMatrix(covMatrix);
}

void STTrackFitter::propagateCovMatrix(const std::vector<double>& stateKF,
                                       TMatrixT<double>& covMatrix)
{
  double a = 2.99792458 * 1e-1; // [MeV/c] * [Tesla^-1] * [mm^-1]

  // parameters for derivatives
  // see calculation of the jacobian d(x_k)/d(x_k-1) in arXiv: 0511177v1

  double tx = stateKF[2];
  double ty = stateKF[3];
  double bx = fMagFieldHomogeneous[0];
  double by = fMagFieldHomogeneous[1];
  double bz = fMagFieldHomogeneous[2];
  double Axf = std::sqrt(1. + tx * tx + ty * ty) * (ty * (tx * bx + bz) - (1. + tx * tx) * by);

  TMatrixT<double> extJac(5, 5);
  for (uint i = 0; i < 5; i++) {
    extJac(i, i) = 1.;
  }
  extJac(0, 2) = fDz;
  extJac(1, 3) = fDz;
  extJac(0, 4) = 0.5 * fDz * fDz * Axf * a;
  extJac(2, 4) = fDz * Axf * a;

  TMatrixT<double> extJacT = extJac;
  extJacT.T();

  // todo : implement noise (Bremsstrahlung...)
  // TMatrixT<double> noiseMatrix(5, 5);

  // calculating new covariance matrix
  TMatrixT<double> nCovMatrixTmp1 = covMatrix * extJacT;
  TMatrixT<double> nCovMatrixTmp2 = extJac * nCovMatrixTmp1;
  TMatrixT<double> nCovMatrix = nCovMatrixTmp2; // + noiseMatrix;

  // update covariance matrix
  for (uint i = 0; i < 5; i++) {
    for (uint j = 0; j < 5; j++) {
      covMatrix(i, j) = nCovMatrix(i, j);
    }
  }
}

/**
 * Fit track using KF
 * input vectors must be of the same size
 * @param inS1Track : track with initial parameters
 * @param vZcoords : coordinates to extrapolate to
 * @param vMeas : 2d vector of measurements
 * @param covMatrixMeas : covariance matrix for measurements
 * @param filteredStates : filtered states at each measurement
 * @param outS1Track : fitted track
 */
void STTrackFitter::fitTrack(STTrack& inS1Track,
                             const std::vector<double>& vZcoords,
                             const std::vector<std::vector<double>>& vMeas,
                             const TMatrixT<double>& covMatrixMeas,
                             std::vector<std::vector<double>>& filteredStates,
                             STTrack& outS1Track,
                             bool isRefit)
{
  fSumMomLoss = 0;

  if (fDebugLevel > 0) {
    printf("LOG(INFO): STTrackFitter::fitTrack(): in\n");
  }

  std::vector<double> inState(5, 0);
  inS1Track.getStateKF(inState);

  TMatrixT<double> covMatrix(5, 5);
  inS1Track.getCovMatrix(covMatrix);

  if (fDebugLevel > 0) {
    printf("\nLOG(INFO): STTrackFitter::fitTrack(): input covMatrix:\n");
    covMatrix.Print();
  }

  // unity matrix
  TMatrixT<double> I(5, 5);
  for (uint i = 0; i < 5; i++) {
    for (uint j = 0; j < 5; j++) {
      I(i, j) = i == j ? 1. : 0;
    }
  }

  // initializing state
  if (!isRefit) {
    double inTx = vMeas[0][0] / vZcoords[0];
    double inTy = vMeas[0][1] / vZcoords[0];
    outS1Track.setStateFromKF(inS1Track.getZ(), std::vector<double>{inState[0], inState[1], inTx, inTy, inState[4]});
  } else {
    outS1Track.setStateFromKF(inS1Track.getZ(), std::vector<double>{inState[0], inState[1], inState[2], inState[3], inState[4]});
  }
  outS1Track.setCovMatrix(covMatrix);
  for (uint iMeas = 0; iMeas < vMeas.size(); iMeas++) {
    // propagating state and covariance matrix to z
    extrapolateTrackToZ(outS1Track, vZcoords[iMeas]);

    // get extrapolated parameters
    std::vector<double> extState(5, 0);
    outS1Track.getStateKF(extState);
    outS1Track.getCovMatrix(covMatrix);

    if (fDebugLevel > 1) {
      printf("\nLOG(INFO): STTrackFitter::fitTrack(): covMatrix:\n");
      covMatrix.Print();
    }

    if (fDebugLevel > 1) {
      printf("LOG(INFO): STTrackFitter::fitTrack(): extrapolated state: vZcoords[%d] = %.3f, extState[x] = %.3f, extState[y] = %.3f, extState[tx] = %.3f, extState[ty] = %.3f, extState[q/p] = %.3f\n",
             iMeas, vZcoords[iMeas], extState[0], extState[1], extState[2], extState[3], extState[4]);
    }

    // calculating model of measurement H
    TMatrixT<double> measModelH(2, 5);
    measModelH(0, 0) = 1.;
    measModelH(1, 1) = 1.;
    TMatrixT<double> measModelHT = measModelH;
    measModelHT.T();

    // calculating residual vector zeta
    TMatrixT<double> residVecZ(2, 1);
    TMatrixT<double> extStateTMat(5, 1);
    for (uint i = 0; i < 5; i++) {
      extStateTMat(i, 0) = extState[i];
    }

    TMatrixT<double> Hr = measModelH * extStateTMat;
    for (uint i = 0; i < 2; i++) {
      residVecZ(i, 0) = vMeas[iMeas][i] - Hr(i, 0);
    }

    // calculating gain matrix K
    // K = C * H^T * (V + H * C * H^T)^-1
    TMatrixT<double> CHt = covMatrix * measModelHT;
    TMatrixT<double> HCHt = measModelH * CHt;
    TMatrixT<double> VHCHt = covMatrixMeas + HCHt;
    VHCHt.Invert();
    TMatrixT<double> preK = measModelHT * VHCHt;
    TMatrixT<double> gainMatrixK = covMatrix * preK;

    // filtering state vector
    TMatrixT<double> Kz(5, 1);
    Kz = gainMatrixK * residVecZ;
    std::vector<double> outState(5, 0);
    for (uint i = 0; i < 5; i++) {
      outState[i] = extState[i] + Kz(i, 0);
    }

    if (fDebugLevel > 1) {
      printf("LOG(INFO): STTrackFitter::fitTrack(): updated state: vZcoords[%d] = %.3f, outState[x] = %.3f, outState[y] = %.3f, outState[tx] = %.3f, outState[ty] = %.3f, outState[q/p] = %.3f\n",
             iMeas, vZcoords[iMeas], outState[0], outState[1], outState[2], outState[3], outState[4]);
    }

    // updating state vector
    outS1Track.setStateFromKF(vZcoords[iMeas], outState);

    // saving intermediate states
    filteredStates.push_back(outState);

    // filtering covariance matrix
    TMatrixT<double> KH = gainMatrixK * measModelH;
    TMatrixT<double> filtCovMatrixTmp1 = I - KH;
    TMatrixT<double> filtCovMatrix = filtCovMatrixTmp1 * covMatrix;

    // update covariance matrix
    outS1Track.setCovMatrix(filtCovMatrix);

    TMatrixT<double> residVecZT = residVecZ;
    residVecZT.T();

    // updating chi2
    // chi2 = chi2_prev + zeta^T * (V + H * C * H^T)^-1 * zeta
    TMatrixT<double> chi2tmp1 = VHCHt * residVecZ;
    TMatrixT<double> chi2 = residVecZT * chi2tmp1;
    outS1Track.setChi2(outS1Track.getChi2() + chi2(0, 0));
  }
}

void STTrackFitter::clearTrajectory()
{
  delete fTrajectoryZXrk4;
  fTrajectoryZXrk4 = new TGraph();
  fTrajectoryZXrk4->SetMarkerStyle(21);
  fTrajectoryZXrk4->SetMarkerSize(0.5);
  fTrajectoryZXrk4->SetLineWidth(0);
  fTrajectoryZXrk4->SetMarkerColor(kBlue);
}
