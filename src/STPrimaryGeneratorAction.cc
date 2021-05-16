#include "STPrimaryGeneratorAction.hh"
#include "STSimulationPars.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "TRandom.h"
#include "TRandomGen.h"
#include "TVector3.h"

using namespace sim_pars;

STPrimaryGeneratorAction::STPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(new G4ParticleGun(1))
{
}

STPrimaryGeneratorAction::~STPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void STPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  gRandom = new TRandomMT64();
  gRandom->SetSeed(std::time(nullptr));
  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("mu+"));
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  eParticleMom = 1000; // [MeV]
  // + gRandom->Gaus(0, 100);
  //  double theta = gRandom->Uniform(-0.25, 0.25);
  //  double phi = gRandom->Uniform(-0.25, 0.25);
  double theta = 0.;
  double phi = 0.;
  TVector3 p3;
  p3.SetMagThetaPhi(eParticleMom, theta, phi);
  G4ThreeVector p(p3.x(), p3.y(), p3.z());

  fParticleGun->SetParticleMomentum(p);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  printf("Particle momentum = %.3f\n", p.mag());

  delete gRandom;
}
