#ifndef S1PrimaryGeneratorAction_h
#define S1PrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4ParticleGun;
class G4Event;

class STPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  STPrimaryGeneratorAction();
  virtual ~STPrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event*);

 private:
  G4ParticleGun* fParticleGun;
};
#endif
