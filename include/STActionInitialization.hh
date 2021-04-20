#ifndef S1ActionInitialization_h
#define S1ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class STActionInitialization : public G4VUserActionInitialization
{
 public:
  STActionInitialization();
  virtual ~STActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
};
#endif
