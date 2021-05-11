#include "STEventAction.hh"
#include "STActionInitialization.hh"
#include "STEventAction.hh"
#include "STPrimaryGeneratorAction.hh"
#include "STRunAction.hh"
#include "STSteppingAction.hh"

STActionInitialization::STActionInitialization()
  : G4VUserActionInitialization()
{
}

STActionInitialization::~STActionInitialization()
{
}

void STActionInitialization::BuildForMaster() const
{
  auto* runAction = new STRunAction();
  SetUserAction(runAction);
}

void STActionInitialization::Build() const
{
  SetUserAction(new STPrimaryGeneratorAction());

  auto* runAction = new STRunAction();
  SetUserAction(runAction);

  auto* eventAction = new STEventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new STSteppingAction(eventAction));
}
