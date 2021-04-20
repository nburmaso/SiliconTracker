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
  STRunAction* runAction = new STRunAction;
  SetUserAction(runAction);
}

void STActionInitialization::Build() const
{
  SetUserAction(new STPrimaryGeneratorAction);

  STRunAction* runAction = new STRunAction;
  SetUserAction(runAction);

  STEventAction* eventAction = new STEventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new STSteppingAction(eventAction));
}
