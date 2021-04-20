#include "STRunAction.hh"
#include "STDigi.hh" // To read digitization parameters
#include "G4Run.hh"

#include "TList.h"
#include "TH2D.h"
#include "TFile.h"

STRunAction::STRunAction()
  : G4UserRunAction(), fList()
{
  hXResiduals = new TH1F("hXResiduals", ";x_{reco} - x_{mc} [mm]; counts", 10000, -50, 50);
  hYResiduals = new TH1F("hYResiduals", ";y_{reco} - y_{mc} [mm]; counts", 10000, -50, 50);
  hTxResiduals = new TH1F("hTxResiduals", ";tx_{reco} - tx_{mc}; counts", 20000, -1, 1);
  hTyResiduals = new TH1F("hTyResiduals", ";ty_{reco} - ty_{mc}; counts", 20000, -1, 1);
  hXPulls = new TH1F("hXPulls", ";pull(x_{reco} - x_{mc}); counts", 20000, -20, 20);
  hYPulls = new TH1F("hYPulls", ";pull(y_{reco} - y_{mc}); counts", 20000, -20, 20);
  hTxPulls = new TH1F("hTxPulls", ";pull(tx_{reco} - tx_{mc}); counts", 20000, -20, 20);
  hTyPulls = new TH1F("hTyPulls", ";pull(ty_{reco} - ty_{mc}); counts", 20000, -20, 20);
  hPResiduals = new TH1F("hPResiduals", ";q/p_{reco} - q/p_{mc} [1/MeV]; counts", 100000, -0.1f, 0.1f);
  hPPulls = new TH1F("hPPulls", ";pull(q/p_{reco} - q/p_{mc}); counts", 20000, -20, 20);
  hPSigmas = new TH1F("hPSigmas", ";#sigma_{q/p}; counts", 100000, -0.05f, 0.05f);
  hChi2s = new TH1F("hChi2s", ";#chi^{2}; counts", 60000, -10, 50);
}

STRunAction::~STRunAction()
{
  delete fList;
  delete hXResiduals;
  delete hYResiduals;
  delete hTxResiduals;
  delete hTyResiduals;
  delete hPResiduals;
  delete hXPulls;
  delete hYPulls;
  delete hTxPulls;
  delete hTyPulls;
  delete hPPulls;
  delete hPSigmas;
  delete hChi2s;
}

void STRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  fList = new TList();
  TH2D* hDigis4 = new TH2D("hDigis4", "Layer 4; x (mm); y (mm)", gnRowsX, -500., 500., gnRowsY, -500., 500.);
  TH1D* hP = new TH1D("hP", ";p(MeV);Entries", 100, 0, 2000.);
  fList->Add(hDigis4);
  fList->Add(hP);
}

void STRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  TFile resultsFile("results.root", "recreate");
  fList->Write();
  resultsFile.Close();

  TFile pullsFile("pulls.root", "recreate");
  hXResiduals->Write();
  hYResiduals->Write();
  hTxResiduals->Write();
  hTyResiduals->Write();
  hTxPulls->Write();
  hTyPulls->Write();
  hPResiduals->Write();
  hXPulls->Write();
  hYPulls->Write();
  hPPulls->Write();
  hPSigmas->Write();
  hChi2s->Write();
  pullsFile.Close();
}
