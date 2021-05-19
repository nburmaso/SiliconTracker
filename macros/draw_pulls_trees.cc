#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include <string>
#include <cmath>

int draw_pulls_trees()
{
  gStyle->SetLineScalePS(2.3);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTitleOffset(1.25, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextSize(0.06);
  gStyle->SetTextFont(42);

  auto* pdg = new TDatabasePDG();

  TFile* fMCEvents = TFile::Open("mc_tracks.root", "r");
  auto* tParticles = (TTree*)fMCEvents->Get("mc_tracks");
  int mcEventID;
  int mcTrackID;
  int pdgID;
  int motherID;
  double mcvx;
  double mcvy;
  double mcvz;
  double mcpx;
  double mcpy;
  double mcpz;
  tParticles->SetBranchAddress("mcEventID", &mcEventID);
  tParticles->SetBranchAddress("mcTrackID", &mcTrackID);
  tParticles->SetBranchAddress("pdgID", &pdgID);
  tParticles->SetBranchAddress("motherID", &motherID);
  tParticles->SetBranchAddress("vx", &mcvx);
  tParticles->SetBranchAddress("vy", &mcvy);
  tParticles->SetBranchAddress("vz", &mcvz);
  tParticles->SetBranchAddress("px", &mcpx);
  tParticles->SetBranchAddress("py", &mcpy);
  tParticles->SetBranchAddress("pz", &mcpz);

  // ----------------------------------------------

  TFile* fEvents = TFile::Open("tracks.root", "r");
  auto* tTracks = (TTree*)fEvents->Get("reco_tracks");
  int eventID;
  int trackID;
  int mcLabel;
  double x;
  double y;
  double z;
  double tx;
  double ty;
  double tz;
  double qp;
  double covM[15];
  tTracks->SetBranchAddress("eventID", &eventID);
  tTracks->SetBranchAddress("trackID", &trackID);
  tTracks->SetBranchAddress("mcLabel", &mcLabel);
  tTracks->SetBranchAddress("z", &z);
  tTracks->SetBranchAddress("x", &x);
  tTracks->SetBranchAddress("y", &y);
  tTracks->SetBranchAddress("tx", &tx);
  tTracks->SetBranchAddress("ty", &ty);
  tTracks->SetBranchAddress("qp", &qp);
  tTracks->SetBranchAddress("covM", &covM);

  // ----------------------------------------------

  auto* hXResiduals = new TH1F("hXResiduals", ";x_{reco} - x_{mc} [mm]; counts", 10000, -50, 50);
  auto* hYResiduals = new TH1F("hYResiduals", ";y_{reco} - y_{mc} [mm]; counts", 10000, -50, 50);
  auto* hTxResiduals = new TH1F("hTxResiduals", ";tx_{reco} - tx_{mc}; counts", 20000, -1, 1);
  auto* hTyResiduals = new TH1F("hTyResiduals", ";ty_{reco} - ty_{mc}; counts", 20000, -1, 1);
  auto* hXPulls = new TH1F("hXPulls", ";pull(x_{reco} - x_{mc}); counts", 20000, -20, 20);
  auto* hYPulls = new TH1F("hYPulls", ";pull(y_{reco} - y_{mc}); counts", 20000, -20, 20);
  auto* hTxPulls = new TH1F("hTxPulls", ";pull(tx_{reco} - tx_{mc}); counts", 20000, -20, 20);
  auto* hTyPulls = new TH1F("hTyPulls", ";pull(ty_{reco} - ty_{mc}); counts", 20000, -20, 20);
  auto* hPResiduals = new TH1F("hPResiduals", ";q/p_{reco} - q/p_{mc} [1/MeV]; counts", 100000, -0.1f, 0.1f);
  auto* hPPulls = new TH1F("hPPulls", ";pull(q/p_{reco} - q/p_{mc}); counts", 20000, -20, 20);

  // ----------------------------------------------

  for (int iTrack = 0; iTrack < tTracks->GetEntries(); iTrack++) {
    tTracks->GetEntry(iTrack);
    for (int iParticle = 0; iParticle < tParticles->GetEntries(); iParticle++) {
      tParticles->GetEntry(iParticle);
      if (eventID == mcEventID && mcLabel == mcTrackID && motherID == 0) {
        auto particle = pdg->GetParticle(pdgID);
        double charge = particle->Charge() / 3.; // root???
        double mcp = std::sqrt(mcpx * mcpx + mcpy * mcpy + mcpz * mcpz);
        hXResiduals->Fill(x - mcvx);
        hYResiduals->Fill(y - mcvy);
        hTxResiduals->Fill(tx); // assuming initial direction to be (0,0,0)
        hTyResiduals->Fill(ty); // assuming initial direction to be (0,0,0)
        hPResiduals->Fill(qp - charge / mcp);
        hXPulls->Fill((x - mcvx) / std::sqrt(covM[0]));
        hYPulls->Fill((y - mcvy) / std::sqrt(covM[2]));
        hTxPulls->Fill((tx) / std::sqrt(covM[5]));
        hTyPulls->Fill((ty) / std::sqrt(covM[9]));
        hPPulls->Fill((qp - charge / mcp) / std::sqrt(covM[14]));
      }
    }
  }

  // ----------------------------------------------
  hXResiduals->Rebin(10);
  hXResiduals->GetXaxis()->SetRangeUser(-2., 2.);
  auto hXRes_fit = hXResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  hYResiduals->Rebin(10);
  hYResiduals->GetXaxis()->SetRangeUser(-2., 2.);
  auto hYRes_fit = hYResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  hTxResiduals->Rebin(10);
  hTxResiduals->GetXaxis()->SetRangeUser(-0.04, 0.04);
  auto hTxRes_fit = hTxResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  hTyResiduals->Rebin(10);
  hTyResiduals->GetXaxis()->SetRangeUser(-0.04, 0.04);
  auto hTyRes_fit = hTyResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  hXPulls->Rebin(200);
  hXPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hXPulls_fit = hXPulls->Fit("gaus", "s");
  // ----------------------------------------------
  hYPulls->Rebin(200);
  hYPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hYPulls_fit = hYPulls->Fit("gaus", "s");
  // ----------------------------------------------
  hTxPulls->Rebin(200);
  hTxPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hTxPulls_fit = hTxPulls->Fit("gaus", "s");
  // ----------------------------------------------
  hTyPulls->Rebin(200);
  hTyPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hTyPulls_fit = hTyPulls->Fit("gaus", "s");
  // ----------------------------------------------
  hPResiduals->Rebin(16);
  hPResiduals->GetXaxis()->SetRangeUser(-0.001, 0.001);
  auto hPResiduals_fit = hPResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  hPPulls->Rebin(200);
  hPPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hPPulls_fit = hPPulls->Fit("gaus", "s");

  auto* canvas = new TCanvas("canvas", "pulls", 1920, 1080);
  canvas->Divide(4, 3);

  canvas->cd(1);
  hXResiduals->Draw("");
  hXRes_fit->Draw("same");

  canvas->cd(2);
  hYResiduals->Draw("");
  hYRes_fit->Draw("same");

  canvas->cd(3);
  hTxResiduals->Draw("");
  hTxRes_fit->Draw("same");

  canvas->cd(4);
  hTyResiduals->Draw("");
  hTyRes_fit->Draw("same");

  canvas->cd(5);
  hXPulls->Draw("");
  hXPulls_fit->Draw("same");

  canvas->cd(6);
  hYPulls->Draw("");
  hYPulls_fit->Draw("same");

  canvas->cd(7);
  hTxPulls->Draw("");
  hTxPulls_fit->Draw("same");

  canvas->cd(8);
  hTyPulls->Draw("");
  hTyPulls_fit->Draw("same");

  canvas->cd(9);
  hPResiduals->Draw("");
  hPResiduals_fit->Draw("same");

  canvas->cd(10);
  hPPulls->Draw("");
  hPPulls_fit->Draw("same");

  canvas->Print("pulls.png");

  return 0;
}
