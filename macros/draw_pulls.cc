#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <string>

int draw_pulls()
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

  std::string filename("pulls");
  TFile* f = new TFile((filename + ".root").c_str());
  f->ls();
  // ----------------------------------------------
  auto* hXResiduals = (TH1D*)f->Get("hXResiduals");
  hXResiduals->Rebin(10);
  hXResiduals->GetXaxis()->SetRangeUser(-2., 2.);
  auto hXRes_fit = hXResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hYResiduals = (TH1D*)f->Get("hYResiduals");
  hYResiduals->Rebin(10);
  hYResiduals->GetXaxis()->SetRangeUser(-2., 2.);
  auto hYRes_fit = hYResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hTxResiduals = (TH1D*)f->Get("hTxResiduals");
  hTxResiduals->Rebin(10);
  hTxResiduals->GetXaxis()->SetRangeUser(-0.04, 0.04);
  auto hTxRes_fit = hTxResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hTyResiduals = (TH1D*)f->Get("hTyResiduals");
  hTyResiduals->Rebin(10);
  hTyResiduals->GetXaxis()->SetRangeUser(-0.04, 0.04);
  auto hTyRes_fit = hTyResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hXPulls = (TH1D*)f->Get("hXPulls");
  hXPulls->Rebin(200);
  hXPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hXPulls_fit = hXPulls->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hYPulls = (TH1D*)f->Get("hYPulls");
  hYPulls->Rebin(200);
  hYPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hYPulls_fit = hYPulls->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hTxPulls = (TH1D*)f->Get("hTxPulls");
  hTxPulls->Rebin(200);
  hTxPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hTxPulls_fit = hTxPulls->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hTyPulls = (TH1D*)f->Get("hTyPulls");
  hTyPulls->Rebin(200);
  hTyPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hTyPulls_fit = hTyPulls->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hPResiduals = (TH1D*)f->Get("hPResiduals");
  hPResiduals->Rebin(16);
  hPResiduals->GetXaxis()->SetRangeUser(-0.001, 0.001);
  auto hPResiduals_fit = hPResiduals->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hPPulls = (TH1D*)f->Get("hPPulls");
  hPPulls->Rebin(200);
  hPPulls->GetXaxis()->SetRangeUser(-12., 12.);
  auto hPPulls_fit = hPPulls->Fit("gaus", "s");
  // ----------------------------------------------
  auto* hChi2s = (TH1D*)f->Get("hChi2s");
  hChi2s->Rebin(1200);
  hChi2s->GetXaxis()->SetRangeUser(0., 50.);

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
  canvas->cd(11);
  hChi2s->Draw("");

  canvas->Print((filename + ".png").c_str());

  return 0;
}
