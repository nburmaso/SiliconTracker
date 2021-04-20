#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
void draw(){
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  TFile* fResults = new TFile("results.root");
  fResults->ls();
  TH2D* hDigis4 = (TH2D*) fResults->Get("hDigis4");
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  hDigis4->SetTitleOffset(1.5,"Y");
  hDigis4->DrawClone("colz");
  
  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  hDigis4->GetYaxis()->SetRangeUser(-100.,100);
  hDigis4->GetXaxis()->SetRangeUser(   0.,200);
  hDigis4->DrawClone("colz");
}
