#include "TH1F.h"
#include "TCanvas.h"

void SetStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

void myMacro()
{
  TH1F *h = new TH1F("h", "histo", 100, -5., 5.);

  h->SetFillColor(kBlue);
  h->SetLineColor(kBlack);
  h->SetMarkerStyle(kFullCircle);
  h->SetMarkerSize(0.8);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleOffset(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitle("x");
  h->GetYaxis()->SetTitle("Entries");

  h->FillRandom("gaus", 1E6);

  TCanvas *c1 = new TCanvas("c1", "test histo with predef function", 200, 10, 600, 400);
  c1->cd();
  h->Draw("H");
  h->Draw("E, P, SAME");

  c1->Print("testHisto.gif");
  c1->Print("testHisto.C");
  c1->Print("testHisto.root");
}