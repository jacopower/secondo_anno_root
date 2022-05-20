#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>

void SetStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

void myMacro(Int_t nGen = 1E5)
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

  h->FillRandom("gaus", nGen);

  TCanvas *c1 = new TCanvas("c1", "test histo with predef function", 200, 10, 600, 400);
  c1->cd();
  h->Draw("H");
  h->Draw("E, P, SAME");

  //c1->Print("testHisto.gif");
  //c1->Print("testHisto.C");
  //c1->Print("testHisto.root");

  Double_t hEntries = h->GetEntries();
  Double_t hMean = h->GetMean();
  Double_t hRMS = h->GetRMS();
  Double_t hMeanErr = h->GetMeanError();
  Double_t hRMSErr = h->GetRMSError();

  std::cout << "Total number of entries: " << hEntries << '\n';
  std::cout << "Mean: " << hMean << '\n';
  std::cout << "RMS: " << hRMS << '\n';
  std::cout << "Location of Maximum: " << h->GetBinCenter(h->GetMaximumBin()) << '\n';
  std::cout << "Maximum of the Distribution: " << h->GetMaximum() << '\n';
}