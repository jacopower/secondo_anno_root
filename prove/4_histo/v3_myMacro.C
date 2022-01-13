//ERRORI IN COMPILAZIONE!!!


#include <iostream>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"

void SetStyle()
{
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  
}

Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void myMacro(Int_t nGen = 1E5)
{
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  char *histName = new char[10];
  TH1F *h[2];

  for (Int_t i = 0; i < 2; ++i)
  {
    sprintf(histName, "h%d", i);
    h[i] = new TH1F(histName, "test histogram", 100, -5., 5.);

    h[i]->SetMarkerStyle(20);
    h[i]->SetMarkerSize(0.5);
    h[i]->SetLineColor(1);
    h[i]->GetYaxis()->SetTitleOffset(1.2);
    h[i]->GetXaxis()->SetTitleSize(0.04);
    h[i]->GetYaxis()->SetTitleSize(0.04);
    h[i]->GetXaxis()->SetTitle("x");
    h[i]->GetYaxis()->SetTitle("Entries");
  }

  // filling histogram with predefined gaussian function
  h[0]->SetFillColor(4);
  h[1]->SetFillColor(2);

  h[0]->FillRandom("gaus", nGen); //gaus predefined function (G(0,1))

  TF1 *f = new TF1("f", myFunction, -10, 10, 3); //user defined function
  f->SetParameter(0, 1);
  f->SetParameter(1, 0);
  f->SetParameter(2, 1);
  h[1]->FillRandom("f", nGen);

  //Drawing
  TCanvas *c1 = new TCanvas("c1", "histo examples", 200, 10, 600, 400);
  c1->Divide(1, 2);

  for (Int_t i = 0; i < 2; ++i)
  {
    c1->cd(i + 1);
    h[i]->Draw("H");
    h[i]->Draw("E, P, SAME");
    //c1->Print("testHisto.gif");
    //c1->Print("testHisto.C");
    //c1->Print("testHisto.root");
  }

  TCanvas *c2 = new TCanvas("c2", "histo examples", 200, 10, 600, 400);
  h[0]->Draw("H");
  h[0]->Draw("E, P, SAME");

  TCanvas *c3 = new TCanvas("c3", "histo examples", 200, 10, 600, 400);
  h[1]->Draw("H");
  h[1]->Draw("E, P, SAME");
}