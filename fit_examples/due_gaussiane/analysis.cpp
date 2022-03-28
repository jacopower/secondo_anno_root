#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include <fstream>
#include <iostream>

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(1111);
}

Double_t myGaus(Double_t *x, Double_t *par)
{
  // par[0] = h, par[1] = x0, par[2] = sigma
  Double_t xx = x[0];
  Double_t val = (xx > 20 && xx < 25) * par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void computeRMS()
{
  TH1F *rumoreHisto = new TH1F("rumoreHisto", "Rumore", 100, -5, 5);
  std::ifstream inRumore;
  inRumore.open("rumore.txt");
  Float_t rumore;
  while (1)
  {
    inRumore >> rumore;
    if (!inRumore.good())
    {
      break;
    }
    rumoreHisto->Fill(rumore);
  }
  inRumore.close();

  Double_t stdDev = rumoreHisto->GetRMS();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE *****" << '\n'
            << "Deviazione Standard: " << stdDev << '\n'
            << "**********" << '\n';

  TCanvas *rumoreCanvas = new TCanvas("rumoreCanvas");
  rumoreCanvas->cd();
  rumoreHisto->Draw();
}

void analyse()
{
  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Due Gaussiana; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  TF1 *g1 = new TF1("g1", "gaus", 20, 25);
  TF1 *g2 = new TF1("g2", "gaus", 27, 30);

  TF1 *total = new TF1("total", "gaus(0)+gaus(3)", 20, 30);
  total->SetLineColor(kRed);

  graph->Fit(g1, "R");
  graph->Fit(g2, "R+");

  Double_t par[6];
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);
  total->SetParameters(par);
  graph->Fit(total, "R+");

  TCanvas *c = new TCanvas("c");
  graph->Draw("APE");
}