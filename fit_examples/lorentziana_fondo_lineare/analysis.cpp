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

Double_t myLorentz(Double_t *x, Double_t *par)
{
  // par[0] = centro, par[1] = semiampiezza
  Double_t xx = x[0];
  Double_t val = (1 / TMath::Pi()) * par[1] / ( (xx - par[0]) * (xx - par[0]) + par[1] * par[1]  );
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
  graph->SetTitle("Test Lorentz; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  TF1 *fitFunc = new TF1("fitFunc", myLorentz, -50, -40, 2);
  fitFunc->SetParameters(-45, 1);
  fitFunc->SetLineColor(kRed);
  graph->Fit(fitFunc, "R");

  TCanvas *canvas = new TCanvas("canvas");
  canvas->cd();
  graph->Draw("APE");
}