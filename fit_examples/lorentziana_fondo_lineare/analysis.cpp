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
  Double_t val = (1 / TMath::Pi()) * par[1] / ((xx - par[0]) * (xx - par[0]) + par[1] * par[1]);
  return val;
}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par)
{
  return par[0] + par[1] * x[0];
}

// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par)
{
  return (0.5 * par[0] * par[1] / TMath::Pi()) / TMath::Max(1.e-10,
                                                            (x[0] - par[2]) * (x[0] - par[2]) + .25 * par[1] * par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par)
{
  return background(x, par) + lorentzianPeak(x, &par[3]);
}

void computeRMS()
{
  constexpr int N = 100;
  // ***** CALCOLO DEVIAZIONE STANDARD DEL RUMORE *****
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", N, -1, 1);
  std::ifstream in;
  in.open("rumore.txt");
  Float_t rumore;
  while (1)
  {
    in >> rumore;
    if (!in.good())
    {
      break;
    }
    stdHisto->Fill(rumore);
  }
  in.close();

  Double_t stdDev = stdHisto->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE *****" << '\n'
            << "Deviazione Standard: " << stdDev << '\n'
            << "Underflows: " << stdHisto->GetBinContent(0) << '\n'
            << "Overflows: " << stdHisto->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  TCanvas *cRMS = new TCanvas("cRMS");
  stdHisto->Draw();
}

void analyse()
{
  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Lorentz; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);
  /*
    TF1 *fitFunc = new TF1("fitFunc", myLorentz, -50, -40, 2);
    fitFunc->SetParameters(-45, 1);
    fitFunc->SetLineColor(kRed);
    graph->Fit(fitFunc, "R");
    */

  // create a TF1 with the range from 0 to 3 and 6 parameters
  TF1 *fitFcn = new TF1("fitFcn", fitFunction, -50, -40, 5);

  // first try without starting values for the parameters
  // this defaults to 1 for each param.
  graph->Fit("fitFcn", "R");
  // this results in an ok fit for the polynomial function however
  // the non-linear part (Lorentzian

  TCanvas *canvas = new TCanvas("canvas");
  canvas->cd();
  graph->Draw("APE");
}