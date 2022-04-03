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

void computeRMS()
{
  constexpr int N = 100;
  // ***** CALCOLO DEVIAZIONE STANDARD DEL RUMORE *****
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", N, -0.5, 0.5);
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

double background(double *x, double *par)
{
  return par[0] + par[1] * x[0];
}
/*
double lorentzianPeak(double *x, double *par)
{
  return (0.5 * par[0] * par[1] / TMath::Pi()) /
         TMath::Max(1.e-10, (x[0] - par[2]) * (x[0] - par[2]) + .25 * par[1] * par[1]);
}
*/
double lorentzianPeak(double *x, double *par)
{
  return (1 / TMath::Pi()) * par[1] / ((x[0] - par[0]) * (x[0] - par[0]) + par[1] * par[1]);
}

double fitFunction(double *x, double *par)
{
  return background(x, par) + lorentzianPeak(x, &par[2]);
}

void analyse()
{
  TCanvas *c1 = new TCanvas("c1", "Fitting Demo", 10, 10, 700, 500);
  c1->SetFillColor(33);
  c1->SetFrameFillColor(41);
  c1->SetGrid();

  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Lorentz; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  //TF1 *fitFunc = new TF1("fitFunc", background, -50, -40, 2);

  TF1 *fitFunc = new TF1("fitFunc", fitFunction, -50, -40, 5);
   fitFunc->SetNpx(500);
  fitFunc->SetLineWidth(4);
  fitFunc->SetLineColor(kMagenta);
  fitFunc->SetParameter(0, -30);
  fitFunc->SetParameter(1, 1.2);
  //fitFunc->SetParameter(2, -45);
  fitFunc->SetParameter(3, 0.2);

  //fitFunc->SetParameter(0, -45);
  //fitFunc->SetParameter(2, -45);
  //fitFunc->SetParameter(3, 0.2);


  graph->Fit("fitFunc");

   graph->Draw("APE");
  //fitFunc->Draw();
}
