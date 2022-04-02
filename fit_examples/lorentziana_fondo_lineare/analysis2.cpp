#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
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

// Quadratic background function
double background(double *x, double *par)
{
  return par[0] + par[1] * x[0];
}

// Lorenzian Peak function
double lorentzianPeak(double *x, double *par)
{
  return (0.5 * par[0] * par[1] / TMath::Pi()) /
         TMath::Max(1.e-10, (x[0] - par[2]) * (x[0] - par[2]) + .25 * par[1] * par[1]);
}

// Sum of background and peak function
double fitFunction(double *x, double *par)
{
  return background(x, par) + lorentzianPeak(x, &par[2]); // 3!
}

void analyse()
{
  TCanvas *c1 = new TCanvas("c1", "Fitting Demo", 10, 10, 700, 500);

  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Lorentz; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  TF1 *fitFcn = new TF1("fitFcn", fitFunction, -50, -40, 5);
  // fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(4);
  fitFcn->SetLineColor(kMagenta);

  fitFcn->SetParameters(1, 1, 1, 1, 1);
  graph->Fit("fitFcn", "0");

  fitFcn->SetParameter(3, 0.2); // width
  fitFcn->SetParameter(4, 1);   // peak

  graph->Fit("fitFcn", "V+", "ep");
}