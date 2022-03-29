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
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", 100, -100, 100);
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
            << "**********" << '\n';

  TCanvas *cRMS = new TCanvas("cRMS");
  cRMS->cd();
  stdHisto->Draw();
}

void analyse()
{
  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Sinusoide; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  TF1 *fitFunc = new TF1("fitFunc", "[0] * x * sin([1] * x)", -50, 940);
  fitFunc->SetParameters(400, 200, 1);
  fitFunc->SetLineColor(kRed);
  graph->Fit("fitFunc", "R");

  TCanvas *canvas = new TCanvas("canvas");
  graph->Draw("APE");
  // graph->Fit(fitFunc);
}