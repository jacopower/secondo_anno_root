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

/*
***** PROBLEMI *****
IL FIT VISIVAMENTE BUONO  HA CHIQUADRO ALTRISSIMO, PROBLEMA CON CALCOLO DELL'ERRORE?
*/

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
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", N, -0.2, 0.2);
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
            << "Overflows: " << stdHisto->GetBinContent(N+1) << '\n'
            << "**********" << '\n';

  TCanvas *cRMS = new TCanvas("cRMS");
  stdHisto->Draw();
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