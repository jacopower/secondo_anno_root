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

// ***** USER-DEFINED GAUSS *****
Double_t myGaus(Double_t *x, Double_t *par)
{
  // par[0] = h, par[1] = x0, par[2] = sigma
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void myMacro()
{
  constexpr int Npoints = 100;

  // ***** CALCOLO DEVIAZIONE STANDARD DEL RUMORE *****
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", 100, -0.2, 0.2);
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

  // ***** CREO GRAFICO DEI PUNTI SPERIMENTALI E FACCIO FIT *****
  // VA AGGIUNTO ERRORE ATTRAVERSO LA VARIABILE stdDev
  TGraphErrors *graph = new TGraphErrors("data.txt", "%lg %lg %lg");
  graph->SetTitle("Test Gaussiana; x(UDM); y(UDM)");
  graph->SetMarkerStyle(kOpenCircle);
  graph->SetMarkerColor(kBlue);
  graph->SetFillColor(0);

  TF1 *fitFunc = new TF1("fitFunc", myGaus, -10, 10, 3);
  fitFunc->SetParameters(1, 1, 1); // settati a caso
  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineStyle(2);
  graph->Fit(fitFunc);

  // ***** CREO GRAFICO DEL RESIDUO *****
  in.open("data.txt");
  Float_t x, y, err;
  Double_t residuo = 0;
  Float_t xRes[Npoints], yRes[Npoints];

  for (Int_t i = 0; i < Npoints; ++i)
  {
    in >> x >> y >> err;
    residuo = y - (fitFunc->Eval(x));
    xRes[i] = x;
    yRes[i] = residuo;
  }
  in.close();
  TGraph *residueGraph = new TGraph(Npoints, xRes, yRes);

  // ***** LEGENDA E CANVAS *****
  TLegend *leg = new TLegend(.1, .7, .3, .9, "Dati Laboratorio");
  leg->SetFillColor(0);
  leg->AddEntry(graph, "Punti Sperimentali");
  leg->AddEntry(fitFunc, "Fit Gaussiano");
  leg->AddEntry(residueGraph, "Residuo");

  TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  // ***** PAD 1 - UPPER PLOT *****
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  graph->Draw();
  leg->Draw("Same");

  // ***** PAD 2 - LOWER PLOT *****
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd(); // pad2 becomes the current pad
  residueGraph->Draw();
}
