// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// 1. Si definiscono 3 istogrammi monodimensionali di 100 bin in un range da 0 a 5.
// 2. Si riempie il primo istogramma con 10^4 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) e distribuite secondo una distribuzione esponenziale decrescente con media =0.5
// 3. Si riempie il secondo istogramma con 10^6 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. atHisfortraverso gRandom) e distribuite secondo una distribuzione gaussiana con media =2.5 e deviazione standard =0.25.
// 4. Si riempie il terzo istogramma con 10^5 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) uniforme in [0,5]
// 5. Si fa la somma dei tre istogrammi. Si stampano a schermo i parametri corrispondenti alla media e alla deviazione standard della gaussiana che risulta dal fit, con relativo errore, e il chi quadro ridotto.

#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
}

// Gaus + A * exp{-Bx} + C
Double_t myFuncion(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3] * TMath::Exp(-par[4] * xx) + par[5];
  return val;
}

void myMacro()
{
  TH1F *h1 = new TH1F("h1", "histo1", 100, 0., 5.);
  //h1->Sumw2();
  TH1F *h2 = new TH1F("h2", "histo2", 100, 0., 5.);
  //h2->Sumw2();
  TH1F *h3 = new TH1F("h3", "histo3", 100, 0., 5.);
  //h3->Sumw2();

  Double_t x1 = 0;
  for (Int_t i = 0; i < 1E4; ++i)
  {
    x1 = gRandom->Exp(0.5);
    h1->Fill(x1);
  }

  Double_t x2 = 0;
  for (Int_t i = 0; i < 1E6; ++i)
  {
    x2 = gRandom->Gaus(2.5, 0.25);
    h2->Fill(x2);
  }

  Double_t x3 = 0;
  for (Int_t i = 0; i < 1E5; ++i)
  {
    x3 = gRandom->Uniform(0., 5.);
    h3->Fill(x3);
  }

  TH1F *hSum12 = new TH1F(*h1);
  TH1F *hSum = new TH1F(*h1);
  hSum12->Add(h1, h2, 1, 1);
  hSum->Add(hSum12, h3, 1, 1);

  TF1 *f = new TF1("f", myFuncion, 0., 5., 6);
  f->SetParameters(1, 1, 1, 1, 1, 1);
  f->SetLineWidth(1);
  f->SetLineColor(kRed);

  TCanvas *canvas = new TCanvas("canvas");
  hSum->Fit("f");
  hSum->Draw();

  TF1 *fitFunc = hSum->GetFunction("f");
  Double_t media = fitFunc->GetParameter(1);
  Double_t mediaErr = fitFunc->GetParError(1);
  Double_t stdDev = fitFunc->GetParameter(2);
  Double_t stdDevErr = fitFunc->GetParError(2);
  Double_t chisquare = fitFunc->GetChisquare();
  Double_t ndof = fitFunc->GetNDF();

  std::cout << "***** PARAMETRI FIT *****" << '\n';
  std::cout << "Media: " << media << " +/- " << mediaErr << '\n';
  std::cout << "StdDev: " << stdDev << " +/- " << stdDevErr << '\n';
  std::cout << "Chi/NDOF: " << chisquare / ndof << '\n';
  std::cout << "**********" << '\n';
}