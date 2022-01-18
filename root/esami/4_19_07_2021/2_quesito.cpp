// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// Si definisce un istogramma unidimensionale di 1000 bin in un range da 0 a 5.
// Si riempie il primo istogramma con 106 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) e distribuite secondo un’esponenziale decrescente con media µ=1.
// Si riempie il secondo istogramma con 107 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) e distribuite secondo una distribuzione gaussiana con media µ=2.5 e deviazione standard σ=0.25.
// Si fa la somma dei due istogrammi, e si effettua il fit dell’istogramma risultante secondo una forma funzionale consistente con una somma di una gaussiana e di un esponenziale.
// Si stampano a schermo la media e la deviazione standard della gaussiana che risultano dal fit, con relativo errore, il χ2 ridotto e la probabilità del fit.

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

// Gaus + A * exp{- B * x}
Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3] * TMath::Exp(-par[4] * xx);
  return val;
}

void myMacro()
{
  TH1F *h1 = new TH1F("h1", "histo1", 1000, 0., 5.);
  TH1F *h2 = new TH1F("h2", "histo2", 1000, 0., 5.);

  Double_t x1 = 0;
  for (Int_t i = 0; i < 1E6; ++i)
  {
    x1 = gRandom->Exp(1);
    h1->Fill(x1);
  }

  Double_t x2 = 0.;
  for (Int_t i = 0; i < 1E7; ++i)
  {
    x2 = gRandom->Gaus(2.5, 0.25);
    h2->Fill(x2);
  }

  TH1F *hSum = new TH1F(*h1);
  hSum->SetTitle("Sum Histo");
  hSum->SetName("hSum");
  hSum->Add(h1, h2, 1, 1);

  TF1 *f = new TF1("f", myFunction, 0., 5., 5);
  f->SetParameters(1, 1, 1, 1, 1);
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
  Double_t chiquadro = fitFunc->GetChisquare();
  Double_t Ndof = fitFunc->GetNDF();
  Double_t prob = fitFunc->GetProb();

  std::cout << "***** PARAMETRI FIT *****" << '\n';
  std::cout << "Media: " << media << " +/- " << mediaErr << '\n';
  std::cout << "StdDev: " << stdDev << " +/- " << stdDevErr << '\n';
  std::cout << "Chi/NDOF: " << chiquadro / Ndof << '\n';
  std::cout << "Prob: " << prob << '\n';
  std::cout << "**********" << '\n';
}