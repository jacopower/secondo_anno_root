// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// 1. Si definiscono 2 istogrammi monodimensionali di 200 bin in un range da 0 a 2.
// 2. Si riempe il primo istogramma con 106 occorrenze generate esplicitamente e singolarmente di una gaussiana con media μ=0.89 e deviazione standard σ=0.05.
// 3. Si riempe il secondo istogramma con 105 occorrenze generate esplicitamente e singolarmente e distribuite uniformemente nel range [0,2]
// 4. Si fa la somma dei due istogrammi, e si effettua il Fit dell’istogramma somma secondo una forma funzionale consistente con la somma di una gaussiana e un polinomio di grado 0.
// 5. Si stampa a schermo il valore dei parametri dopo il fit, con relativa incertezza, e il χ2 ridotto

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

// Gaus + A
Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3];
  return val;
}

void myMacro()
{
  TH1F *h1 = new TH1F("h1", "histo1", 200, 0., 2.);
  TH1F *h2 = new TH1F("h2", "histo2", 200, 0., 2.);

  Double_t x1 = 0;
  for (Int_t i = 0; i < 1E6; ++i)
  {
    x1 = gRandom->Gaus(0.89, 0.05);
    h1->Fill(x1);
  }

  Double_t x2 = 0;
  for (Int_t i = 0; i < 1E5; ++i)
  {
    x2 = gRandom->Uniform(0., 2.);
    h2->Fill(x2);
  }

  TH1F *hSum = new TH1F(*h1);
  hSum->SetTitle("Sum Histo");
  hSum->SetName("hSum");
  hSum->Add(h1, h2, 1, 1);

  TF1 *f = new TF1("f", myFunction, 0., 2., 4);
  f->SetParameters(1, 1, 1, 1);
  f->SetLineWidth(1);
  f->SetLineColor(kRed);

  TCanvas *canvas = new TCanvas("canvas");
  hSum->Fit("f");
  hSum->Draw();

  // POI C'E' LA STAMPA MA NON HO VOGLIA
}
