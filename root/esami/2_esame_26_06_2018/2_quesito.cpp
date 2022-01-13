// Scrivi una Macro in cui
// 1) Si definiscono 2 istogrammi monodimensionali di 500 bin in un range da 0 a 5.
// 2) Si riempe il primo istogramma con 10^6 occorrenze generate esplicitamente e singolarmente secondo una distribuzione gaussiana con
// media 2 e deviazione standard 0.5.
// 3) Si riempe il secondo istogramma con 105 occorrenze generate esplicitamente e singolarmente secondo una distribuzione
// esponenziale decrescente con media 1
// 4) Si fa la somma dei due istogrammi, e si effettua il Fit dell’istogramma somma secondo una forma funzionale consistente di
// una gaussiana (3 parametri) e un esponenziale (2 parametri), per un totale di 5 parametri liberi.
// 5) Si stampa a schermo il valore dei parametri dopo il fit, con relativo errore, e il ChiQuadro ridotto

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

// Gaus + A * exp{-B * x}
Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3] * TMath::Exp(-par[4] * xx);
  return val;
}

void myMacro()
{
  constexpr Int_t nGen1 = 1E6;
  constexpr Int_t nGen2 = 1E5;

  TH1F *h1 = new TH1F("h1", "histo1", 500, 0., 5.);
  TH1F *h2 = new TH1F("h2", "histo2", 500, 0., 5.);

  Double_t x1 = 0;
  for (Int_t i = 0; i < nGen1; ++i)
  {
    x1 = gRandom->Gaus(2., 0.5);
    h1->Fill(x1);
  }

  Double_t x2 = 0;
  for (Int_t i = 0; i < nGen2; ++i)
  {
    x2 = gRandom->Exp(1);
    h2->Fill(x2);
  }

  TH1F *hSum = new TH1F(*h1);
  hSum->SetTitle("Istogramma Somma");
  hSum->SetName("hSum");
  hSum->Add(h1, h1, 1, 1);

  TF1 *f = new TF1("f", myFunction, 0., 5., 5);
  f->SetParameters(1, 1, 1, 1, 1);
  f->SetLineColor(kRed);
  f->SetLineWidth(1);

  TCanvas *canvas = new TCanvas("canvas");
  hSum->Fit("f");
  hSum->Draw();

  TCanvas *c = new TCanvas("c");
  c->Divide(1, 2);
  c->cd(1);
  h1->Draw();
  c->cd(2);
  h2->Draw();

  TF1 *fitFunc = hSum->GetFunction("f");
  Double_t ampiezza = fitFunc->GetParameter(0);
  Double_t media = fitFunc->GetParameter(1);
  Double_t stdDev = fitFunc->GetParameter(2);
  Double_t A = fitFunc->GetParameter(3);
  Double_t B = fitFunc->GetParameter(4);
  Double_t ampiezzaErr = fitFunc->GetParError(0);
  Double_t mediaErr = fitFunc->GetParError(1);
  Double_t stdDevErr = fitFunc->GetParError(2);
  Double_t AErr = fitFunc->GetParError(3);
  Double_t BErr = fitFunc->GetParError(4);
  Double_t chiSquare = fitFunc->GetChisquare();
  Double_t NDOF = fitFunc->GetNDF();

  std::cout << "***** PARAMETRI FIT *****" << '\n';
  std::cout << "Ampiezza: " << ampiezza << " +/- " << ampiezzaErr << '\n';
  std::cout << "Media: " << media << " +/- " << mediaErr << '\n';
  std::cout << "StdDev: " << stdDev << " +/- " << stdDevErr << '\n';
  std::cout << "A: " << A << " +/- " << AErr << '\n';
  std::cout << "B: " << B << " +/- " << BErr << '\n';
  std::cout << "Chi/NDOF: " << chiSquare / NDOF << '\n';
  std::cout << "**********" << '\n';
}
