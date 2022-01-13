// Scrivi una macro in cui
// 1) Si definiscono 2 istogrammi monodim di 500 bin tra 0 e 5
// 2) Si riempie il primo istogramma con 10^6 occorrenze di una variabile casuale x generata esplicitamente e singolarmente e distribuita
// secondo una gaussiana con media 2.5 e dev std 0.25
// 3) Si riempie il secondo istogramma con 10^5 occorrenze di una variabile casuale x generata esplicitamente e singolarmente e distribuita
// secondo una esponenziale decrescente con media 1
// 4) Si fa la somma dei due istogrammi, e si effettua il fit dell'istogramma somma secondo una forma funzionale consistente di una gaussiana
// (3 parametri: ampiezza, media e stddev) e un esponenziale (2 parametri), per un totale di 5 parametri liberi
// 5) Stampa a scherma i parametri dopo il fit con relativo errore, e i X^2 ridotto

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
}

Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3] * TMath::Exp(par[4] * xx);
  return val;
}

void myMacro()
{
  constexpr Int_t nGen1 = 1E6;
  constexpr Int_t nGen2 = 1E5;

  TH1F *h1 = new TH1F("h1", "histo1", 500, 0., 5.);
  h1->Sumw2();
  TH1F *h2 = new TH1F("h2", "histo2", 500, 0., 5.);
  h2->Sumw2();

  Double_t x1, x2 = 0;

  for (Int_t i = 0; i < nGen1; ++i)
  {
    x1 = gRandom->Gaus(2.5, 0.25);
    h1->Fill(x1);
  }

  for (Int_t i = 0; i < nGen2; ++i)
  {
    x2 = gRandom->Exp(1);
    h2->Fill(x2);
  }

  TH1F *hSum = new TH1F(*h1);
  hSum->Add(h1, h2, 1, 1);

  TF1 *fitFunc = new TF1("fitFunc", myFunction, 0, 5, 5);
  fitFunc->SetParameters(1, 1, 1, 1, 1);
  hSum->Fit("fitFunc");

  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineWidth(2);

  TCanvas *canvas = new TCanvas("canvas");
  canvas->cd();
  hSum->Draw("H");
  hSum->Draw("E, SAME");
}