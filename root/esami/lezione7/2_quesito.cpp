// Scrivi Macro in cui

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"

// 1) Si definiscono 2 istogrammi monodimensionali di 100 bin tra 0 e 10
// 2) Si riempie il primo istogramma con 10^5 occorrenze di una variabile casuale x generata esplicitamente e singolarmente e distribuita
// secondo una gaussiana con media 5 e dev std 1
// 3) Su tali occorrenze, si simula (attraverso "hit or miss") un'efficienza di rivelazione dipendente dalla variabile casuale x secondo
// E(x) = x / 10. Riempire il secondo istogramma con le occorrenze accettate
// 4) Si effettui la divisione fra i due istogrammi, usando Divide, e usando l'opportuna opzione per la valutazione degli errori secondo
// la statistica binomiale
// 5) Stampa i tre istorgammi in una unica canvas divisa in tre pad

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
}

void myMacro()
{
  gRandom->SetSeed();
  constexpr Int_t nGen = 1E5;

  TH1F *h1 = new TH1F("h1", "histo1", 100, 0, 10);
  TH1F *h2 = new TH1F("h2", "histo2", 100, 0, 10);

  TF1 *func = new TF1("func", "x / 10.");

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Gaus(5., 1.); // Mean 5 e STD 1
    xRND = gRandom->Rndm();

    h1->Fill(x);

    if (xRND < func->Eval(x))
    {
      h2->Fill(x);
    }
  }

  // Efficienza = Accettate (h2) / Generate (h1)
  TH1F *hEff = new TH1F(*h1);
  hEff->SetTitle("hEff");
  hEff->SetName("hEff");
  hEff->Divide(h2, h1, 1, 1, "B");

  TCanvas *canvas = new TCanvas("canvas");
  canvas->Divide(2, 2);
  canvas->cd(1);
  h1->Draw();
  canvas->cd(2);
  h2->Draw();
  canvas->cd(3);
  hEff->Draw();
}