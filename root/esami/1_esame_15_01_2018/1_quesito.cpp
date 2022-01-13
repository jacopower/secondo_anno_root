// Scrivi Macro in cui
// 1) si definiscono 2 istogrammi monodimensionali di 1000 bin tra 0 e 5
// 2) Si riempie il primo con 10^7 occorrenze di una variabile x generata esplicitamente e singolarmente (gRandom) e distribuita
// secondo una esponenziale decrescente con media 1
// 3) Su tali occorrenze si simula (hit or miss) un'efficienza dipendente da x secondo E(x) = x / 5. Riempire il secondo istogramma
// con le occorrenze accettate
// 4) Si effettua la divisione fra i due istogrammi per ottenere l’efficienza di rivelazione osservata,
// utilizzando il metodo Divide e inserendo l’opportuna opzione per la valutazione degli errori secondo la statistica binomiale.
// 5) Si disegna l’istogramma dell’ efficienza visualizzando le incertezze sui contenuti dei bin.

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  constexpr Int_t nGen = 1E7;

  TH1F *hGen = new TH1F("hGen", "generazione", 1000, 0., 5.);
  TH1F *hAcc = new TH1F("hAcc", "accettati", 1000, 0., 5.);

  TF1 *f = new TF1("f", "x/5.", 0., 5.);

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Exp(1);
    hGen->Fill(x);

    xRND = gRandom->Rndm();
    if (xRND < f->Eval(x))
    {
      hAcc->Fill(x);
    }
  }

  TH1F *hEff = new TH1F(*hGen);
  hEff->SetTitle("0bserved Efficiency");
  hEff->SetName("hEff");
  hEff->Divide(hAcc, hGen, 1, 1, "B"); // B = statistica binomiale

  TCanvas *canvas = new TCanvas("canvas");
  hEff->Draw("H");
  hEff->Draw("E, SAME"); // E = incertezze sui bin

  TCanvas *c = new TCanvas("c");
  c->Divide(1, 2);
  c->cd(1);
  hGen->Draw();
  c->cd(2);
  hAcc->Draw();
}