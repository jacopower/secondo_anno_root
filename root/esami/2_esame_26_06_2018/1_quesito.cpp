// Scrivi una Macro in cui
// 1) Si definiscono 2 istogrammi monodimensionali di 1000 bin in un range da 0 a 10.
// 2) Si riempe il primo istogramma con 10^7 occorrenze di una variabile casuale x generate esplicitamente e singolarmente e distribuite
// secondo una gaussiana con media 5 e deviazione standard 1 (campione totale).
// 3) Su tali occorrenze, si simula (“hit or miss”) un’efficienza di rivelazione dipendente dalla variabile casuale x secondo la forma:
// E(x) = 0.1 * x * e ^ -x. Riempire il secondo istogramma con le occorrenze accettate (campione accettato).
// 4) Si effettua la divisione fra i due istogrammi per ottenere l’efficienza di rivelazione osservata, utilizzando Divide
// e inserendo l’opportuna opzione per la valutazione degli errori secondo la statistica binomiale.
// 5) Si disegna l’istogramma dell’ efficienza visualizzando le incertezze sui contenuti dei bin.

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  gRandom->SetSeed();
  constexpr Int_t nGen = 1E7;

  TH1F *hGen = new TH1F("hGen", "generated", 1000, 0., 10.);
  TH1F *hAcc = new TH1F("hAcc", "accepted", 1000, 0., 10.);
  TF1 *f = new TF1("f", "0.1 * x * e ^ -x", 0., 10.);

  Double_t x, xRND = 0.;
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Gaus(5., 1.);
    hGen->Fill(x);

    xRND = gRandom->Rndm();
    if (xRND < f->Eval(x))
    {
      hAcc->Fill(x);
    }
  }

  // Efficiency = Accepted / Generated
  TH1F *hEff = new TH1F(*hGen);
  hEff->SetTitle("Observed Efficiency");
  hEff->SetName("hEff");
  hEff->Divide(hAcc, hGen, 1, 1, "B");

  TCanvas *canvas = new TCanvas("canvas");
  hEff->Draw("H");
  hEff->Draw("E, SAME");

  TCanvas *c = new TCanvas("c");
  c->Divide(1, 2);
  c->cd(1);
  hGen->Draw();
  c->cd(2);
  hAcc->Draw();
}