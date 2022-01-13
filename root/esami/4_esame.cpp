// Scrivi Macro in cui si monitora, con la classe TBenchmark, il tempo di CPU rispettivamente impiegato per fare le due seguenti operazioni,
// utilizzando anche l'opportuno metodo per stampare a schermo i tempi di esecuzione
// OP1: Generare 10^5 occorrenze generate espliicitamente e singolarmente di una gaussiana con media 0 e std dev 1, riempiendo un istogramma h1
// OP2: Generare 10^5 occorrenze di una gaussiana G(0,1) con il metodo FillRandom, riempiendo un istogramma h2

#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TBenchmark.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  constexpr Int_t nGen = 1E5;

  TH1F *h1 = new TH1F("h1", "histo1", 100, -5, 5);
  TH1F *h2 = new TH1F("h2", "histo2", 100, -5, 5);

  // Riempio h1 (OP1)
  Double_t xRND = 0;

  gBenchmark->Start("With TH1::Gaus invocation");
  for (Int_t i = 0; i < nGen; ++i)
  {
    xRND = gRandom->Gaus(0, 1);
    h1->Fill(xRND);
  }
  gBenchmark->Show("With TH1::Gaus invocation");

  // Riempio h2 (OP2)
  gBenchmark->Start("With TH1::FillRandom");
  h2->FillRandom("gaus", nGen);
  gBenchmark->Show("With TH1::FillRandom");

  TCanvas *canvas = new TCanvas("canvas");
  canvas->Divide(1, 2);
  canvas->cd(1);
  h1->Draw();
  canvas->cd(2);
  h2->Draw();
}