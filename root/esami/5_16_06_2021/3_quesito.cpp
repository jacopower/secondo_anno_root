// Si scriva la parte rilevante e autoconsistente di codice di una macro di ROOT rivolto a monitorare con la classe TBenchmark il tempo di CPU rispettivamente impiegato per fare le due seguenti operazioni, utilizzando anche l’opportuno metodo per stampare a schermo i tempi di esecuzione:
// 1. Operazione 1: Generare 105 occorrenze generate esplicitamente e singolarmente di una gaussiana con media μ=0 e deviazione standard σ=1, riempendo un istogramma h1 che si assume già opportunamente definito.
// 2. Operazione 2: Generare 105 occorrenze di una Gaussiana G(0,1) con il metodo FillRandom, riempendo un istogramma h2 che si assume già opportunamente definito.

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TBenchmark.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  constexpr Int_t nGen = 1E5;
  TH1F *h1 = new TH1F("h1", "test1", 1000, -5., 5.);
  TH1F *h2 = new TH1F("h2", "test2", 1000, -5., 5.);

  Double_t x = 0.;
  gBenchmark->Start("With::Gaus Invocation");
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Gaus(0., 1.);
    h1->Fill(x);
  }
  gBenchmark->Show("With::Gaus Invocation");

  gBenchmark->Start("With::Fill Random");
  h2->FillRandom("gaus", nGen);
  gBenchmark->Show("With::Fill Random");
}