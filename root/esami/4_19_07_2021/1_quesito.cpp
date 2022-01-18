// Si scriva la parte rilevante ed autoconsistente del codice di una macro di ROOT in cui:
// Si definiscono 2 istogrammi unidimensionali di 500 bin in un range da 0 a 5.
// Si riempie il primo istogramma con 107 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) e distribuite uniformemente nel range [0,5] (campione totale).
// Su tali occorrenze, si simula (attraverso un criterio di reiezione di tipo “hit or miss”) un’efficienza di rivelazione dipendente dalla variabile casuale x secondo la forma:
//ε(x)= e-x per 0<=x<=3
//ε(x)=0.05 per x>3
// Riempire il secondo istogramma con le occorrenze accettate (campione accettato).
// Si effettua la divisione fra i due istogrammi per ottenere l’efficienza di rivelazione osservata, utilizzando il metodo Divide della classe degli istogrammi e inserendo l’opportuna opzione per la valutazione degli errori secondo la statistica binomiale.
// Si disegna l’istogramma dell’efficienza visualizzando le incertezze sui contenuti dei bin.
// Si stampa a schermo il valore dell’efficienza integrale.

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
  TH1F *hGen = new TH1F("hGen", "Generated", 500, 0., 5.);
  TH1F *hAcc = new TH1F("hAcc", "Accepted", 500, 0., 5.);
  TF1 *f = new TF1("f", "(x >= 0 && x <= 3) * e^-x + (x > 3) * 0.05");

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Uniform(0., 5.);
    hGen->Fill(x);

    xRND = gRandom->Rndm();
    if (xRND < f->Eval(x))
    {
      hAcc->Fill(x);
    }
  }

  TH1F *hEff = new TH1F(*hGen);
  hEff->SetTitle("Observed Efficiency");
  hEff->SetName("hEff");
  hEff->Divide(hAcc, hGen, 1, 1, "B");

  TCanvas *canvas = new TCanvas("canvas");
  hEff->Draw("H");
  hEff->Draw("E, SAME");

  std::cout << "***** INTEGRAL EFFICIENCY *****" << '\n';
  std::cout << "Integral efficiency: " << hAcc->Integral() / hGen->Integral() << '\n';
  std::cout << "**********" << '\n';
}