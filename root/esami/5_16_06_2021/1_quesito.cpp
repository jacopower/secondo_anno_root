// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// 1. Si definiscono 2 istogrammi monodimensionali di 500 bin in un range da 0 a 5.
// 1. Si riempie il primo istogramma con 107 occorrenze di una variabile casuale x generate esplicitamente e singolarmente (i.e. attraverso gRandom) e distribuite secondo una distribuzione esponenziale decrescente con media μ=1 (campione totale).
// 2. Su tali occorrenze, si simula (attraverso un criterio di reiezione di tipo “hit or miss”) un’efficienza di rivelazione dipendente dalla variabile casuale x secondo la forma funzionale:
//• ε(x)=x2/2 per 0<=x<=1
//• ε(x)=0.5 per 1<x<5
//• ε(x)=0. per ogni altro valore di x
// Riempire il secondo istogramma con le occorrenze accettate (campione accettato).
// 3. Si effettua la divisione fra i due istogrammi per ottenere l’efficienza di rivelazione osservata, utilizzando il metodo Divide della classe degli istogrammi e inserendo l’opportuna opzione per la valutazione degli errori secondo la statistica binomiale.
// 4. Si disegna l’istogramma dell’ efficienza visualizzando anche le incertezze sui contenuti dei bin.
// 5. Si valuta, usando il metodo che restituisce il numero totale di ingressi, l’efficienza integrale (occorrenze accettate/occorrenze generate)

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <iostream>

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  TH1F *hGen = new TH1F("hGen", "Generated", 500, 0., 5.);
  TH1F *hAcc = new TH1F("hAcc", "Accepted", 500, 0., 5.);
  TF1 *f = new TF1("f", "(x >= 0 && x <= 1) * x * x / 2. + (x > 1 && x < 5) * 0.5 + (x >= 5) * 0.");

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < 1E7; ++i)
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
  hEff->SetTitle("Observed Efficiency");
  hEff->SetName("hEff");
  hEff->Divide(hAcc, hGen, 1, 1, "B");

  TCanvas *canvas = new TCanvas("canvas");
  hEff->Draw("H");
  hEff->Draw("E, SAME");

  Double_t genEntries = hGen->GetEntries();
  Double_t accEntries = hAcc->GetEntries();
  std::cout << "***** INTEGRAL EFFICIENCY *****" << '\n';
  std::cout << "Integral Efficiency: " << accEntries / genEntries << '\n';
  std::cout << "Prova: " << hAcc->Integral() / hGen->Integral() << '\n';
  std::cout << "**********" << '\n';
}