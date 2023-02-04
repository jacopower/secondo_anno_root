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
  gROOT->SetStyle("Plain");   // the Plain style gives"conventional" Postscript output (the only other style is Default)
  gStyle->SetPalette(57);     // sets the palette called kBird (kCool is num 109)
  gStyle->SetOptStat(2210);   // prints on screen (graph legend): name; entries; mean; rms(stdev); overflows; underflows; integral; skewness; kurtosis
}

void myMacro()
{
  constexpr Int_t nGen = 1E7;

  TH1F *hGen = new TH1F("hGen", "generazione", 1000, 0., 5.);     // creating a histogram (pointer to an histo) of 1E3 bins in [0,5]
  TH1F *hAcc = new TH1F("hAcc", "accettati", 1000, 0., 5.);       // creating another one

  TF1 *f = new TF1("f", "x/5.", 0., 5.);    // creates a functional form (function = x/5) in [0,5]

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < nGen; ++i)        // generating x and filling hGen 1E7 times
  {
    x = gRandom->Exp(1);                  // x is randomly set following a exponential descending distribution
    hGen->Fill(x);                        // a bin of hGen is filled with x

    xRND = gRandom->Rndm();               /// xRND is random from uniform distribution in [0,1]
    if (xRND < f->Eval(x))                // if xRND < x/5 (=f(x))
    {
      hAcc->Fill(x);
    }
  }

  TH1F *hEff = new TH1F(*hGen);
  hEff->SetTitle("0bserved Efficiency");
  hEff->SetName("hEff");
  hEff->Divide(hAcc, hGen, 1, 1, "B");  // B = statistica binomiale //  hEff = h1/h2 (following B?)

  TCanvas *canvas = new TCanvas("canvas");        // creating new canvas (where the drawing goes)
  hEff->Draw("H");                                // draws the histo with a continuous line
  hEff->Draw("E, SAME");                          // draws hEff error (uncertanties) bars, on the same graph as before

  TCanvas *c = new TCanvas("c");                  // (wasnt requested...but ok)
  c->Divide(1, 2);                                // divide the canvas in a 1x2 grid    
  c->cd(1);                                       // selects the 1st part of the grid
  hGen->Draw();                                   // draw hGen on c, in the grid space 1 ???
  c->cd(2);
  hAcc->Draw();
}
