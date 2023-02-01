// Scrivi una Macro in cui
// 1) Si definisce un istogramma monodimensionale di 100 bin in un range da 0 a 10.
// 2) Si riempe l’istogramma con 10^7 occorrenze di una variabile casuale x distribuita secondo la p.d.f.
// f(x) = sqrt(x) + x^2 nel range [0,10] , utilizzando il metodo FillRandom(const char* f,Int_t N) della classe di istogrammi.

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
}

void myMacro()
{
  constexpr Int_t nGen = 1E7;

  TH1F *h = new TH1F("h", "histo", 100, 0., 10.);       // creating histo of 100 bin, in [0,10]
  TF1 *f = new TF1("f", "sqrt(x) + x ^2", 0., 10.);     // crating a functional form (function=sqrt(x) + x ^2) in [0.10]

  h->FillRandom("f", nGen);                             // fill h 1E7 times with random data from the p.d.f. defined in f

  TCanvas *canvas = new TCanvas("canvas");              // creates a canvas
  h->Draw();                                            // draws h (not needed)
}
