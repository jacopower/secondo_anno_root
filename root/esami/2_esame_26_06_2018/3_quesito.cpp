// Scrivi una Macro in cui
// 1) Si definisce un istogramma monodimensionale di 100 bin in un range da 0 a 10.
// 2) Si riempe l’istogramma con 10^5 occorrenze di una variabile casuale x distribuita secondo la p.d.f.
// f(x) = sin(x) + x^2 nel range [0,10] , utilizzando il metodo FillRandom(const char* f,Int_t N).

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
  constexpr Int_t nGen = 1E5;

  TH1F *h = new TH1F("h", "histo", 100, 0., 10.);
  TF1 *f = new TF1("f", "sin(x) + x^2", 0., 10.);

  h->FillRandom("f", nGen);

  TCanvas *canvas = new TCanvas("canvas");
  h->Draw();
}