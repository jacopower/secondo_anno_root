// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// 1. Si definisce un istogramma monodimensionale di 1000 bin in un range da 0 a 5.
// 2. Si riempie l’istogramma con 10^8 occorrenze di una variabile casuale x distribuita secondo la p.d.f. (definire esplicitamente una unica funzione che descriva la p.d.f. in esame):
// f(x)=x² /8 per 0<=x <2
// f(x)=0.5 per 2<=x<=5
// nel range [0,5] , utilizzando il metodo FillRandom(const char* f,Int_t N) della classe di istogrammi.

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
  constexpr Int_t nGen = 1E8;
  TH1F *h = new TH1F("h", "histo", 1000, 0., 5.);
  TF1 *f = new TF1("f", "(x >= 0 && x < 2) * x * x / 8. + (x >= 2 && x <= 5) * 0.5", 0., 5.);
  h->FillRandom("f", nGen);

  TCanvas *canvas = new TCanvas("canvas");
  h->Draw();
}