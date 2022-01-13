// 1) Si definisce un istogramma monodimensionale di 100 bin in un range da 0 a 10
// 2) Si riempie l'isotgramma con 10^5 occorrenze di una variabile casuale x distribuita secondo la pdf f(x) = x nel range [0,10], utilizzando il metodo
// FillRandom(const char* f, Int_t N)

#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2211);
  gStyle->SetOptFit(1111);
}

void myMacro()
{
  constexpr Int_t nGen = 1E5;
  TH1F *histo = new TH1F("histo", "istogramma", 100, 0., 10.);

  TF1 *func = new TF1("func", "x", 0, 10);

  histo->FillRandom("func", nGen);

  TCanvas *canvas = new TCanvas();
  // TF1 *fitFunc = new TF1("fitFunc", "[0]*x", 0, 10);
  // fitFunc->SetParameter(0, 1);
  // histo->Fit("fitFunc");
  histo->Draw();
}

// VARIANTE: pdf definita come
// f(x) = x/5. per 0<=x<5
// f(x) = 1.   per 5<=x<=10

void myMacro2()
{
  constexpr Int_t nGen = 1E5;
  TH1F *histo = new TH1F("histo", "istogramma", 100, 0., 10.);

  TF1 *func = new TF1("func", "(x < 5 && x >= 0) * x / 5. + (x >= 5 && x <= 10)", 0, 10);
  histo->FillRandom("func", nGen);

  TCanvas *canvas = new TCanvas();
  // TF1 *fitfunc = new TF1("fitfunc", "[0] * (x < 5 && x >= 0) * x / 5. + (x >= 5 && x <= 10)", 0, 10);
  // fitfunc->SetParameter(0, 1);
  // histo->Fit("fitfunc");
  histo->Draw();
}

// VARIANTE
// Su tali occorrenze, si simula ("hit or miss") un'efficienza di rivelazione dipendente da x secondo
// E(x) = 0.3 per 0 <= x <= 3
// E(x) = 0.7 per x > 3

void myMacro3()
{
  gRandom->SetSeed();

  constexpr Int_t nGen = 1E5;
  TH1F *histo = new TH1F("histo", "istogramma", 100, 0., 10.);

  TF1 *func = new TF1("func", "x", 0, 10);
  TF1 *effFunc = new TF1("effFunc", "(x <= 3 && x >= 0) * 0.3 + (x > 3) * 0.7");

  Double_t x, xRND = 0;
  for (Int_t i = 0; i < nGen; ++i)
  {
    x = func->GetRandom();
    xRND = gRandom->Rndm();

    if (xRND < effFunc->Eval(x))
    {
      histo->Fill(x);
    }
  }

  // TH1F *h = new TH1F("histo", "istogramma", 100, 0., 10.);
  // h->FillRandom("func", nGen);

  TCanvas *canvas = new TCanvas();
  histo->Draw();
  // h->Draw();
}
