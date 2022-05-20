// Scrivi macro in cui
// 1) Si definisce un istogramma monodimensionale di 4 bin in un range da 0 a 4
// 2) Genera una popolazione di 10^6 elementi appartenenti a 4 categorie diverse, e riempi l'istogramma, con probabilità:
// caso 0: 60%
// caso 1: 30%
// caso 2: 9%
// caso 3: 1%
// 3) Stampa a schermo il contenuto dei quattro bin e la relativa incertezza, usando i metodi espliciti degli istogrammi GetBinContent() e GetBinError()

#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
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
  constexpr Int_t nGen = 1E6;

  TH1F *h = new TH1F("h", "histo", 4, 0., 4.);

  Double_t xRND = 0.;

  for (Int_t i = 0; i < nGen; ++i)
  {
    xRND = gRandom->Rndm();

    if (xRND < 0.6) // 60 %
    {
      h->Fill(0);
    }
    else if (xRND < 0.9) // 30%
    {
      h->Fill(1);
    }
    else if (xRND < 0.99) // 9%
    {
      h->Fill(2);
    }
    else // 1%
    {
      h->Fill(3);
    }
  }

  TCanvas *canvas = new TCanvas("canvas");
  h->Draw();

  std::cout << "***** CONTENUTI DEI BIN *****" << '\n';
  std::cout << "Unferflows: " << h->GetBinContent(0) << " +/- " << h->GetBinError(0) << '\n';
  std::cout << "Bin 1: " << h->GetBinContent(1) << " +/- " << h->GetBinError(1) << '\n';
  std::cout << "Bin 2: " << h->GetBinContent(2) << " +/- " << h->GetBinError(2) << '\n';
  std::cout << "Bin 3: " << h->GetBinContent(3) << " +/- " << h->GetBinError(3) << '\n';
  std::cout << "Bin 4: " << h->GetBinContent(4) << " +/- " << h->GetBinError(4) << '\n';
  std::cout << "Overflows: " << h->GetBinContent(5) << " +/- " << h->GetBinError(5) << '\n';
}