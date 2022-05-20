// Si scriva la parte rilevante e autoconsistente del codice di una macro di ROOT in cui:
// 1. Si definisce un istogramma monodimensionale di 4 bin in un range da 0 a 4, finalizzato a visualizzare il numero di conteggi osservati nelle categorie che contraddistinguono la popolazione di cui a punto successivo .
// 2. Si genera una popolazione di N_TOT=10^5 elementi appartenenti a 4 categorie distinte, con rispettive probabilità di occorrenza:
// o caso “0”: 1/8
// o caso “1”: 1/4
// o caso “2”: 1/2
// o caso “3”: 1/8
// e si riempie l’istogramma con le occorrenze osservate nei diversi casi
// 3. si stampano a schermo le frazioni di popolazione osservate nei quattro bin dell’istogramma, con relativa incertezza, usando i metodi espliciti degli istogrammi per estrarre il numero di ingressi totali dell’istogramma, il contenuto dei bin e il relativo errore.

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <iostream>

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  constexpr Int_t Ntot = 1E5;
  TH1F *hObs = new TH1F("hObs", "Osservati", 4, 0., 4.);

  Double_t xRND = 0;
  for (Int_t i = 0; i < Ntot; ++i)
  {
    xRND = gRandom->Rndm();
    if (xRND < 0.125) //(1/8)%
    {
      hObs->Fill(0);
    }
    else if (xRND < 0.375) //(1/4)%
    {
      hObs->Fill(1);
    }
    else if (xRND < 0.875) //(1/2)%
    {
      hObs->Fill(2);
    }
    else //(1/8)%
    {
      hObs->Fill(3);
    }
  }

  TCanvas *canvas = new TCanvas("canvas");
  hObs->Draw();

  Double_t entries = hObs->GetEntries();
  Double_t underflows = hObs->GetBinContent(0);
  Double_t case0 = hObs->GetBinContent(1);
  Double_t case0Err = hObs->GetBinError(1);
  Double_t case1 = hObs->GetBinContent(2);
  Double_t case1Err = hObs->GetBinError(2);
  Double_t case2 = hObs->GetBinContent(3);
  Double_t case2Err = hObs->GetBinError(3);
  Double_t case3 = hObs->GetBinContent(4);
  Double_t case3Err = hObs->GetBinError(4);
  Double_t overflows = hObs->GetBinContent(5);

  std::cout << "***** BIN CONTENTS *****" << '\n';
  std::cout << "Entries: " << entries << '\n';
  std::cout << "Underflows: " << underflows << '\n';
  std::cout << "Overflows: " << overflows << '\n';
  std::cout << "Case 0: " << case0 / entries << " +/- " << case0Err / entries << '\n';
  std::cout << "Case 1: " << case1 / entries << " +/- " << case1Err / entries << '\n';
  std::cout << "Case 2: " << case2 / entries << " +/- " << case2Err / entries << '\n';
  std::cout << "Case 3: " << case3 / entries << " +/- " << case3Err / entries << '\n';
  std::cout << "**********" << '\n';
}