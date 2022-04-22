#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include <fstream>
#include <iostream>

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(1111);
}

void computeRMS()
{
  constexpr int N = 1000;
  // ***** CALCOLO DEVIAZIONE STANDARD DEL RUMORE *****
  TH1F *stdHisto = new TH1F("stdHisto", "Rumore", N, -5.0, -4.8);
  TH1F *ampiezza4kHisto = new TH1F("ampiezza4kHisto", "Rumore", N, 4.9, 5.1);
  TH1F *ampiezza10kHisto = new TH1F("ampiezza10kHisto", "Rumore", N, 4.9, 5.1);
  std::ifstream in;
  in.open("rumore.txt");
  Float_t rumore;
  while (1)
  {
    in >> rumore;
    if (!in.good())
    {
      break;
    }
    stdHisto->Fill(rumore);
  }
  in.close();

  in.open("ampiezza4k.txt");
  Float_t ampiezza4k;
  while (1)
  {
    in >> ampiezza4k;
    if (!in.good())
    {
      break;
    }
    ampiezza4kHisto->Fill(ampiezza4k);
  }
  in.close();

  in.open("ampiezza10k.txt");
  Float_t ampiezza10k;
  while (1)
  {
    in >> ampiezza10k;
    if (!in.good())
    {
      break;
    }
    ampiezza10kHisto->Fill(ampiezza10k);
  }
  in.close();

  Double_t stdDev = stdHisto->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE *****" << '\n'
            << "Deviazione Standard: " << stdDev << '\n'
            << "Underflows: " << stdHisto->GetBinContent(0) << '\n'
            << "Overflows: " << stdHisto->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezza4kDev = ampiezza4kHisto->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD AMPIEZZA 4k *****" << '\n'
            << "Deviazione Standard: " << ampiezza4kDev << '\n'
            << "Underflows: " << ampiezza4kHisto->GetBinContent(0) << '\n'
            << "Overflows: " << ampiezza4kHisto->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezza10kDev = ampiezza10kHisto->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD AMPIEZZA 10k *****" << '\n'
            << "Deviazione Standard: " << ampiezza10kDev << '\n'
            << "Underflows: " << ampiezza10kHisto->GetBinContent(0) << '\n'
            << "Overflows: " << ampiezza10kHisto->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  TCanvas *cRMS = new TCanvas("cRMS");
  cRMS->Divide(1, 2);
  cRMS->cd(1);
  stdHisto->Draw();
  cRMS->cd(2);
  ampiezza4kHisto->Draw();

  TCanvas *c2 = new TCanvas("c2");
  c2->cd();
  ampiezza10kHisto->Draw();
}