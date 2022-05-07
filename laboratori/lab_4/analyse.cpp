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

// R = 149.83 Ohm
// C = 158.4 nF
// L = 10.43 mH


void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(1111);
}

Double_t freq_time_resistenza(Double_t *x, Double_t *par)
{
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[4] * par[1] * par[3]))) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
}

Double_t freq_time_induttanza(Double_t *x, Double_t *par)
{
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[4] * par[2] * par[0] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[1])) + TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
}

Double_t freq_time_condensatore(Double_t *x, Double_t *par)
{
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[0] / par[4] / par[3] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[1])) - TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
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
            << " ***** DEVIAZIONE STANDARD RUMORE 1k *****" << '\n'
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

void amplitude_sweep()
{
  TGraphErrors *graphResistenza = new TGraphErrors("sweep_freq_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; x(UDM); y(UDM)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("sweep_freq_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; x(UDM); y(UDM)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("sweep_freq_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; x(UDM); y(UDM)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("sweep_freq_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; x(UDM); y(UDM)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  TCanvas *c1 = new TCanvas();
  c1->Divide(2, 2);
  c1->cd(1);
  graphResistenza->Draw("APE");
  c1->cd(2);
  graphInduttanza->Draw("APE");
  c1->cd(3);
  graphCondensatore->Draw("APE");
  c1->cd(4);
  graphTotale->Draw("APE");
}

void phase_sweep()
{
  TGraphErrors *graphResistenza = new TGraphErrors("sweep_phase_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; x(UDM); y(UDM)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("sweep_phase_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; x(UDM); y(UDM)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("sweep_phase_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; x(UDM); y(UDM)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("sweep_phase_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; x(UDM); y(UDM)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  TCanvas *c1 = new TCanvas();
  c1->Divide(2, 2);
  c1->cd(1);
  graphResistenza->Draw("APE");
  c1->cd(2);
  graphInduttanza->Draw("APE");
  c1->cd(3);
  graphCondensatore->Draw("APE");
  c1->cd(4);
  graphTotale->Draw("APE");
}

void tensioni_tempo()
{
  TGraphErrors *graph1 = new TGraphErrors("tensione1.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; x(UDM); y(UDM)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);
}