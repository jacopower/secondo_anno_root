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
#include "TMultiGraph.h"

// R = 149.83 Ohm
// C = 158.4 nF
// L = 10.43 mH

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(1111);
}

Double_t amp_time_resistenza(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[4] * par[1] * par[3]))) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
}

Double_t amp_time_induttanza(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[4] * par[2] * par[0] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[1] * par[4] * par[3])) + TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
}

Double_t amp_time_condensatore(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = W
  Double_t xx = x[0];
  Double_t val = par[0] / par[4] / par[3] * TMath::Cos(par[4] * xx + TMath::ATan((1 - par[4] * par[4] * par[2] * par[3]) / (par[1] * par[4] * par[3])) - TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * par[2] - 1 / (par[4] * par[3])) * (par[4] * par[2] - 1 / (par[4] * par[3])));
  return val;
}

Double_t amp_freq_resistenza(Double_t *x, Double_t *par)
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] / TMath::Sqrt(par[1] * par[1] + (xx * par[2] - 1 / (xx * par[3])) * (xx * par[2] - 1 / (xx * par[3])));
  return val;
}

Double_t amp_freq_induttanza(Double_t *x, Double_t *par)
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = xx * par[2] * par[0] / TMath::Sqrt(par[1] * par[1] + (xx * par[2] - 1 / (xx * par[3])) * (xx * par[2] - 1 / (xx * par[3])));
  return val;
}

Double_t amp_freq_condensatore(Double_t *x, Double_t *par)
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = par[0] / xx / par[3] / TMath::Sqrt(par[1] * par[1] + (xx * par[2] - 1 / (xx * par[3])) * (xx * par[2] - 1 / (xx * par[3])));
  return val;
}

Double_t phase_freq_resistenza(Double_t *x, Double_t *par)
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * xx * par[1] * par[2]) / (xx * par[0] * par[2]));
  return val;
}

Double_t phase_freq_induttanza(Double_t *x, Double_t *par)
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * xx * par[1] * par[2]) / (xx * par[0] * par[2])) + TMath::PiOver2();
  return val;
}

Double_t phase_freq_condensatore(Double_t *x, Double_t *par)
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * xx * par[1] * par[2]) / (xx * par[0] * par[2])) - TMath::PiOver2();
  return val;
}

void rumore() // CALCOLO DEVIAZIONE STANDARD DAL RUMORE
{
  constexpr int N = 1000; // VEDI QUESTO
  TH1F *histoOndaQuadra = new TH1F("histoOndaQuadra", "Rumore con Onda Quadra", N, 4.9, 5.1);
  TH1F *histo1k = new TH1F("histo1k", "Rumore a 1kHz", N, 4.9, 5.1);
  TH1F *histo4k = new TH1F("histo4k", "Rumore a 4kHz", N, 4.9, 5.1);
  TH1F *histo10k = new TH1F("histo10k", "Rumore a 10kHz", N, 4.9, 5.1);
  TH1F *histoFase = new TH1F("histoFase", "Rumore Fase", N, 4.9, 5.1); // ATTENZIONE A NUMERO BIN E ESTREMI

  std::ifstream in; // VA BENE SE APRO E CHIUDO QUESTO?

  // LEGGO RUMORE A 1K
  in.open("data/rumore/ampiezza1k.txt");
  Float_t ampiezza1k;
  while (1)
  {
    in >> ampiezza1k;
    if (!in.good())
    {
      break;
    }
    histo1k->Fill(ampiezza1k);
  }
  in.close();

  // LEGGO RUMORE A 4K
  in.open("data/rumore/ampiezza4k.txt");
  Float_t ampiezza4k;
  while (1)
  {
    in >> ampiezza4k;
    if (!in.good())
    {
      break;
    }
    histo4k->Fill(ampiezza4k);
  }
  in.close();

  // LEGGO RUMORE A 10K
  in.open("data/rumore/ampiezza10k.txt");
  Float_t ampiezza10k;
  while (1)
  {
    in >> ampiezza10k;
    if (!in.good())
    {
      break;
    }
    histo10k->Fill(ampiezza10k);
  }
  in.close();

  // LEGGO RUMORE ONDA QUADRA
  in.open("data/rumore/rumore_onda_quadra.txt");
  Float_t ampiezzaOndaQuadra;
  while (1)
  {
    in >> ampiezzaOndaQuadra;
    if (!in.good())
    {
      break;
    }
    histoOndaQuadra->Fill(ampiezzaOndaQuadra);
  }
  in.close();

  // LEGGO RUMORE FASE
  in.open("data/rumore/rumore_fase.txt");
  Float_t ampiezzaFase;
  while (1)
  {
    in >> ampiezzaFase;
    if (!in.good())
    {
      break;
    }
    histoFase->Fill(ampiezzaFase);
  }
  in.close();

  // STAMPO LE DEVIAZIONI STANDARD
  Double_t ampiezza1kDev = histo1k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE 1kHz *****" << '\n'
            << "Deviazione Standard: " << ampiezza1kDev << '\n'
            << "Underflows: " << histo1k->GetBinContent(0) << '\n'
            << "Overflows: " << histo1k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezza4kDev = histo4k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE 4kHz *****" << '\n'
            << "Deviazione Standard: " << ampiezza4kDev << '\n'
            << "Underflows: " << histo4k->GetBinContent(0) << '\n'
            << "Overflows: " << histo4k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezza10kDev = histo10k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE 10kHz *****" << '\n'
            << "Deviazione Standard: " << ampiezza10kDev << '\n'
            << "Underflows: " << histo10k->GetBinContent(0) << '\n'
            << "Overflows: " << histo10k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaOndaQuadraDev = histoOndaQuadra->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE ONDA QUADRA (CIRCUITO VUOTO) *****" << '\n'
            << "Deviazione Standard: " << ampiezzaOndaQuadraDev << '\n'
            << "Underflows: " << histoOndaQuadra->GetBinContent(0) << '\n'
            << "Overflows: " << histoOndaQuadra->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev = histoFase->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev << '\n'
            << "Underflows: " << histoFase->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  // PLOT ISTOGRAMMI
  TCanvas *c = new TCanvas;
  c->Divide(2, 2);
  c->cd(1);
  histo1k->Draw();
  c->cd(2);
  histo4k->Draw();
  c->cd(3);
  histo10k->Draw();
  c->cd(4);
  histoOndaQuadra->Draw();

  TCanvas *c2 = new TCanvas;
  c2->cd();
  histoFase->Draw();
}

// CHE ERRORE SU Y ASSOCIARE QUI?
void amplitude_sweep()
{
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; x(UDM); y(UDM)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; x(UDM); y(UDM)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_ampiezza/sweep_freq_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; x(UDM); y(UDM)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg");
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

  TCanvas *multiCanvas = new TCanvas();
  multiCanvas->cd();
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA DA LP?
  multiCanvas->BuildLegend();; // Vedi cosa fa
}

// CHE ERRORE SU Y ASSOCIARE QUI?
void phase_sweep()
{
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_fase/sweep_phase_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; x(UDM); y(UDM)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_fase/sweep_phase_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; x(UDM); y(UDM)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_fase/sweep_phase_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; x(UDM); y(UDM)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_fase/sweep_phase_totale.txt", "%lg %lg %lg");
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

  TCanvas *multiCanvas = new TCanvas();
  multiCanvas->cd();
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA DA LP?
  multiCanvas->BuildLegend();; // Vedi cosa fa
}