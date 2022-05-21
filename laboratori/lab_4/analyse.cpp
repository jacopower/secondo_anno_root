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

/*
// FORMULA f = W / (2*pi)
// R = 149.83 Ohm
// C = 158.4 nF
// L = 10.43 mH

Double_t amp_freq_resistenza(Double_t *x, Double_t *par)
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] / TMath::Sqrt(par[1] * par[1] + (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])) * (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])));
  return val;
}

Double_t amp_freq_induttanza(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = xx * TMath::Pi() * 2 * par[2] * par[0] / TMath::Sqrt(par[1] * par[1] + (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])) * (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])));
  return val;
}

Double_t amp_freq_condensatore(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = par[0] / xx * TMath::Pi() * 2 / par[3] / TMath::Sqrt(par[1] * par[1] + (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])) * (xx * TMath::Pi() * 2 * par[2] - 1 / (xx * TMath::Pi() * 2 * par[3])));
  return val;
}
*/

// VEDI QUESTI PARAMETRI
constexpr Double_t V0_mis = 3.21514;
constexpr Double_t R_mis = 150.47;
constexpr Double_t L_mis = 11.46 * 1E-3;
constexpr Double_t C_mis = 157.8 * 1E-9;
constexpr Double_t R_tot = 200;

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(111111);
}

Double_t amp_time_resistenza(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = f
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] * TMath::Cos(par[4] * TMath::Pi() * 2 * xx + TMath::ATan((1 - par[4] * TMath::Pi() * 2 * par[4] * TMath::Pi() / 2 * par[2] * par[3]) / (par[4] * TMath::Pi() * 2 * par[1] * par[3]))) / TMath::Sqrt(par[1] * par[1] + (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])) * (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])));
  return val;
}

Double_t amp_time_induttanza(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = f
  Double_t xx = x[0];
  Double_t val = par[4] * TMath::Pi() * 2 * par[2] * par[0] * TMath::Cos(par[4] * TMath::Pi() * 2 * xx + TMath::ATan((1 - par[4] * TMath::Pi() * 2 * par[4] * TMath::Pi() * 2 * par[2] * par[3]) / (par[1] * par[4] * TMath::Pi() * 2 * par[3])) + TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])) * (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])));
  return val;
}

Double_t amp_time_condensatore(Double_t *x, Double_t *par)
{
  // 5 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = f
  Double_t xx = x[0];
  Double_t val = par[0] / par[4] * TMath::Pi() * 2 / par[3] * TMath::Cos(par[4] * TMath::Pi() * 2 * xx + TMath::ATan((1 - par[4] * TMath::Pi() * 2 * par[4] * TMath::Pi() * 2 * par[2] * par[3]) / (par[1] * par[4] * TMath::Pi() * 2 * par[3])) - TMath::PiOver2()) / TMath::Sqrt(par[1] * par[1] + (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])) * (par[4] * TMath::Pi() * 2 * par[2] - 1 / (par[4] * TMath::Pi() * 2 * par[3])));
  return val;
}

Double_t amp_freq_resistenza(Double_t *x, Double_t *par)
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  // par[4] = Rtot
  Double_t xx = x[0];
  Double_t val = par[1] * par[0] / (TMath::Sqrt(par[4] * par[4] + (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3])) * (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3]))));
  return val;
}

Double_t amp_freq_induttanza(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = TMath::TwoPi() * xx * par[2] * par[0] / (TMath::Sqrt(par[1] * par[1] + (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3])) * (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3]))));
  return val;
}

Double_t amp_freq_condensatore(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 4 PARAMETRI
  // par[0] = V0
  // par[1] = R
  // par[2] = L
  // par[3] = C
  Double_t xx = x[0];
  Double_t val = par[0] / (TMath::TwoPi() * xx * par[3]) / (TMath::Sqrt(par[1] * par[1] + (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3])) * (TMath::TwoPi() * xx * par[2] - 1 / (TMath::TwoPi() * xx * par[3]))));
  return val;
}

Double_t phase_freq_resistenza(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * TMath::Pi() * 2 * xx * TMath::Pi() * 2 * par[1] * par[2]) / (xx * TMath::Pi() * 2 * par[0] * par[2]));
  return val;
}

Double_t phase_freq_induttanza(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * TMath::Pi() * 2 * xx * TMath::Pi() * 2 * par[1] * par[2]) / (xx * TMath::Pi() * 2 * par[0] * par[2])) + TMath::PiOver2();
  return val;
}

Double_t phase_freq_condensatore(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
{
  // 3 PARAMETRI
  // par[0] = R
  // par[1] = L
  // par[2] = C
  Double_t xx = x[0];
  Double_t val = TMath::ATan((1 - xx * TMath::Pi() * 2 * xx * TMath::Pi() * 2 * par[1] * par[2]) / (xx * TMath::Pi() * 2 * par[0] * par[2])) - TMath::PiOver2();
  return val;
}

void rumore()
{
  constexpr int N = 1000;
  TH1F *histoOndaQuadra = new TH1F("histoOndaQuadra", "Rumore con Onda Quadra", N, 4.9, 5.1);
  TH1F *histo1k = new TH1F("histo1k", "Rumore a 1kHz", N, 4.9, 5.1);
  TH1F *histo4k = new TH1F("histo4k", "Rumore a 4kHz", N, 4.9, 5.1);
  TH1F *histo10k = new TH1F("histo10k", "Rumore a 10kHz", N, 4.9, 5.1);
  TH1F *histo15k = new TH1F("histo15k", "Rumore a 15kHz", N, 4.9, 5.1);
  TH1F *histo20k = new TH1F("histo20k", "Rumore a 20kHz", N, 4.9, 5.1);
  TH1F *histoFase1k = new TH1F("histoFase1k", "Rumore Fase a 1k", N, 0.20, 0.21);
  TH1F *histoFase4k = new TH1F("histoFase4k", "Rumore Fase a 4k", N, 0.25, 0.26);
  TH1F *histoFase10k = new TH1F("histoFase10k", "Rumore Fase a 10k", N, 0.20, 0.22);
  TH1F *histoFase15k = new TH1F("histoFase15k", "Rumore Fase a 15k", N, 0.21, 0.22);
  TH1F *histoFase20k = new TH1F("histoFase20k", "Rumore Fase a 20k", N, 0.21, 0.23);

  std::ifstream in;

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

  // LEGGO RUMORE A 15K
  in.open("data/rumore/ampiezza15k.txt");
  Float_t ampiezza15k;
  while (1)
  {
    in >> ampiezza15k;
    if (!in.good())
    {
      break;
    }
    histo15k->Fill(ampiezza15k);
  }
  in.close();

  // LEGGO RUMORE A 20K
  in.open("data/rumore/ampiezza20k.txt");
  Float_t ampiezza20k;
  while (1)
  {
    in >> ampiezza20k;
    if (!in.good())
    {
      break;
    }
    histo20k->Fill(ampiezza20k);
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

  // LEGGO RUMORE FASE 1k
  in.open("data/rumore/fase1k.txt");
  Float_t ampiezzaFase1k;
  while (1)
  {
    in >> ampiezzaFase1k;
    if (!in.good())
    {
      break;
    }
    histoFase1k->Fill(ampiezzaFase1k);
  }
  in.close();

  // LEGGO RUMORE FASE 4k
  in.open("data/rumore/fase4k.txt");
  Float_t ampiezzaFase4k;
  while (1)
  {
    in >> ampiezzaFase4k;
    if (!in.good())
    {
      break;
    }
    histoFase4k->Fill(ampiezzaFase4k);
  }
  in.close();

  // LEGGO RUMORE FASE 10k
  in.open("data/rumore/fase10k.txt");
  Float_t ampiezzaFase10k;
  while (1)
  {
    in >> ampiezzaFase10k;
    if (!in.good())
    {
      break;
    }
    histoFase10k->Fill(ampiezzaFase10k);
  }
  in.close();

  // LEGGO RUMORE FASE 15k
  in.open("data/rumore/fase15k.txt");
  Float_t ampiezzaFase15k;
  while (1)
  {
    in >> ampiezzaFase15k;
    if (!in.good())
    {
      break;
    }
    histoFase15k->Fill(ampiezzaFase15k);
  }
  in.close();

  // LEGGO RUMORE FASE 20k
  in.open("data/rumore/fase20k.txt");
  Float_t ampiezzaFase20k;
  while (1)
  {
    in >> ampiezzaFase20k;
    if (!in.good())
    {
      break;
    }
    histoFase20k->Fill(ampiezzaFase20k);
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

  Double_t ampiezza15kDev = histo15k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE 15kHz *****" << '\n'
            << "Deviazione Standard: " << ampiezza15kDev << '\n'
            << "Underflows: " << histo15k->GetBinContent(0) << '\n'
            << "Overflows: " << histo15k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezza20kDev = histo20k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE 20kHz *****" << '\n'
            << "Deviazione Standard: " << ampiezza20kDev << '\n'
            << "Underflows: " << histo20k->GetBinContent(0) << '\n'
            << "Overflows: " << histo20k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaOndaQuadraDev = histoOndaQuadra->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE ONDA QUADRA (CIRCUITO VUOTO) *****" << '\n'
            << "Deviazione Standard: " << ampiezzaOndaQuadraDev << '\n'
            << "Underflows: " << histoOndaQuadra->GetBinContent(0) << '\n'
            << "Overflows: " << histoOndaQuadra->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev1k = histoFase1k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE A 1K *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev1k << '\n'
            << "Underflows: " << histoFase1k->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase1k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev4k = histoFase4k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE A 4K *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev4k << '\n'
            << "Underflows: " << histoFase4k->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase4k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev10k = histoFase10k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE A 10K *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev10k << '\n'
            << "Underflows: " << histoFase10k->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase10k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev15k = histoFase15k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE A 15K *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev15k << '\n'
            << "Underflows: " << histoFase15k->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase15k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  Double_t ampiezzaFaseDev20k = histoFase20k->GetStdDev();
  std::cout << '\n'
            << " ***** DEVIAZIONE STANDARD RUMORE FASE A 20K *****" << '\n'
            << "Deviazione Standard: " << ampiezzaFaseDev20k << '\n'
            << "Underflows: " << histoFase20k->GetBinContent(0) << '\n'
            << "Overflows: " << histoFase20k->GetBinContent(N + 1) << '\n'
            << "**********" << '\n';

  // FACCIO FIT LINEARE SUGLI ERRORI DELL'AMPIEZZA
  const Int_t nPoints = 5;
  Double_t yAmpiezza[nPoints] = {ampiezza1kDev, ampiezza4kDev, ampiezza10kDev, ampiezza15kDev, ampiezza20kDev};
  Double_t xAmpiezza[nPoints] = {1E3, 4E3, 10E3, 15E3, 20E3};
  TGraph *graphAmpiezza = new TGraph(nPoints, xAmpiezza, yAmpiezza);
  graphAmpiezza->SetTitle("Incertezze sulla misura dell'Ampiezza con Extract Single Tone");
  graphAmpiezza->GetXaxis()->SetTitle("Frequenza Generatore");
  graphAmpiezza->GetYaxis()->SetTitle("Incertezza su Ampiezza");
  graphAmpiezza->SetMarkerStyle(kFullCircle);
  graphAmpiezza->SetMarkerColor(kBlue);

  TF1 *f1 = new TF1("f1", "[0] + [1] * x", 1E3, 20E3);
  f1->SetParNames("q", "m");
  f1->SetLineColor(kRed);
  f1->SetLineWidth(2);
  graphAmpiezza->Fit(f1, "R");

  // FACCIO FIT LINEARE SUGLI ERRORI DELLA FASE
  Double_t yFase[nPoints] = {ampiezzaFaseDev1k, ampiezzaFaseDev4k, ampiezzaFaseDev10k, ampiezzaFaseDev15k, ampiezzaFaseDev20k};
  Double_t xFase[nPoints] = {1E3, 4E3, 10E3, 15E3, 20E3};
  TGraph *graphFase = new TGraph(nPoints, xFase, yFase);
  graphFase->SetTitle("Incertezze sulla misura della Fase con Extract Single Tone");
  graphFase->GetXaxis()->SetTitle("Frequenza Generatore");
  graphFase->GetYaxis()->SetTitle("Incertezza su Fase");
  graphFase->SetMarkerStyle(kFullCircle);
  graphFase->SetMarkerColor(kBlue);

  TF1 *f2 = new TF1("f2", "[0] + [1] * x", 1E3, 20E3);
  f2->SetParNames("q", "m");
  f2->SetLineColor(kRed);
  f2->SetLineWidth(2);
  graphFase->Fit(f2, "R");

  // PLOT ISTOGRAMMI
  TCanvas *cOndaQuadra = new TCanvas;
  cOndaQuadra->cd();
  histoOndaQuadra->Draw();

  TCanvas *cAmpiezza = new TCanvas;
  cAmpiezza->Divide(2, 3);
  cAmpiezza->cd(1);
  histo1k->Draw();
  cAmpiezza->cd(2);
  histo4k->Draw();
  cAmpiezza->cd(3);
  histo10k->Draw();
  cAmpiezza->cd(4);
  histo15k->Draw();
  cAmpiezza->cd(5);
  histo20k->Draw();

  TCanvas *cFase = new TCanvas;
  cFase->Divide(2, 3);
  cFase->cd(1);
  histoFase1k->Draw();
  cFase->cd(2);
  histoFase4k->Draw();
  cFase->cd(3);
  histoFase10k->Draw();
  cFase->cd(4);
  histoFase15k->Draw();
  cFase->cd(5);
  histo20k->Draw();

  TCanvas *cErroreAmpiezza = new TCanvas;
  cErroreAmpiezza->cd();
  graphAmpiezza->Draw("APE");
  // STAMPO

  TCanvas *cErroreFase = new TCanvas;
  cErroreFase->cd();
  graphFase->Draw("APE");
  // STAMPO
}

void amplitude_sweep()
{
  // ***** LEGGO DATI INPUT *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_resistenza.txt", "%lg %lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_induttanza.txt", "%lg %lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kAzure);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_ampiezza/sweep_freq_condensatore.txt", "%lg %lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kAzure);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; Frequency (Hz); Amplitude (V)");
  graphTotale->SetMarkerStyle(kPlus);
  graphTotale->SetMarkerColor(kPlus);
  graphTotale->SetFillColor(0);

  // ***** FIT SULLA RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_freq_resistenza, 2E3, 6E3, 5);
  funcResistenza->SetParameters(V0_mis, R_mis, L_mis, C_mis, R_tot);
  // funcResistenza->SetParLimits(0, 3.20, 3.22);
  funcResistenza->SetParLimits(1, 145, 155);
  funcResistenza->SetParLimits(3, 155E-9, 159E-9);
  funcResistenza->SetParNames("V0", "R", "L", "C", "Rtot");
  funcResistenza->SetLineWidth(2);
  funcResistenza->SetLineColor(kRed);
  graphResistenza->Fit(funcResistenza, "R");

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_freq_induttanza, 3E3, 10E3, 4);
  funcInduttanza->SetParameters(V0_mis, R_mis, L_mis, C_mis);
  funcInduttanza->SetParNames("V0", "R", "L", "C");
  funcInduttanza->SetParLimits(3, 155E-9, 159E-9);
  funcInduttanza->SetLineWidth(2);
  funcInduttanza->SetLineColor(kRed);
  graphInduttanza->Fit(funcInduttanza, "R");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcResistenza", amp_freq_condensatore, 1.5E3, 5E3, 4);
  funcCondensatore->SetParameters(V0_mis, R_mis, L_mis, C_mis);
  funcCondensatore->SetParNames("V0", "R", "L", "C");
  funcCondensatore->SetParLimits(3, 155E-9, 159E-9);
  funcCondensatore->SetLineWidth(2);
  funcCondensatore->SetLineColor(kRed);
  graphCondensatore->Fit(funcCondensatore, "R");

  // ***** FIT SU GENERATORE *****

  // ***** PLOTTO GRAFICI *****
  TCanvas *cResistenza = new TCanvas();
  graphResistenza->Draw("APE"); // ALP
  TCanvas *cInduttanza = new TCanvas();
  graphInduttanza->Draw("APE"); // ALP
  TCanvas *cCondensatore = new TCanvas();
  graphCondensatore->Draw("APE"); // ALP
  TCanvas *cTotale = new TCanvas();
  graphTotale->Draw("APE");

  TCanvas *multiCanvas = new TCanvas();
  multiCanvas->cd();
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();

/*
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
*/
}

// CHE ERRORE SU Y ASSOCIARE QUI?
void phase_sweep()
{
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_fase/sweep_phase_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Phase (RAD)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_fase/sweep_phase_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Phase (RAD)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_fase/sweep_phase_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Phase (RAD)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_fase/sweep_phase_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; Frequency (Hz); Phase (RAD)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  // ***** CREO LE FUNZIONI DI FIT *****
  TF1 *funcResistenza = new TF1("funcResistenza", phase_freq_resistenza, 0, 2E4, 3);     // LIMITI
  TF1 *funcInduttanza = new TF1("funcInduttanza", phase_freq_induttanza, 0, 2E4, 3);     // LIMITI
  TF1 *funcCondensatore = new TF1("funcResistenza", phase_freq_condensatore, 0, 2E4, 3); // LIMITI

  funcResistenza->SetParameters(R_mis, L_mis, C_mis);
  funcResistenza->SetParNames("R", "L", "C");
  funcResistenza->SetLineColor(kRed);
  funcResistenza->SetLineStyle(2);

  funcInduttanza->SetParameters(R_mis, L_mis, C_mis);
  funcInduttanza->SetParNames("R", "L", "C");
  funcInduttanza->SetLineColor(kRed);
  funcInduttanza->SetLineStyle(2);

  funcCondensatore->SetParameters(R_mis, L_mis, C_mis);
  funcCondensatore->SetParNames("R", "L", "C");
  funcCondensatore->SetLineColor(kRed);
  funcCondensatore->SetLineStyle(2);
  // MANCA IL TOTALE, POI VEDIAMO COME SI FA

  graphResistenza->Fit(funcResistenza, "R");     // OPZIONI R
  graphInduttanza->Fit(funcInduttanza, "R");     // OPZIONI R
  graphCondensatore->Fit(funcCondensatore, "R"); // OPZIONI R
  // MANCA FIT CON POL0 DEL TOTALE
  // ***** FINE PARTE FIT *****

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
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Phase Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
  // Vedi cosa fa
}

void amplitude_time_sotto_risonanza()
{
  constexpr Double_t f_mis = 0.; // INSERISCIIII

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  // ***** CREO LE FUNZIONI DI FIT *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.03, 5);       // LIMITI
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.03, 5);       // LIMITI
  TF1 *funcCondensatore = new TF1("funcResistenza", amp_time_condensatore, 0.03, 2E4, 5); // LIMITI

  funcResistenza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcResistenza->SetParNames("V0", "R", "L", "C", "f");
  funcResistenza->SetLineColor(kRed);
  funcResistenza->SetLineStyle(2);

  funcInduttanza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcInduttanza->SetParNames("V0", "R", "L", "C", "f");
  funcInduttanza->SetLineColor(kRed);
  funcInduttanza->SetLineStyle(2);

  funcCondensatore->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcCondensatore->SetParNames("V0", "R", "L", "C", "f");
  funcCondensatore->SetLineColor(kRed);
  funcCondensatore->SetLineStyle(2);
  // MANCA IL TOTALE, POI VEDIAMO COME SI FA

  graphResistenza->Fit(funcResistenza, "R");     // OPZIONI R
  graphInduttanza->Fit(funcInduttanza, "R");     // OPZIONI R
  graphCondensatore->Fit(funcCondensatore, "R"); // OPZIONI R
  // MANCA FIT CON POL0 DEL TOTALE
  // ***** FINE PARTE FIT *****

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
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Ampiezze e Tempo - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
  // Vedi cosa fa
}

void amplitude_time_in_risonanza()
{
  constexpr Double_t f_mis = 3530; // INSERISCIIII

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/in_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/in_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/in_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/in_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  // ***** CREO LE FUNZIONI DI FIT *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.03, 5);       // LIMITI
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.03, 5);       // LIMITI
  TF1 *funcCondensatore = new TF1("funcResistenza", amp_time_condensatore, 0.03, 2E4, 5); // LIMITI

  funcResistenza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcResistenza->SetParNames("V0", "R", "L", "C", "f");
  funcResistenza->SetLineColor(kRed);
  funcResistenza->SetLineStyle(2);

  funcInduttanza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcInduttanza->SetParNames("V0", "R", "L", "C", "f");
  funcInduttanza->SetLineColor(kRed);
  funcInduttanza->SetLineStyle(2);

  funcCondensatore->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcCondensatore->SetParNames("V0", "R", "L", "C", "f");
  funcCondensatore->SetLineColor(kRed);
  funcCondensatore->SetLineStyle(2);
  // MANCA IL TOTALE, POI VEDIAMO COME SI FA

  graphResistenza->Fit(funcResistenza, "R");     // OPZIONI R
  graphInduttanza->Fit(funcInduttanza, "R");     // OPZIONI R
  graphCondensatore->Fit(funcCondensatore, "R"); // OPZIONI R
  // MANCA FIT CON POL0 DEL TOTALE
  // ***** FINE PARTE FIT *****

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
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Ampiezze e Tempo - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
  // Vedi cosa fa
}

void amplitude_time_sopra_risonanza()
{
  constexpr Double_t f_mis = 0.; // INSERISCIIII

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  // ***** CREO LE FUNZIONI DI FIT *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.03, 5);       // LIMITI
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.03, 5);       // LIMITI
  TF1 *funcCondensatore = new TF1("funcResistenza", amp_time_condensatore, 0.03, 2E4, 5); // LIMITI

  funcResistenza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcResistenza->SetParNames("V0", "R", "L", "C", "f");
  funcResistenza->SetLineColor(kRed);
  funcResistenza->SetLineStyle(2);

  funcInduttanza->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcInduttanza->SetParNames("V0", "R", "L", "C", "f");
  funcInduttanza->SetLineColor(kRed);
  funcInduttanza->SetLineStyle(2);

  funcCondensatore->SetParameters(V0_mis, R_mis, L_mis, C_mis, f_mis);
  funcCondensatore->SetParNames("V0", "R", "L", "C", "f");
  funcCondensatore->SetLineColor(kRed);
  funcCondensatore->SetLineStyle(2);
  // MANCA IL TOTALE, POI VEDIAMO COME SI FA

  graphResistenza->Fit(funcResistenza, "R");     // OPZIONI R
  graphInduttanza->Fit(funcInduttanza, "R");     // OPZIONI R
  graphCondensatore->Fit(funcCondensatore, "R"); // OPZIONI R
  // MANCA FIT CON POL0 DEL TOTALE
  // ***** FINE PARTE FIT *****

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
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Ampiezze e Tempo - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
  // Vedi cosa fa
}