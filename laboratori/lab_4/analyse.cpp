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
#include "TList.h"
#include "TFitResultPtr.h"
#include "TMatrixD.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"

// VEDI QUESTI PARAMETRI
constexpr Double_t V0_mis = 5;
constexpr Double_t L_mis = 11.46 * 1E-3;

constexpr Double_t R_mis = 150.47;
constexpr Double_t R_agg = 20;
constexpr Double_t R_tot = R_mis + 50. + R_agg;

constexpr Double_t C_mis = 157.8 * 1E-9;
constexpr Double_t C_tot = 177 * 1E-9;
constexpr Double_t C_agg = C_tot - C_mis;

constexpr Double_t width = 1280;
constexpr Double_t height = 720;

void setStyle()
{
  gROOT->SetStyle("Modern");
  //  VEDI TUTTE LE OPZIONI DI FONT
  // gStyle->SetStatFont()
  // gStyle->SetPalette(57); // NON CAPISCO CHE FA
  gStyle->SetOptTitle(1);
  // gStyle->SetOptStat(112211); // 1=Integral 1=Overf 1=Underf 2=RMS 2=Mean 1=Entries 1=Name
  gStyle->SetOptFit(1111); // 1=Prob 1=Chi 1=Err 1=Param
}

void setGraphicsGraph(TGraphErrors *graph)
{
  // ***** GRAFICO *****
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kAzure - 1);
  // graph->SetMarkerSize(8);
  // graph->SetFillColor(0);

  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();

  asseX->SetTitleOffset(1.4);
  asseX->SetTitleSize(0.03);
  // asseX->SetTitleColor(kBlue);
  asseX->SetTitleFont(1);

  // asseX->SetAxisColor(kRed);
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  // asseX->SetLabelColor(kGreen);
  asseX->SetLabelFont(1);
  asseX->SetLabelSize(0.025);
  asseX->SetLabelOffset(0.01);

  asseX->SetNdivisions(1010); // 10 divsioni secondarie, 30 divisioni primarie
  // asseX->SetTickSize(0.03);
  // asseX->SetTickLength(0.03);

  // ***** ASSE Y *****
  TAxis *asseY = graph->GetYaxis();

  asseY->SetTitleOffset(1.4);
  asseY->SetTitleSize(0.03);
  // asseY->SetTitleColor(kBlue);
  asseY->SetTitleFont(1);

  // asseY->SetAxisColor(kRed);
  asseY->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseY->SetNoExponent(); // no exp su assi

  // asseY->SetLabelColor(kGreen);
  asseY->SetLabelFont(1);
  asseY->SetLabelSize(0.025);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(1010); // 10 divsioni secondarie, 30 divisioni primarie
  // asseY->SetTickSize(0.03);
  // asseY->SetTickLength(0.03);
}

void setGraphicsFit(TF1 *func)
{
  // func->SetLineStyle(2);
  func->SetLineWidth(2);
  func->SetLineColor(kRed);
}

void setGraphicsCanvas(TCanvas *c)
{
  c->SetWindowSize(width + (width - c->GetWw()), height + (height - c->GetWh()));
  c->SetGridx();
  c->SetGridy();

  // c->SetFrameFillColor(21);
  // c->SetFrameLineColor(kRed);
  c->SetFrameLineWidth(2);
}

Double_t amp_time_resistenza(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];
  Double_t f = par[3];
  Double_t V = par[4];
  Double_t xx = x[0];
  Double_t denominatore = R * R + (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C)) * (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C));
  Double_t fase = TMath::ATan((1 - TMath::TwoPi() * TMath::TwoPi() * f * f * L * C) / (TMath::TwoPi() * f * R * C)) + par[5];
  Double_t result = V * R / sqrt(denominatore) * TMath::Sin(TMath::TwoPi() * f * xx + fase);
  return result;
}

Double_t amp_time_induttanza(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];
  Double_t f = par[3];
  Double_t V = par[4];
  Double_t xx = x[0];
  Double_t denominatore = R * R + (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C)) * (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C));
  Double_t fase = TMath::ATan((1 - TMath::TwoPi() * TMath::TwoPi() * f * f * L * C) / (TMath::TwoPi() * f * R * C)) + par[5];
  Double_t result = V * TMath::TwoPi() * f * L / sqrt(denominatore) * TMath::Sin(TMath::TwoPi() * f * xx + fase + TMath::PiOver2());
  return result;
}

Double_t amp_time_condensatore(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];
  Double_t f = par[3];
  Double_t V = par[4];
  Double_t xx = x[0];
  Double_t denominatore = R * R + (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C)) * (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C));
  Double_t fase = TMath::ATan((1 - TMath::TwoPi() * TMath::TwoPi() * f * f * L * C) / (TMath::TwoPi() * f * R * C)) + par[5];
  Double_t result = V / (TMath::TwoPi() * f * C) / sqrt(denominatore) * TMath::Sin(TMath::TwoPi() * f * xx + fase - TMath::PiOver2());
  return result;
}

Double_t amp_time_totale(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];
  Double_t f = par[3];
  Double_t r = 50.;
  Double_t V = 5.;
  Double_t xx = x[0];
  Double_t numeratore = (R - r) * (R - r) + (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C)) * (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C));
  Double_t denominatore = R * R + (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C)) * (TMath::TwoPi() * f * L - 1 / (TMath::TwoPi() * f * C));
  Double_t result = V * sqrt(numeratore / denominatore) * TMath::Sin(TMath::TwoPi() * f * xx);
  return result;
}

Double_t amp_freq_resistenza(Double_t *x, Double_t *par)
{
  // 4 parametri
  Double_t L = par[1];
  Double_t C = par[2];
  Double_t R = par[3];
  Double_t Rtot = R + 50. + par[0];

  Double_t xx = x[0];

  Double_t denominatore = Rtot * Rtot + (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C)) * (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C));
  Double_t result = 5. * R / sqrt(denominatore);
  return result;
}

Double_t amp_freq_induttanza(Double_t *x, Double_t *par)
{
  // 3 parametri
  Double_t Rtot = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t denominatore = Rtot * Rtot + (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C)) * (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C));
  Double_t result = TMath::TwoPi() * xx * L * 5.0 / sqrt(denominatore);
  return result;
}

Double_t amp_freq_condensatore(Double_t *x, Double_t *par)
{
  Double_t Rtot = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t denominatore = Rtot * Rtot + (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C)) * (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C));
  Double_t result = 5.0 / (TMath::TwoPi() * xx * C) / sqrt(denominatore);
  return result;
}

Double_t amp_freq_totale(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t numeratore = (R - 50.) * (R - 50.) + (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C)) * (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C));
  Double_t denominatore = R * R + (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C)) * (TMath::TwoPi() * xx * L - 1 / (TMath::TwoPi() * xx * C));
  Double_t result = 5.0 * sqrt(numeratore / denominatore);
  return result;
}

Double_t phase_freq_resistenza(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t fase = (1 - TMath::TwoPi() * TMath::TwoPi() * xx * xx * L * C) / (R * TMath::TwoPi() * xx * C);
  Double_t result = TMath::ATan(fase);
  return result;
}

Double_t phase_freq_induttanza(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t fase = (1 - TMath::TwoPi() * TMath::TwoPi() * xx * xx * L * C) / (R * TMath::TwoPi() * xx * C);
  Double_t result = TMath::ATan(fase) + TMath::PiOver2();
  return result;
}

Double_t phase_freq_condensatore(Double_t *x, Double_t *par)
{
  Double_t R = par[0];
  Double_t L = par[1];
  Double_t C = par[2];

  Double_t xx = x[0];

  Double_t fase = (1 - TMath::TwoPi() * TMath::TwoPi() * xx * xx * L * C) / (R * TMath::TwoPi() * xx * C);
  Double_t result = TMath::ATan(fase) - TMath::PiOver2();
  return result;
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
  setStyle();
  // ***** LEGGO DATI INPUT *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Amplitude (V)");
  setGraphicsGraph(graphResistenza);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Amplitude (V)");
  setGraphicsGraph(graphInduttanza);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_ampiezza/sweep_freq_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Amplitude (V)");
  setGraphicsGraph(graphCondensatore);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; Frequency (Hz); Amplitude (V)");
  setGraphicsGraph(graphTotale);

  // ***** FIT SULLA RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_freq_resistenza, 2E3, 6E3, 4);
  funcResistenza->SetParameters(R_agg, L_mis, C_mis, R_mis);
  funcResistenza->SetParNames("R_agg", "L", "C", "R");
  setGraphicsFit(funcResistenza);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_freq_induttanza, 3E3, 10E3, 3);
  funcInduttanza->SetParameters(R_tot, L_mis, C_tot);
  funcInduttanza->SetParNames("R_tot", "L", "C_tot");
  setGraphicsFit(funcInduttanza);
  TFitResultPtr rInduttanza = graphInduttanza->Fit(funcInduttanza, "REMSQ");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcCondensatore", amp_freq_condensatore, 2E3, 4E3, 3);
  funcCondensatore->SetParameters(R_tot, L_mis, C_mis);
  funcCondensatore->SetParNames("R_tot", "L", "C_mis");
  setGraphicsFit(funcCondensatore);
  TFitResultPtr rCondensatore = graphCondensatore->Fit(funcCondensatore, "REMSQ");

  // ***** FIT SU GENERATORE *****
  TF1 *funcTotale = new TF1("funcTotale", amp_freq_totale, 2E3, 6E3, 3);
  funcTotale->SetParameters(R_tot, L_mis, C_mis);
  funcTotale->SetParNames("Rtot", "L", "C");
  setGraphicsFit(funcTotale);
  TFitResultPtr rTotale = graphTotale->Fit(funcTotale, "REMSQ");

  // ***** FATTORE DI QUALITA' *****
  Double_t maxResistenza = funcResistenza->GetMaximum();
  Double_t f_min = funcResistenza->GetX(maxResistenza / sqrt(2), 1E3, 3E3);
  Double_t f_max = funcResistenza->GetX(maxResistenza / sqrt(2), 4E3, 7E3);
  Double_t f_centro = funcResistenza->GetMaximumX(3E3, 4E3);
  Double_t Q_resistenza = f_centro / (f_max - f_min);

  std::cout << "***** CALCOLO FATTORE DI QUALITA' *****" << '\n'
            << "Massimo: " << maxResistenza << '\n'
            << "F MIN: " << f_min << '\n'
            << "F MAX: " << f_max << '\n'
            << "F0: " << f_centro << '\n'
            << "FATTORE DI QUALITA'. Q = " << Q_resistenza << '\n';

  // ***** FREQUENZA RISONANZA *****
  Double_t fRisResistenza = funcResistenza->GetMaximumX(3E3, 4E3);
  Double_t fRisInduttanza = funcInduttanza->GetMaximumX(3E3, 4E3);
  Double_t fRisCondensatore = funcCondensatore->GetMaximumX(3E3, 4E3);
  Double_t fRisTotale = funcTotale->GetMinimumX(3E3, 4E3);

  std::cout << "***** CALCOLO FREQUENZA RISONANZA *****" << '\n'
            << "Da Resistenza: " << fRisResistenza << '\n'
            << "Da Induttanza: " << fRisInduttanza << '\n'
            << "Da Condensatore: " << fRisCondensatore << '\n'
            << "Da Totale: " << fRisTotale << '\n';

  // ***** PLOTTO GRAFICI *****
  TCanvas *cResistenza = new TCanvas("cResistenza", "Sweep Ampiezza Resistenza", width, height);
  setGraphicsCanvas(cResistenza);
  // cResistenza->SetLogx();
  graphResistenza->Draw("ALP"); // L=polyline C=SmoothCurve E=ErrorBar

  // TPaveText *boxResistenza = new TPaveText(1., 1., .7, .7, "NDC, NB"); // NDC=CoordinateRelative NB=noBorders
  // boxResistenza->AddText("A TPaveText can contain severals line of text.");
  // boxResistenza->AddText("They are added to the pave using the AddText method.");
  // boxResistenza->AddLine(.0, .5, 1., .5);
  // boxResistenza->AddText("Even complex TLatex formulas can be added:");
  // TText *t1 = boxResistenza->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
  // t1->SetTextColor(kBlue);
  // TText *t2 = boxResistenza->GetLineWith("Even");
  // t2->SetTextColor(kOrange + 1);
  // boxResistenza->Draw();

  TCanvas *cInduttanza = new TCanvas("cInduttanza", "Sweep Ampiezza Induttanza", width, height);
  setGraphicsCanvas(cInduttanza);
  graphInduttanza->Draw("ALP"); // L=polyline C=SmoothCurve E=ErrorBar

  TCanvas *cCondensatore = new TCanvas("cCondensatore", "Sweep Ampiezza Condensatore", width, height);
  setGraphicsCanvas(cCondensatore);
  graphCondensatore->Draw("ALP"); // L=polyline C=SmoothCurve E=ErrorBar

  TCanvas *cTotale = new TCanvas("cTotale", "Sweep Ampiezza Totale", width, height);
  setGraphicsCanvas(cTotale);
  graphTotale->Draw("ALP"); // L=polyline C=SmoothCurve E=ErrorBar

  // TCanvas *multiCanvas = new TCanvas();
  // multiCanvas->cd();
  // setGraphicsCanvas(multiCanvas);
  // TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep - Risultati finali");
  // multiGraph->Add(graphResistenza);
  // multiGraph->Add(graphInduttanza);
  // multiGraph->Add(graphCondensatore);
  // multiGraph->Add(graphTotale);
  // multiGraph->Draw("ALP"); // COSA FA LP?
  // multiCanvas->BuildLegend();
}

void phase_sweep()
{
  /////////////////////////////////////////////////////////////////////
  //    ATTENZIONE!                                                  //
  //                                                                 //
  //    VEDI LA FASE DEL GENERATORE, NON E' FISSA A ZERO, SI PUO'    //
  //    ELIMINARE OFFSET SULLE MISURE SOTTRAENDO AI GRAFICI          //
  //    I VALORI DELLA FASE DEL GENERATORE                           //
  //                                                                 //
  //    NOTA                                                         //
  //                                                                 //
  //    VEDI PERCHE? E' COSI'                                        //
  //                                                                 //
  ////////////////////////////////////////////////////////////////////

  // ***** LEGGO DATI INPUT *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_fase/sweep_phase_resistenza.txt", "%lg %lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Phase (RAD)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_fase/sweep_phase_induttanza.txt", "%lg %lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Phase (RAD)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kAzure);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_fase/sweep_phase_condensatore.txt", "%lg %lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Phase (RAD)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kAzure);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_fase/sweep_phase_totale.txt", "%lg %lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; Frequency (Hz); Phase (RAD)");
  graphTotale->SetMarkerStyle(kPlus);
  graphTotale->SetMarkerColor(kAzure);
  graphTotale->SetFillColor(0);

  // ***** FIT SU RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", phase_freq_resistenza, 1E3, 6E3, 3); // LIMITI
  funcResistenza->SetParameters(R_tot, L_mis, C_tot);
  funcResistenza->SetParNames("R", "L", "C");
  funcResistenza->SetLineWidth(2);
  funcResistenza->SetLineColor(kRed);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");
  TMatrixD covResistenza = rResistenza->GetCovarianceMatrix();
  std::cout << '\n'
            << "***** MATRICE COVARIANZA FIT RESISTENZA *****" << '\n';
  Double_t ChiResistenza = funcResistenza->GetChisquare();
  Double_t NdofResistenza = funcResistenza->GetNDF();
  std::cout << "ChiSquare: " << ChiResistenza << '\n'
            << "NDF: " << NdofResistenza << '\n';
  covResistenza.Print();

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", phase_freq_induttanza, 500, 5E3, 3); // LIMITI
  funcInduttanza->SetParameters(R_tot, L_mis, C_tot);
  funcInduttanza->SetParNames("R", "L", "C");
  funcInduttanza->SetLineWidth(2);
  funcInduttanza->SetLineColor(kRed);
  graphInduttanza->Fit(funcInduttanza, "REMSQ");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcResistenza", phase_freq_condensatore, 500, 5.5E3, 3); // LIMITI
  funcCondensatore->SetParameters(R_tot, L_mis, C_tot);
  funcCondensatore->SetParNames("R", "L", "C");
  funcCondensatore->SetLineWidth(2);
  funcCondensatore->SetLineColor(kRed);
  TFitResultPtr rCondensatore = graphCondensatore->Fit(funcCondensatore, "REMSQ");
  TMatrixD covCondensatore = rCondensatore->GetCovarianceMatrix();
  std::cout << '\n'
            << "***** MATRICE COVARIANZA FIT CONDENSATORE *****" << '\n';
  Double_t ChiCondensatore = funcCondensatore->GetChisquare();
  Double_t NdofCondensatore = funcCondensatore->GetNDF();
  std::cout << "ChiSquare: " << ChiCondensatore << '\n'
            << "NDF: " << NdofCondensatore << '\n';
  covCondensatore.Print();

  // ***** CALCOLO RISONANZA *****
  Double_t fRisResistenza = funcResistenza->GetX(0., 3E3, 4E3);
  // ROOT::Math::RootFinder finder;
  // finder.Solve(*funcResistenza, 3E3, 4E3);
  // Double_t fRisResistenza = finder.Root();
  Double_t fRisInduttanza = funcInduttanza->GetX(TMath::PiOver2(), 3E3, 4E3);
  Double_t fRisCondensatore = funcCondensatore->GetX(-TMath::PiOver2(), 3E3, 4E3);

  std::cout << '\n'
            << "***** FREQUENZA RISONANZA *****" << '\n'
            << "Da Resistenza: " << fRisResistenza << '\n'
            << "Da Induttanza: " << fRisInduttanza << '\n'
            << "Da Condensatore: " << fRisCondensatore << '\n';

  // ***** PLOTTO GRAFICI *****
  TCanvas *cResistenza = new TCanvas();
  graphResistenza->Draw("APE"); // ALP
  TCanvas *cInduttanza = new TCanvas();
  graphInduttanza->Draw("APE"); // ALP
  TCanvas *cCondensatore = new TCanvas();
  graphCondensatore->Draw("APE"); // ALP
  // TCanvas *cTotale = new TCanvas();
  // graphTotale->Draw("APE");

  TCanvas *multiCanvas = new TCanvas();
  multiCanvas->cd();
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Phase Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
}

void phase_offset()
{
  // ***** LEGGO DATI INPUT *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_fase_traslato/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Phase (RAD)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_fase_traslato/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Phase (RAD)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kRed);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_fase_traslato/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Phase (RAD)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kGreen);
  graphCondensatore->SetFillColor(0);

  // TLine *line1 = new TLine(500, 0, 500, 1);

  // ***** MULTIPLOT *****
  TCanvas *multiCanvas = new TCanvas();
  multiCanvas->cd();
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Phase Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Draw("ALP");
  multiCanvas->BuildLegend();
}

void amplitude_time_sotto_risonanza() // FREQUENZA = 2 KHz
{
  constexpr Double_t f_mis = 2E3;

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kAzure);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kAzure);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kPlus);
  graphTotale->SetMarkerColor(kAzure);
  graphTotale->SetFillColor(0);

  // ***** FIT SU RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.00479, 6);
  funcResistenza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcResistenza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcResistenza->SetNpx(10000);
  funcResistenza->FixParameter(3, 2E3);
  funcResistenza->FixParameter(4, 5.0);
  funcResistenza->SetLineWidth(2);
  funcResistenza->SetLineColor(kRed);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.00479, 6);
  funcInduttanza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcInduttanza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcInduttanza->SetNpx(10000);
  funcInduttanza->FixParameter(3, 2E3);
  funcInduttanza->FixParameter(4, 5.0);
  funcInduttanza->SetLineWidth(2);
  funcInduttanza->SetLineColor(kRed);
  TFitResultPtr rInduttanza = graphInduttanza->Fit(funcInduttanza, "REMSQ");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcCondensatore", amp_time_condensatore, 0, 0.00479, 6);
  funcCondensatore->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcCondensatore->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcCondensatore->SetNpx(10000);
  funcCondensatore->FixParameter(3, 2E3);
  funcCondensatore->FixParameter(4, 5.0);
  funcCondensatore->SetLineWidth(2);
  funcCondensatore->SetLineColor(kRed);
  TFitResultPtr rCondensatore = graphCondensatore->Fit(funcCondensatore, "REMSQ");

  // ***** FIT SU TOTALE *****
  TF1 *funcTotale = new TF1("funcTotale", amp_time_totale, 0, 0.00479, 4);
  funcTotale->SetParameters(R_mis, L_mis, C_mis, f_mis);
  funcTotale->SetParNames("R", "L", "C", "f");
  funcTotale->SetNpx(10000);
  funcTotale->FixParameter(3, 2E3);
  funcTotale->SetLineWidth(2);
  funcTotale->SetLineColor(kRed);
  TFitResultPtr rTotale = graphTotale->Fit(funcTotale, "REMSQ");

  //***** PLOTTO I GRAFICI *****
  TCanvas *cResistenza = new TCanvas();
  graphResistenza->Draw("APE");
  TCanvas *cInduttanza = new TCanvas();
  graphInduttanza->Draw("APE");
  TCanvas *cCondensatore = new TCanvas();
  graphCondensatore->Draw("APE");
  TCanvas *cTotale = new TCanvas();
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
}

void amplitude_time_in_risonanza() // FREQUENZA = 3.5301 KHz
{
  constexpr Double_t f_mis = 3530.1;

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/in_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/in_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kAzure);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/in_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kAzure);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/in_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kPlus);
  graphTotale->SetMarkerColor(kAzure);
  graphTotale->SetFillColor(0);

  // ***** FIT SU RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.00279, 6);
  funcResistenza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcResistenza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcResistenza->SetNpx(10000);
  funcResistenza->FixParameter(3, 3530.1);
  funcResistenza->FixParameter(4, 5.0);
  funcResistenza->SetLineWidth(2);
  funcResistenza->SetLineColor(kRed);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.00279, 6);
  funcInduttanza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcInduttanza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcInduttanza->SetNpx(10000);
  funcInduttanza->FixParameter(3, 3530.1);
  funcInduttanza->FixParameter(4, 5.0);
  funcInduttanza->SetLineWidth(2);
  funcInduttanza->SetLineColor(kRed);
  TFitResultPtr rInduttanza = graphInduttanza->Fit(funcInduttanza, "REMSQ");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcCondensatore", amp_time_condensatore, 0, 0.00279, 6);
  funcCondensatore->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcCondensatore->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcCondensatore->SetNpx(10000);
  funcCondensatore->FixParameter(3, 3530.1);
  funcCondensatore->FixParameter(4, 5.0);
  funcCondensatore->SetLineWidth(2);
  funcCondensatore->SetLineColor(kRed);
  TFitResultPtr rCondensatore = graphCondensatore->Fit(funcCondensatore, "REMSQ");

  // ***** FIT SU TOTALE *****
  TF1 *funcTotale = new TF1("funcTotale", amp_time_totale, 0, 0.00279, 4);
  funcTotale->SetParameters(R_mis, L_mis, C_mis, f_mis);
  funcTotale->SetParNames("R", "L", "C", "f");
  funcTotale->SetNpx(10000);
  funcTotale->FixParameter(3, 3530.1);
  funcTotale->SetLineWidth(2);
  funcTotale->SetLineColor(kRed);
  TFitResultPtr rTotale = graphTotale->Fit(funcTotale, "REMSQ");

  //***** PLOTTO I GRAFICI *****
  TCanvas *cResistenza = new TCanvas();
  graphResistenza->Draw("APE");
  TCanvas *cInduttanza = new TCanvas();
  graphInduttanza->Draw("APE");
  TCanvas *cCondensatore = new TCanvas();
  graphCondensatore->Draw("APE");
  TCanvas *cTotale = new TCanvas();
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
}

void amplitude_time_sopra_risonanza() // FREQUENZA 10kHz
{
  constexpr Double_t f_mis = 10000;

  TGraphErrors *graphResistenza = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Ampiezza Resistenza; time (s); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kPlus);
  graphResistenza->SetMarkerColor(kAzure);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Ampiezza Induttanza; time (s); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kPlus);
  graphInduttanza->SetMarkerColor(kAzure);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Ampiezza Condensatore; time (s); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kPlus);
  graphCondensatore->SetMarkerColor(kAzure);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Ampiezza Totale; time (s); Amplitude (V)");
  graphTotale->SetMarkerStyle(kPlus);
  graphTotale->SetMarkerColor(kAzure);
  graphTotale->SetFillColor(0);

  // ***** FIT SU RESISTENZA *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.00119, 6);
  funcResistenza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcResistenza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcResistenza->SetNpx(10000);
  funcResistenza->FixParameter(3, 10000);
  funcResistenza->FixParameter(4, 5.0);
  funcResistenza->SetLineWidth(2);
  funcResistenza->SetLineColor(kRed);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  // ***** FIT SU INDUTTANZA *****
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_time_induttanza, 0, 0.00119, 6);
  funcInduttanza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcInduttanza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcInduttanza->SetNpx(10000);
  funcInduttanza->FixParameter(3, 10000);
  funcInduttanza->FixParameter(4, 5.0);
  funcInduttanza->SetLineWidth(2);
  funcInduttanza->SetLineColor(kRed);
  TFitResultPtr rInduttanza = graphInduttanza->Fit(funcInduttanza, "REMSQ");

  // ***** FIT SU CONDENSATORE *****
  TF1 *funcCondensatore = new TF1("funcCondensatore", amp_time_condensatore, 0, 0.00119, 6);
  funcCondensatore->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcCondensatore->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcCondensatore->SetNpx(10000);
  funcCondensatore->FixParameter(3, 10000);
  funcCondensatore->FixParameter(4, 5.0);
  funcCondensatore->SetLineWidth(2);
  funcCondensatore->SetLineColor(kRed);
  TFitResultPtr rCondensatore = graphCondensatore->Fit(funcCondensatore, "REMSQ");

  // ***** FIT SU TOTALE *****
  TF1 *funcTotale = new TF1("funcTotale", amp_time_totale, 0, 0.00119, 4);
  funcTotale->SetParameters(R_mis, L_mis, C_mis, f_mis);
  funcTotale->SetParNames("R", "L", "C", "f");
  funcTotale->SetNpx(10000);
  funcTotale->FixParameter(3, 10000);
  funcTotale->SetLineWidth(2);
  funcTotale->SetLineColor(kRed);
  TFitResultPtr rTotale = graphTotale->Fit(funcTotale, "REMSQ");

  //***** PLOTTO I GRAFICI *****
  TCanvas *cResistenza = new TCanvas();
  graphResistenza->Draw("APE");
  TCanvas *cInduttanza = new TCanvas();
  graphInduttanza->Draw("APE");
  TCanvas *cCondensatore = new TCanvas();
  graphCondensatore->Draw("APE");
  TCanvas *cTotale = new TCanvas();
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
}