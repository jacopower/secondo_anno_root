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

constexpr Double_t width = 1280;
constexpr Double_t height = 720;

void setStyle()
{
   gROOT->SetStyle("BELLE2"); // BELLE2!!!!!!!!
  // gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPadRightMargin(0.03);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineScalePS(1);
}

void setGraphicsGraph(TGraph *h)
{
  h->SetMarkerStyle(kFullCircle);
  h->SetMarkerColor(kAzure - 1);
  h->SetLineColor(kAzure - 1);
  h->SetMarkerSize(0.9);

  // *****  ASSE X *****
  TAxis *asseX = h->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.2);
  asseX->SetTitleSize(0.04);

  asseX->SetLabelSize(0.03);
  asseX->SetLabelOffset(0.013);

  asseX->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseX->SetTickLength(0.03);

  // ***** ASSE Y *****
  TAxis *asseY = h->GetYaxis();
  asseY->CenterTitle();
  asseY->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseY->SetNoExponent(); // no exp su assi
  asseY->SetDecimals();

  asseY->SetTitleOffset(1.);
  asseY->SetTitleSize(0.04);

  asseY->SetLabelSize(0.03);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}

void setGraphicsHisto(TH1 *h)
{
  h->SetMarkerStyle(kFullCircle);
  h->SetMarkerColor(kAzure - 1);
  h->SetLineColor(kAzure - 1);
  h->SetMarkerSize(0.9);

  // *****  ASSE X *****
  TAxis *asseX = h->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.2);
  asseX->SetTitleSize(0.04);

  asseX->SetLabelSize(0.03);
  asseX->SetLabelOffset(0.013);

  asseX->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseX->SetTickLength(0.03);

  // ***** ASSE Y *****
  TAxis *asseY = h->GetYaxis();
  asseY->CenterTitle();
  asseY->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseY->SetNoExponent(); // no exp su assi
  asseY->SetDecimals();

  asseY->SetTitleOffset(1.);
  asseY->SetTitleSize(0.04);

  asseY->SetLabelSize(0.03);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}


void setGraphicsFit(TF1 *func)
{
  func->SetLineWidth(2);
  func->SetLineColor(kRed);
}


void setGraphicsTitolo(TPaveText *titolo)
{
  titolo->SetMargin(0);
  titolo->SetTextFont(2);
  titolo->SetShadowColor(kRed);
  titolo->SetBorderSize(2);
  titolo->SetFillColor(kWhite);
  titolo->SetTextAlign(22);
}

void setGraphicsBox(TPaveText *box)
{
  box->SetFillColor(kWhite);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
}

void rumore()
{
  setStyle();

  constexpr int N = 100;
  TH1F *histoOndaQuadra = new TH1F("rumore e risultati fit", "Rumore con Onda Quadra", N, 4.978, 4.99);
  histoOndaQuadra->GetXaxis()->SetTitle("Tensione (V)");
  histoOndaQuadra->GetYaxis()->SetTitle("Occorrenze");
  setGraphicsHisto(histoOndaQuadra);

  std::ifstream in;

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

  histoOndaQuadra->Fit("gaus");

  // PLOT ONDA QUADRA
  TCanvas *canvas = new TCanvas("Canvas Onda Quadra", "Sweep Ampiezza", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);
  TPad *pad = new TPad("pad", "Pad", 0., 0., 1., 1., kWhite);
  pad->Draw();
  pad->cd();
  pad->SetGridx();
  pad->SetGridy();
  pad->SetFrameLineWidth(2);

  histoOndaQuadra->SetMaximum(100);
  histoOndaQuadra->Draw();

  TPaveText *titoloOndaQuadra = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloOndaQuadra);
  titoloOndaQuadra->AddText("Rumore misurato");
  titoloOndaQuadra->Draw();

  TPaveText *boxOndaQuadra = new TPaveText(1., 1., .7, .7, "NDC, NB"); // NDC=CoordinateRelative NB=noBorders RB=RightBottom
  setGraphicsBox(boxOndaQuadra);
  boxOndaQuadra->AddText("Parametri istogramma:");
  boxOndaQuadra->AddText("Entries:   982");
  boxOndaQuadra->AddText("Media:    (4.98391 +/- 0.00004)V");
  boxOndaQuadra->AddText("Std Dev: (0.00144 +/- 0.00003)V");
  boxOndaQuadra->AddText("");
  boxOndaQuadra->AddText("Risultati fit:");
  boxOndaQuadra->AddText("Media:    (4.98390 +/- 0.00005)V");
  boxOndaQuadra->AddText("Std Dev: (0.00139 +/- 0.00003)V");

  boxOndaQuadra->Draw();

  canvas->Update();
}

void phase()
{
  TGraph *graph = new TGraph("data/incertezza_fase.txt");
  graph->SetTitle(0);
  graph->GetXaxis()->SetTitle("Frequenza (Hz)");
  graph->GetYaxis()->SetTitle("Incertezza");
  setGraphicsGraph(graph);

  TF1 *f1 = new TF1("f1", "[0] + [1] * x", 4E3, 20E3);
  f1->SetParameters(5.85E-4, 2.47E-7);
  f1->SetParNames("q", "m");
  setGraphicsFit(f1);
  graph->Fit(f1, "REMSQ");

  // PLOT 
  TCanvas *canvas = new TCanvas("Canvas Onda Quadra", "Sweep Ampiezza", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);
  TPad *pad = new TPad("pad", "Pad", 0., 0., 1., 1., kWhite);
  pad->Draw();
  pad->cd();
  pad->SetGridx();
  pad->SetGridy();
  pad->SetFrameLineWidth(2);

  graph->SetMaximum(0.006);
  graph->SetMinimum(0.0010);
  graph->GetXaxis()->SetLimits(0, 22000);
  graph->Draw("ALP");

  TPaveText *titolo = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titolo);
  titolo->AddText("Incertezza sulla fase - fit lineare");
  titolo->Draw();

  canvas->Update();
}

void old()
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