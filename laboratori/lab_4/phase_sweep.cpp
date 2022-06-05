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

//////////////////////////////////
//                              //
// !!!! OPERAZIONI DA FARE !!!! //
//                              //
// FONT SX = 20                 //
// FONT DX = 16                 //
//                              //
//////////////////////////////////

// VEDI QUESTI PARAMETRI
constexpr Double_t V0_mis = 5;
constexpr Double_t L_mis = 11.46 * 1E-3;

constexpr Double_t R_mis = 150.47;
constexpr Double_t R_agg = 20;
constexpr Double_t R_tot = R_mis + 50. + R_agg;

constexpr Double_t C_mis = 157.8 * 1E-9;
constexpr Double_t C_tot = 177 * 1E-9;
constexpr Double_t C_agg = C_tot - C_mis;

constexpr Double_t width = 1440;
constexpr Double_t height = 720;

void setStyle()
{
  gROOT->SetStyle("BELLE2"); // BELLE2!!!!!!!!
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetLineScalePS(1);
}

void setGraphicsGraph(TGraphErrors *graph)
{
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kAzure - 1);
  graph->SetLineColor(kAzure - 1);
  graph->SetMarkerSize(0.9);
}

void setGraphicsMultiplot(TMultiGraph *graph)
{
  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.2);
  asseX->SetTitleSize(0.04);
  asseX->SetTitleFont(2);

  asseX->SetLabelFont(1);
  asseX->SetLabelSize(0.04);
  asseX->SetLabelOffset(0.01);

  asseX->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseX->SetTickLength(0.03);

  // ***** ASSE Y *****
  TAxis *asseY = graph->GetYaxis();
  asseY->CenterTitle();
  asseY->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseY->SetNoExponent(); // no exp su assi
  asseY->SetDecimals();

  asseY->SetTitleOffset(1.3);
  asseY->SetTitleSize(0.04);
  asseY->SetTitleFont(1);

  asseY->SetLabelFont(2);
  asseY->SetLabelSize(0.04);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}

void setGraphicsFit(TF1 *func)
{
  func->SetLineWidth(3);
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

void phase_sweep()
{
  setStyle();

  // ***** GRAFICI FIT *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_fase/sweep_phase_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Fase resistenza; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphResistenza);

  TF1 *funcResistenza = new TF1("funcResistenza", phase_freq_resistenza, 500, 5.5E3, 3); // LIMITI
  funcResistenza->SetParameters(R_tot, L_mis, C_tot);
  funcResistenza->SetParNames("R", "L", "C");
  setGraphicsFit(funcResistenza);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_fase/sweep_phase_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Fase induttanza; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphInduttanza);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_fase/sweep_phase_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Fase condensatore; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphCondensatore);

  // ***** MULTIPLOT FIT *****
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Phase Sweep - Offset");
  multiGraph->SetTitle("Fasi della tensione - Offset; Frequenza (Hz); Fase (rad)");
  multiGraph->Add(graphResistenza);
  // graphResistenza->SetLineColor(kPink + 1);
  // graphResistenza->SetMarkerColor(kPink + 1);
  multiGraph->Add(graphInduttanza);
  graphInduttanza->SetLineColor(kOrange + 1);
  graphInduttanza->SetMarkerColor(kOrange + 1);
  multiGraph->Add(graphCondensatore);
  graphCondensatore->SetLineColor(kSpring - 6);
  graphCondensatore->SetMarkerColor(kSpring - 6);
  setGraphicsMultiplot(multiGraph);
  multiGraph->GetXaxis()->SetTitleSize(0.05);
  multiGraph->GetYaxis()->SetTitleSize(0.05);
  multiGraph->GetXaxis()->SetTitleOffset(1.);
  multiGraph->GetYaxis()->SetTitleOffset(1.);
  multiGraph->GetXaxis()->SetNdivisions(520);
  multiGraph->GetXaxis()->SetLimits(0, 6100);

  // ***** GRAFICI OFFSET *****
  TGraphErrors *graphResistenzaOffset = new TGraphErrors("data/sweep_fase_traslato/resistenza.txt", "%lg %lg %lg");
  graphResistenzaOffset->SetTitle("Fase resistenza; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphResistenzaOffset);

  TGraphErrors *graphInduttanzaOffset = new TGraphErrors("data/sweep_fase_traslato/induttanza.txt", "%lg %lg %lg");
  graphInduttanzaOffset->SetTitle("Fase induttanza; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphInduttanzaOffset);

  TGraphErrors *graphCondensatoreOffset = new TGraphErrors("data/sweep_fase_traslato/condensatore.txt", "%lg %lg %lg");
  graphCondensatoreOffset->SetTitle("Fase condensatore; Frequency (Hz); Phase (RAD)");
  setGraphicsGraph(graphCondensatoreOffset);

  // ***** MULTIPLOT OFFSET *****
  TMultiGraph *multiGraphOffset = new TMultiGraph("multiGraphOffset", "Phase Sweep - Offset");
  multiGraphOffset->SetTitle("Fasi della tensione - Offset; Frequenza (Hz); Fase (rad)");
  multiGraphOffset->Add(graphResistenzaOffset);
  // graphResistenzaOffset->SetLineColor(kPink + 1);
  // graphResistenzaOffset->SetMarkerColor(kPink + 1);
  multiGraphOffset->Add(graphInduttanzaOffset);
  graphInduttanzaOffset->SetLineColor(kOrange + 1);
  graphInduttanzaOffset->SetMarkerColor(kOrange + 1);
  multiGraphOffset->Add(graphCondensatoreOffset);
  graphCondensatoreOffset->SetLineColor(kSpring - 6);
  graphCondensatoreOffset->SetMarkerColor(kSpring - 6);
  setGraphicsMultiplot(multiGraphOffset);
  multiGraphOffset->GetXaxis()->SetLabelSize(0.03);
  multiGraphOffset->GetYaxis()->SetLabelSize(0.03);
  multiGraphOffset->GetXaxis()->SetLimits(0, 21000);

  // ***** PLOTTO GRAFICI *****
  TCanvas *canvas = new TCanvas("canvasSweepFase", "Sweep Fase", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);

  // Creo le Pad
  TPad *padFit = new TPad("padFit", "Fit", 0., 0., 0.40, 1., kWhite);
  TPad *padOffset = new TPad("padOffset", "Offset", 0.4, 0., 1., 1., kWhite);
  padFit->Draw();
  padOffset->Draw();

  // pad FIT
  padFit->cd();
  padFit->SetGridx();
  padFit->SetGridy();
  padFit->SetFrameLineWidth(2);
  multiGraph->SetMaximum(4);
  multiGraph->SetMinimum(-3);
  multiGraph->Draw("APE");

  padFit->BuildLegend();

  TPaveText *titoloFit = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloFit);
  titoloFit->AddText("Fasi tensione misurate");
  titoloFit->Draw();

  // pad OFFSET
  padOffset->cd();
  padOffset->SetGridx();
  padOffset->SetGridy();
  padOffset->SetFrameLineWidth(2);
  multiGraphOffset->SetMaximum(2);
  multiGraphOffset->SetMinimum(-1.5);
  multiGraphOffset->Draw("APE");

  padOffset->BuildLegend();

  TPaveText *titoloOffset = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloOffset);
  titoloOffset->AddText("Fasi tensione sovrapposte");
  titoloOffset->Draw();

  canvas->Update();
}