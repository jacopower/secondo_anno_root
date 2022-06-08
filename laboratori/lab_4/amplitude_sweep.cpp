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

// TPad *padResistenza = new TPad("padResistenza", "Resistenza", 0., 0., 0.50, 1., kRed);
//  LOWER LEFT CORNER = 5% of width from left edge | 50% of hwight from low edge
//  UPPER RIGHT CORNER = 95% of width from left edge | 95% of hwight from low edge

constexpr Double_t V0_mis = 5;
constexpr Double_t L_mis = 11.46 * 1E-3;
constexpr Double_t R_mis = 150.47;
constexpr Double_t R_agg = 20;
constexpr Double_t R_tot = R_mis + 50. + R_agg;
constexpr Double_t C_mis = 157.8 * 1E-9;
constexpr Double_t C_tot = 177 * 1E-9;
constexpr Double_t C_agg = C_tot - C_mis;

constexpr Double_t width = 1440;
constexpr Double_t height = 600;

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

void setStyle()
{
  gROOT->SetStyle("BELLE2"); // BELLE2!!!!!!!!
  //gStyle->SetPadLeftMargin(0.07);
  //gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  //gStyle->SetLineScalePS(1);
}

void setGraphicsGraph(TGraphErrors *graph)
{
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kAzure - 1);
  graph->SetLineColor(kAzure - 1);
  graph->SetMarkerSize(0.9);
  graph->SetLineWidth(2);

  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.2);
  asseX->SetTitleSize(0.04);
  //asseX->SetTitleFont(2);

  //asseX->SetLabelFont(1);
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
  //asseY->SetTitleFont(1);

  //asseY->SetLabelFont(1);
  asseY->SetLabelSize(0.04);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}

void setGraphicsMultiplot(TMultiGraph *graph)
{
  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.3);
  asseX->SetTitleSize(0.04);
  //asseX->SetTitleFont(2);

  //asseX->SetLabelFont(1);
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

  asseY->SetTitleOffset(0.6);
  asseY->SetTitleSize(0.04);
 // asseY->SetTitleFont(1);

 // asseY->SetLabelFont(1);
  asseY->SetLabelSize(0.04);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}

void setGraphicsFit(TF1 *func)
{
  func->SetLineWidth(3);
  func->SetLineColor(kRed);
  func->SetMarkerColor(kRed);
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

void amplitude_sweep()
{
  setStyle();

  // ***** RESISTENZA *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Tensione resistenza; Frequenza (Hz); Tensione (V)");
  setGraphicsGraph(graphResistenza);

  TF1 *funcResistenza = new TF1("funcResistenza", amp_freq_resistenza, 2E3, 6E3, 4);
  funcResistenza->SetParameters(R_agg, L_mis, C_mis, R_mis);
  funcResistenza->SetParNames("R_agg", "L", "C", "R");
  setGraphicsFit(funcResistenza);
  TFitResultPtr rResistenza = graphResistenza->Fit(funcResistenza, "REMSQ");

  // ***** GENERATORE *****
  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Tensione generatore; Frequenza (Hz); Tensione (V)");
  setGraphicsGraph(graphTotale);

  TF1 *funcTotale = new TF1("funcTotale", amp_freq_totale, 2E3, 6E3, 3);
  funcTotale->SetParameters(R_tot, L_mis, C_mis);
  funcTotale->SetParNames("Rtot", "L", "C");
  setGraphicsFit(funcTotale);
  TFitResultPtr rTotale = graphTotale->Fit(funcTotale, "REMSQ");

  // ***** PLOT *****
  // Canvas
  TCanvas *canvas = new TCanvas("canvasFFF", "Resistenza / Generatore -  Sweep ampiezza", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);

  // Creo le Pad
  TPad *padResistenza = new TPad("padResistenza", "Resistenza", 0., 0., 0.50, 1., kWhite);
  TPad *padGeneratore = new TPad("padGeneratore", "Generatore", 0.50, 0., 1., 1., kWhite);
  padResistenza->Draw();
  padGeneratore->Draw();

  // pad Resistenza
  padResistenza->cd();
  padResistenza->SetLogx();
  padResistenza->SetGridx();
  padResistenza->SetGridy();
  padResistenza->SetFrameLineWidth(2);
  graphResistenza->GetXaxis()->SetMoreLogLabels();
  graphResistenza->SetMinimum(0);
  graphResistenza->SetMaximum(4.);
  graphResistenza->Draw("ALP");

  auto legendResistenza = new TLegend(0.1, 0.7, 0.48, 0.9);
  legendResistenza->AddEntry(graphResistenza, "Dati sperimentali");
  legendResistenza->AddEntry(funcResistenza, "Fit", "l");
  legendResistenza->Draw();

  TPaveText *titoloResistenza = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloResistenza);
  titoloResistenza->AddText("Risposta in frequenza - resistenza");
  titoloResistenza->Draw();

  TPaveText *boxResistenza = new TPaveText(1., 1., .7, .7, "NDC, NB"); // NDC=CoordinateRelative NB=noBorders RB=RightBottom
  setGraphicsBox(boxResistenza);
  boxResistenza->AddText("Parametri Fit:");
  boxResistenza->AddText("Rtot = (230.6 +/- 3.8) Ohm");
  boxResistenza->AddText("Rmis = (148.0 +/- 2.4) Ohm");
  boxResistenza->AddText("L = (11.96 +/- 0.20) mH");
  boxResistenza->AddText("C = (170.3 +/- 2.8) nF");
  boxResistenza->Draw();

  // pad Generatore
  padGeneratore->cd();
  padGeneratore->SetGridx();
  padGeneratore->SetGridy();
  padGeneratore->SetFrameLineWidth(2);
  graphTotale->SetMinimum(3.8);
  graphTotale->SetMaximum(5.2);
  graphTotale->Draw("ALP");

  auto legendGeneratore = new TLegend(0.1, 0.7, 0.48, 0.9);
  legendGeneratore->AddEntry(graphTotale, "Dati sperimentali");
  legendGeneratore->AddEntry(funcTotale, "Fit", "l");
  legendGeneratore->Draw();

  TPaveText *titoloTotale = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloTotale);
  titoloTotale->AddText("Risposta in frequenza - generatore");
  titoloTotale->Draw();

  TPaveText *boxTotale = new TPaveText(1., 1., .7, .7, "NDC, NB"); // NDC=CoordinateRelative NB=noBorders RB=RightBottom
  setGraphicsBox(boxTotale);
  boxTotale->AddText("Resistenza dal fit:");
  boxTotale->AddText("R = (225.9 +/- 1.2) Ohm");
  boxTotale->Draw();

  canvas->Update();
}

void multiplot()
{
  setStyle();

  // ***** LEGGO GRAFICI *****
  TGraphErrors *graphResistenza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_resistenza.txt", "%lg %lg %lg");
  graphResistenza->SetTitle("Tensione resistenza; Frequenza (Hz); Tensione (V)");
  setGraphicsGraph(graphResistenza);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Tensione induttanza; Frequenza (Hz); Tensione (V)");
  setGraphicsGraph(graphInduttanza);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_ampiezza/sweep_freq_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Tensione condensatore; Frequenza (Hz); Tensione (V)");
  setGraphicsGraph(graphCondensatore);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Tensione generatore; Frequenza (Hz); Tensione (V)");
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

  // ***** MULTIPLOT FIT *****
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep");
  multiGraph->SetTitle("Risposta in frequenza - tutte le componenti; Frequenza (Hz); Tensione (V)");
  multiGraph->Add(graphResistenza);
  graphResistenza->SetLineColor(kPink + 1);
  graphResistenza->SetMarkerColor(kPink + 1);
  multiGraph->Add(graphInduttanza);
  graphInduttanza->SetLineColor(kOrange + 1);
  graphInduttanza->SetMarkerColor(kOrange + 1);
  multiGraph->Add(graphCondensatore);
  graphCondensatore->SetLineColor(kSpring - 6);
  graphCondensatore->SetMarkerColor(kSpring - 6);
  multiGraph->Add(graphTotale);

  setGraphicsMultiplot(multiGraph);

  TCanvas *canvas = new TCanvas("canvasSweepAmpiezza", "Sweep Ampiezza", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);

  TPad *pad = new TPad("pad", "Pad", 0., 0., 1., 1., kWhite);
  pad->Draw();

  pad->cd();
  pad->SetGridx();
  pad->SetGridy();
  pad->SetFrameLineWidth(2);
  multiGraph->SetMaximum(8);
  multiGraph->SetMinimum(0);
  multiGraph->GetXaxis()->SetLimits(0, 20500);
  multiGraph->Draw("APE");

  auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  legend->AddEntry(graphResistenza, "Tensione resistenza");
  legend->AddEntry(graphInduttanza, "Tensione induttanza");
  legend->AddEntry(graphCondensatore, "Tensione condensatore");
  legend->AddEntry(graphTotale, "Tensione generatore");
  legend->AddEntry(funcResistenza, "Fit", "l");
  legend->Draw();

  TPaveText *titolo = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titolo);
  titolo->AddText("Risposta in frequenza - tutte le componenti");
  titolo->Draw();

  canvas->Update();
}