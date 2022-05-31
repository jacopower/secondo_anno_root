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

constexpr Double_t width = 1280;
constexpr Double_t height = 480;

void setStyle()
{
  gROOT->SetStyle("BELLE2");
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
}

void setGraphicsGraph(TGraphErrors *graph)
{
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kAzure - 1);
  graph->SetLineColor(kAzure - 1);

  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->SetTitleOffset(1.4);
  asseX->SetTitleSize(0.03);
  asseX->SetTitleFont(1);

  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

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
  asseY->SetTitleFont(1);

  asseY->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseY->SetNoExponent(); // no exp su assi

  asseY->SetLabelFont(1);
  asseY->SetLabelSize(0.025);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(1010); // 10 divsioni secondarie, 30 divisioni primarie
  // asseY->SetTickSize(0.03);
  // asseY->SetTickLength(0.03);
}

void setGraphicsFit(TF1 *func)
{
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
  box->SetBorderSize(1);
  box->SetTextAlign(12);
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

void sweep_ampiezza()
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

  // ***** PLOT *****s
  // Canvas
  TCanvas *canvas = new TCanvas("canvas", "Resistenza / Generatore -  Sweep ampiezza", 0, 0, width, height);
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
  // graphResistenza->SetMaximum(maxPlotResistenza);
  graphResistenza->Draw("ALP");

  // pad Generatore
  padGeneratore->cd();
  padGeneratore->SetGridx();
  padGeneratore->SetGridy();
  padGeneratore->SetFrameLineWidth(2);
  // graphTotale->SetMaximum(maxPlotTotale);
  graphTotale->Draw("ALP");

  canvas->Update();
}