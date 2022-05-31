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

constexpr Double_t width = 1280;
constexpr Double_t height = 1440;

void setStyle()
{
  gROOT->SetStyle("BELLE2");
  // gStyle->SetTextFont(4);
  //   VEDI TUTTE LE OPZIONI DI FONT
  //  gStyle->SetStatFont()
  gStyle->SetPalette(57); // NON CAPISCO CHE FA
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(112211); // 1=Integral 1=Overf 1=Underf 2=RMS 2=Mean 1=Entries 1=Name
  gStyle->SetOptFit(1111);    // 1=Prob 1=Chi 1=Err 1=Param
}

void setGraphicsGraph(TGraphErrors *graph)
{
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kAzure - 1);
  graph->SetLineColor(kAzure - 1);
  graph->SetMarkerSize(0.9);

  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.2);
  asseX->SetTitleSize(0.04);
  asseX->SetTitleFont(2);

  asseX->SetLabelFont(2);
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
  asseY->SetTitleFont(2);

  asseY->SetLabelFont(1);
  asseY->SetLabelSize(0.04);
  asseY->SetLabelOffset(0.01);

  asseY->SetNdivisions(510); // 10 divsioni secondarie, 30 divisioni primarie
  asseY->SetTickLength(0.03);
}

void setMultiPlot(TMultiGraph *graph)
{
  // *****  ASSE X *****
  TAxis *asseX = graph->GetXaxis();
  asseX->CenterTitle();
  asseX->SetMaxDigits(3); // massimo numero cifre, dopo notazione scientifica
  asseX->SetNoExponent(); // no exp su assi

  asseX->SetTitleOffset(1.3);
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

  asseY->SetTitleOffset(1.1);
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
  // func->SetLineStyle(2);
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
  box->SetBorderSize(1);
  box->SetTextAlign(12);
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

void time()
{
  setStyle();

  // ***** GRAFICI SOTTO RISONANZA *****
  TGraphErrors *graphResistenza_sotto = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza_sotto->SetTitle("Tensione resistenza; time (s); Amplitude (V)");
  setGraphicsGraph(graphResistenza_sotto);
  TGraphErrors *graphInduttanza_sotto = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza_sotto->SetTitle("Tensione induttanza; time (s); Amplitude (V)");
  setGraphicsGraph(graphInduttanza_sotto);
  TGraphErrors *graphCondensatore_sotto = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore_sotto->SetTitle("Tensione condensatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphCondensatore_sotto);
  TGraphErrors *graphTotale_sotto = new TGraphErrors("data/ampiezza_tempo/sotto_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale_sotto->SetTitle("Tensione generatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphTotale_sotto);

  // ***** GRAFICI IN RISONANZA *****
  TGraphErrors *graphResistenza_in = new TGraphErrors("data/ampiezza_tempo/in_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza_in->SetTitle("Tensione resistenza; time (s); Amplitude (V)");
  setGraphicsGraph(graphResistenza_in);
  TGraphErrors *graphInduttanza_in = new TGraphErrors("data/ampiezza_tempo/in_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza_in->SetTitle("Tensione induttanza; time (s); Amplitude (V)");
  setGraphicsGraph(graphInduttanza_in);
  TGraphErrors *graphCondensatore_in = new TGraphErrors("data/ampiezza_tempo/in_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore_in->SetTitle("Tensione condensatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphCondensatore_in);
  TGraphErrors *graphTotale_in = new TGraphErrors("data/ampiezza_tempo/in_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale_in->SetTitle("Tensione generatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphTotale_in);

  // ***** GRAFICI SOPRA RISONANZA *****
  TGraphErrors *graphResistenza_sopra = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/resistenza.txt", "%lg %lg %lg");
  graphResistenza_sopra->SetTitle("Tensione resistenza; time (s); Amplitude (V)");
  setGraphicsGraph(graphResistenza_sopra);
  TGraphErrors *graphInduttanza_sopra = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/induttanza.txt", "%lg %lg %lg");
  graphInduttanza_sopra->SetTitle("Tensione induttanza; time (s); Amplitude (V)");
  setGraphicsGraph(graphInduttanza_sopra);
  TGraphErrors *graphCondensatore_sopra = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/condensatore.txt", "%lg %lg %lg");
  graphCondensatore_sopra->SetTitle("Tensione condensatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphCondensatore_sopra);
  TGraphErrors *graphTotale_sopra = new TGraphErrors("data/ampiezza_tempo/sopra_risonanza/totale.txt", "%lg %lg %lg");
  graphTotale_sopra->SetTitle("Tensione generatore; time (s); Amplitude (V)");
  setGraphicsGraph(graphTotale_sopra);

  // *****MULTIPLOT SOTTO RISONANZA *****
  TMultiGraph *multiGraph_sotto = new TMultiGraph("multiGraph_sotto", "Amplitude Time");
  multiGraph_sotto->SetTitle("Tensione tempo - sotto risonanza; Tempo (ms); Tensione (V)");
  multiGraph_sotto->Add(graphResistenza_sotto);
  graphResistenza_sotto->SetLineColor(kPink + 1);
  graphResistenza_sotto->SetMarkerColor(kPink + 1);
  multiGraph_sotto->Add(graphInduttanza_sotto);
  graphInduttanza_sotto->SetLineColor(kOrange + 1);
  graphInduttanza_sotto->SetMarkerColor(kOrange + 1);
  multiGraph_sotto->Add(graphCondensatore_sotto);
  graphCondensatore_sotto->SetLineColor(kSpring - 6);
  graphCondensatore_sotto->SetMarkerColor(kSpring - 6);
  multiGraph_sotto->Add(graphTotale_sotto);
  setMultiPlot(multiGraph_sotto);

  // ***** MULTIPLOT IN RISONANZA *****
  TMultiGraph *multiGraph_in = new TMultiGraph("multiGraph_in", "Amplitude Time");
  multiGraph_in->SetTitle("Tensione tempo - in risonanza; Tempo (ms); Tensione (V)");
  multiGraph_in->Add(graphResistenza_in);
  graphResistenza_in->SetLineColor(kPink + 1);
  graphResistenza_in->SetMarkerColor(kPink + 1);
  multiGraph_in->Add(graphInduttanza_in);
  graphInduttanza_in->SetLineColor(kOrange + 1);
  graphInduttanza_in->SetMarkerColor(kOrange + 1);
  multiGraph_in->Add(graphCondensatore_in);
  graphCondensatore_in->SetLineColor(kSpring - 6);
  graphCondensatore_in->SetMarkerColor(kSpring - 6);
  multiGraph_in->Add(graphTotale_in);
  setMultiPlot(multiGraph_in);

  // ***** MULTIPLOT SOPRA RISONANZA *****
  TMultiGraph *multiGraph_sopra = new TMultiGraph("multiGraph_sopra", "Amplitude Time");
  multiGraph_sopra->SetTitle("Tensione tempo - sopra risonanza; Tempo (ms); Tensione (V)");
  multiGraph_sopra->Add(graphResistenza_sopra);
  graphResistenza_sopra->SetLineColor(kPink + 1);
  graphResistenza_sopra->SetMarkerColor(kPink + 1);
  multiGraph_sopra->Add(graphInduttanza_sopra);
  graphInduttanza_sopra->SetLineColor(kOrange + 1);
  graphInduttanza_sopra->SetMarkerColor(kOrange + 1);
  multiGraph_sopra->Add(graphCondensatore_sopra);
  graphCondensatore_sopra->SetLineColor(kSpring - 6);
  graphCondensatore_sopra->SetMarkerColor(kSpring - 6);
  multiGraph_sopra->Add(graphTotale_sopra);
  setMultiPlot(multiGraph_sopra);

  // ***** GRAFICO RESISTENZA CON FIT *****
  TGraphErrors *graphResistenza_fit = new TGraphErrors("data/ampiezza_tempo/resistenza_sotto_risonanza.txt", "%lg %lg %lg");
  graphResistenza_fit->SetTitle("Tensione resistenza; Tempo (ms); Tensione (V)");
  setGraphicsGraph(graphResistenza_fit);

  // ***** PLOTTO *****
  TCanvas *canvas = new TCanvas("canvasAmpiezzaTempo", "Amplitude Time", 0, 0, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));
  canvas->SetFillColor(kWhite);

  constexpr Double_t R_mis = 150.47;
  constexpr Double_t L_mis = 11.46 * 1E-3;
  constexpr Double_t C_mis = 157.8 * 1E-9;
  constexpr Double_t V0_mis = 5;
  constexpr Double_t f_mis = 2E3;

  TF1 *funcResistenza = new TF1("funcResistenza", amp_time_resistenza, 0, 0.00479, 6);
  funcResistenza->SetParameters(R_mis, L_mis, C_mis, f_mis, V0_mis);
  funcResistenza->SetParNames("R", "L", "C", "f", "V", "Offset");
  funcResistenza->SetNpx(10000);
  funcResistenza->FixParameter(3, 2E3);
  funcResistenza->FixParameter(4, 5.0);
  setGraphicsFit(funcResistenza);
  TFitResultPtr rResistenza = graphResistenza_fit->Fit(funcResistenza, "REMSQ");

  // Creo le Pad
  TPad *padSotto = new TPad("padFit", "Sotto", 0., 0.5, 0.5, 1., kWhite);
  TPad *padIn = new TPad("padOffset", "in", 0.5, 0.5, 1., 1., kWhite);
  TPad *padSopra = new TPad("padSopra", "Sopra", 0, 0., 0.5, 0.5, kWhite);
  TPad *padFit = new TPad("padFit", "Fit", 0.5, 0., 1., 0.5, kWhite);
  padSotto->Draw();
  padIn->Draw();
  padSopra->Draw();
  padFit->Draw();

  // pad SOTTO
  padSotto->cd();
  padSotto->SetGridx();
  padSotto->SetGridy();
  padSotto->SetFrameLineWidth(2);
  multiGraph_sotto->SetMaximum(8.);
  multiGraph_sotto->SetMinimum(-7);
  multiGraph_sotto->Draw("ALP");

  padSotto->BuildLegend();

  TPaveText *titoloSotto = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloSotto);
  titoloSotto->AddText("Tensione tempo - sotto risonanza");
  titoloSotto->Draw();

  // pad IN
  padIn->cd();
  padIn->SetGridx();
  padIn->SetGridy();
  padIn->SetFrameLineWidth(2);
  multiGraph_in->SetMaximum(8);
  multiGraph_in->SetMinimum(-7);
  multiGraph_in->Draw("ALP");

  padIn->BuildLegend();

  TPaveText *titoloIn = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloIn);
  titoloIn->AddText("Tensione tempo - in risonanza");
  titoloIn->Draw();

  // pad SOPRA
  padSopra->cd();
  padSopra->SetGridx();
  padSopra->SetGridy();
  padSopra->SetFrameLineWidth(2);
  multiGraph_sopra->SetMaximum(8);
  multiGraph_sopra->SetMinimum(-6);
  multiGraph_sopra->Draw("ALP");

  padSopra->BuildLegend();

  TPaveText *titoloSopra = new TPaveText(0, 1., .5, .94, "NDC BL");
  setGraphicsTitolo(titoloSopra);
  titoloSopra->AddText("Tensione tempo - sopra risonanza");
  titoloSopra->Draw();

  // pad FIT
  padFit->cd();
  padFit->SetGridx();
  padFit->SetGridy();
  padFit->SetFrameLineWidth(2);
  graphResistenza_fit->GetXaxis()->SetMaxDigits(3);
  graphResistenza_fit->GetXaxis()->SetNoExponent(false);
  graphResistenza_fit->SetMaximum(3);
  graphResistenza_fit->SetMinimum(-2);
  graphResistenza_fit->Draw("ALP");

  canvas->Update();
}
