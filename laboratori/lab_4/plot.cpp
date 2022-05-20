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

// FORMULA f = W / (2*pi)

constexpr Double_t V0_mis = 5.;
constexpr Double_t R_mis = 149.83;
constexpr Double_t L_mis = 10.43 * 1E-3;
constexpr Double_t C_mis = 158.4 * 1E-9;

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

Double_t amp_freq_resistenza(Double_t *x, Double_t *par) // E' LA FREQUENZA NON W
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
  graphResistenza->SetTitle("Sweep Resistenza; Frequency (Hz); Amplitude (V)");
  graphResistenza->SetMarkerStyle(kOpenCircle);
  graphResistenza->SetMarkerColor(kBlue);
  graphResistenza->SetFillColor(0);

  TGraphErrors *graphInduttanza = new TGraphErrors("data/sweep_ampiezza/sweep_freq_induttanza.txt", "%lg %lg %lg");
  graphInduttanza->SetTitle("Sweep Induttanza; Frequency (Hz); Amplitude (V)");
  graphInduttanza->SetMarkerStyle(kOpenCircle);
  graphInduttanza->SetMarkerColor(kBlue);
  graphInduttanza->SetFillColor(0);

  TGraphErrors *graphCondensatore = new TGraphErrors("data/sweep_ampiezza/sweep_freq_condensatore.txt", "%lg %lg %lg");
  graphCondensatore->SetTitle("Sweep Condensatore; Frequency (Hz); Amplitude (V)");
  graphCondensatore->SetMarkerStyle(kOpenCircle);
  graphCondensatore->SetMarkerColor(kBlue);
  graphCondensatore->SetFillColor(0);

  TGraphErrors *graphTotale = new TGraphErrors("data/sweep_ampiezza/sweep_freq_totale.txt", "%lg %lg %lg");
  graphTotale->SetTitle("Sweep Totale; Frequency (Hz); Amplitude (V)");
  graphTotale->SetMarkerStyle(kOpenCircle);
  graphTotale->SetMarkerColor(kBlue);
  graphTotale->SetFillColor(0);

  // ***** CREO LE FUNZIONI DI FIT *****
  TF1 *funcResistenza = new TF1("funcResistenza", amp_freq_resistenza, 0, 2E4, 4);
  TF1 *funcInduttanza = new TF1("funcInduttanza", amp_freq_induttanza, 0, 2E4, 4);
  TF1 *funcCondensatore = new TF1("funcResistenza", amp_freq_condensatore, 0, 2E4, 4);

  funcResistenza->SetParameters(V0_mis, R_mis, L_mis, C_mis);
  funcResistenza->SetParNames("V0", "R", "L", "C");
  funcResistenza->SetLineColor(kRed);
  funcResistenza->SetLineStyle(2);

  funcInduttanza->SetParameters(V0_mis, R_mis, L_mis, C_mis);
  funcInduttanza->SetParNames("V0", "R", "L", "C");
  funcInduttanza->SetLineColor(kRed);
  funcInduttanza->SetLineStyle(2);

  funcCondensatore->SetParameters(V0_mis, R_mis, L_mis, C_mis);
  funcCondensatore->SetParNames("V0", "R", "L", "C");
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
  TMultiGraph *multiGraph = new TMultiGraph("multiGraph", "Amplitude Sweep - Risultati finali");
  multiGraph->Add(graphResistenza);
  multiGraph->Add(graphInduttanza);
  multiGraph->Add(graphCondensatore);
  multiGraph->Add(graphTotale);
  multiGraph->Draw("ALP"); // COSA FA LP?
  multiCanvas->BuildLegend();
  // Vedi cosa fa
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
  constexpr Double_t f_mis = 0.; // INSERISCIIII

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