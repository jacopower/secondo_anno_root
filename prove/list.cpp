#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TCanvas.h"

void setStyle()
{
  // gStyle->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

Double_t myFunction1(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

Double_t myFunction2(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t yy = x[1];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) * TMath::Exp(-(yy - par[1]) * (yy - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void myMacro(Int_t nGen = 10000)
{
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  TH1 *h[3];
  char *histName = new char[10];

  for (Int_t i = 0; i < 3; ++i)
  {
    sprintf(histName, "h%d", i);

    if (i < 2)
    {
      h[i] = new TH1F(histName, "Test Histogram", 100, -5., 5.);
    }
    else
    {
      h[2] = new TH2F(histName, "Test Histogram", 100, -5., 5., 100, -5., 5.);
    }

    // cosmetics
    h[i]->SetMarkerStyle(20);
    h[i]->SetMarkerSize(0.5);
    h[i]->SetLineColor(1);
    h[i]->GetYaxis()->SetTitleOffset(1.2);
    h[i]->GetXaxis()->SetTitleSize(0.04);
    h[i]->GetYaxis()->SetTitleSize(0.04);
    h[i]->GetXaxis()->SetTitle("x");
    h[i]->GetYaxis()->SetTitle("Entries");
  }

  h[0]->SetFillColor(4);
  h[1]->SetFillColor(2);

  // Filling the histo with FIllRandom
  h[0]->FillRandom("gaus", nGen);

  TF1 *f1 = new TF1("f1", myFunction1, -10, 10, 3);
  f1->SetParameter(0, 1);
  f1->SetParameter(1, 0);
  f1->SetParameter(2, 1);
  h[1]->FillRandom("f1", nGen);

  TF2 *f2 = new TF2("f2", myFunction2, -10, 10, -10, 10, 3);
  f2->SetParameter(0, 1);
  f2->SetParameter(1, 0);
  f2->SetParameter(2, 1);
  h[2]->FillRandom("f2", nGen);

  // Graph
  int const n_points = 10;
  double x_vals[n_points] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double y_vals[n_points] = {6, 12, 14, 20, 22, 24, 35, 45, 44, 53};
  double x_errs[n_points] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double y_errs[n_points] = {5, 5, 4.7, 4.5, 4.2, 5.1, 2.9, 4.1, 4.8, 5.43};

  TGraphErrors *graph = new TGraphErrors(n_points, x_vals, y_vals, 0, y_errs);
  graph->SetTitle("Verifica della legge di Ohm;DDP (V); Corrente (A)");

  // TList
  TList *list = new TList();

  for (Int_t i = 0; i < 3; ++i)
  {
    list->Add(h[i]);
  }
  list->Add(graph);

  list->Print();

  // Drawing
  TCanvas *canvas = new TCanvas("canvas", "TF1 List", 200, 10, 600, 400);
  canvas->Divide(2, 2);

  for (Int_t i = 0; i < 4; ++i)
  {
    canvas->cd(i + 1);

    if (!list->At(i)->InheritsFrom("TGraph")) // RTTI feature of ROOT
    {
      list->At(i)->Draw("H");
    }
    else
    {
      list->At(i)->Draw("APE");
    }
  }
}