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

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(1111);
}

void myMacro()
{
  auto mg = new TMultiGraph();

  auto gr1 = new TGraph();
  gr1->SetMarkerStyle(20);
  auto gr2 = new TGraph();
  gr2->SetMarkerStyle(21);
  auto gr3 = new TGraph();
  gr3->SetMarkerStyle(23);
  auto gr4 = new TGraph();
  gr4->SetMarkerStyle(24);

  Double_t dx = 6.28 / 100;
  Double_t x = -3.14;

  for (int i = 0; i <= 100; i++)
  {
    x = x + dx;
    gr1->SetPoint(i, x, 2. * TMath::Sin(x));
    gr2->SetPoint(i, x, TMath::Cos(x));
    gr3->SetPoint(i, x, TMath::Cos(x * x));
    gr4->SetPoint(i, x, TMath::Cos(x * x * x));
  }

  mg->Add(gr4, "PL");
  mg->Add(gr3, "PL");
  mg->Add(gr2, "*L");
  mg->Add(gr1, "PL");

  TCanvas *c = new TCanvas();
  c->cd();
  mg->Draw("A pmc plc");
}

void prova()
{
  TMultiGraph *mg = new TMultiGraph("mg", "mg");
  const Int_t n = 10;
   Double_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->SetTitle("TGraphErrors Example");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);

   Double_t x1[100], y1[100];
   Int_t n1 = 20;
   for (Int_t i=0;i<n1;i++) {
     x1[i] = i*0.1;
     y1[i] = 10*sin(x1[i]+0.2);
   }
   TGraph* gr1 = new TGraph(n1,x1,y1);

  mg->Add(gr);
  mg->Add(gr1);

   TCanvas *c1 = new TCanvas();
   c1->Divide(2, 1);
   c1->cd(1);
  gr->Draw("ALP");
   c1->cd(2);
   gr1->Draw("AC*");

   TCanvas *c2 = new TCanvas();
   c2->cd();
   mg->Draw("ALP"); // HO AGGIUNTO STA MERDA
   c2->BuildLegend();



}