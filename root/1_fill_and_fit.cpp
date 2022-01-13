#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TFile.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
}

// user defined-Gauss
Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void myMacro(Int_t nGen = 1E5)
{
  // TFile *file = new TFile("example.root", "RECREATE");
  // Draw() -----> DrawCopy()

  char *histName = new char[10];
  TH1F *h[2];

  for (Int_t i = 0; i < 2; ++i)
  {
    sprintf(histName, "h%d", i);
    h[i] = new TH1F(histName, "test histo", 100, -5., 5.);

    h[i]->SetLineColor(kBlack);
    h[i]->SetMarkerStyle(kFullCircle);
    h[i]->SetMarkerSize(0.5);
    h[i]->GetYaxis()->SetTitleOffset(1.2);
    h[i]->GetXaxis()->SetTitleSize(0.04);
    h[i]->GetYaxis()->SetTitleSize(0.04);
    h[i]->GetXaxis()->SetTitle("x");
    h[i]->GetYaxis()->SetTitle("Entries");
  }

  h[0]->SetFillColor(kBlue);
  h[1]->SetFillColor(kRed);

  h[0]->FillRandom("gaus", nGen); // predef gauss (0,1)

  TF1 *f = new TF1("f", myFunction, -10, 10, 3); // 3 = #parameters
  f->SetParameter(0, 1);
  f->SetParameter(1, 0);
  f->SetParameter(2, 1);
  h[1]->FillRandom("f", nGen); // user-def gauss

  TCanvas *c = new TCanvas("c", "fill random canvas", 200, 10, 600, 400);
  c->Divide(1, 2);

  for (Int_t i = 0; i < 2; ++i)
  {
    c->cd(i + 1);
    h[i]->Draw("H");
    h[i]->Draw("E, P, SAME");
  }

  Double_t hEntries = h[0]->GetEntries();
  Double_t hMean = h[0]->GetMean();
  Double_t hRMS = h[0]->GetRMS();
  Double_t hMeanErr = h[0]->GetMeanError();
  Double_t hRMSErr = h[0]->GetRMSError();
  std::cout << "Location of Maximum: " << h[0]->GetBinCenter(h[0]->GetMaximumBin()) << '\n';
  std::cout << "Maximum of the distribution: " << h[0]->GetMaximum() << '\n';

  TCanvas *c1 = new TCanvas("c1", "fit1 canvas", 200, 10, 600, 400);
  h[0]->Fit("gaus");
  h[0]->Draw("H");
  h[0]->Draw("E, P, SAME");

  TCanvas *c2 = new TCanvas("c2", "fit2 canvas", 200, 10, 600, 400);
  f->SetParameter(0, 4000);
  f->SetParameter(1, 0);
  f->SetParameter(2, 1);
  f->SetLineColor(kCyan);

  TFitResultPtr fRes = h[1]->Fit(f, "S");

  h[1]->Draw("H");
  h[1]->Draw("E, P, SAME");

  TLegend *leg = new TLegend(.1, .7, .3, .9, "test fit");
  leg->SetFillColor(0);
  leg->AddEntry(h[1], "Punti sperimentali");
  leg->AddEntry(f, "Fit Gaussiano");
  leg->Draw("S");

  h[1]->GetListOfFunctions()->Print();
  TF1 *fitFunc = h[1]->GetFunction("f");

  Double_t fMean = fitFunc->GetParameter(1);
  Double_t fMeanErr = fitFunc->GetParError(1);
  Double_t fSigma = fitFunc->GetParameter(2);
  Double_t fSigmaErr = fitFunc->GetParError(2);
  Double_t chiSquare = fitFunc->GetChisquare();
  Double_t nDOF = fitFunc->GetNDF();
  Double_t Prob = fitFunc->GetProb();

  TMatrixD cov = fRes->GetCovarianceMatrix();
  TMatrixD cor = fRes->GetCorrelationMatrix();
  cov.Print();
  cor.Print();

  // file->Write();
  // file->Close();
}