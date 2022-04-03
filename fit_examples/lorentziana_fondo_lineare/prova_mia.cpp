#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

double background(double *x, double *par)
{
  Double_t xx = x[0];
  Double_t result = par[0] + par[1] * xx + par[2] * xx * xx;
  return result;
}

double lorentzianPeak(double *x, double *par)
{
  Double_t xx = x[0];
  Double_t result = (0.5 * par[0] * par[1] / TMath::Pi()) /
                    TMath::Max(1.e-10, (xx - par[2]) * (xx - par[2]) + .25 * par[1] * par[1]);
  return result;
}

double fitFunction(double *x, double *par)
{
  return background(x, par) + lorentzianPeak(x, &par[3]);
}

void fittingDemo()
{
  const int nBins = 60;

  double data[nBins] = {6, 1, 10, 12, 6, 13, 23, 22, 15, 21,
                        23, 26, 36, 25, 27, 35, 40, 44, 66, 81,
                        75, 57, 48, 45, 46, 41, 35, 36, 53, 32,
                        40, 37, 38, 31, 36, 44, 42, 37, 32, 32,
                        43, 44, 35, 33, 33, 39, 29, 41, 32, 44,
                        26, 39, 29, 35, 32, 21, 21, 15, 25, 15};
  TCanvas *c1 = new TCanvas("c1", "Fitting Demo", 10, 10, 700, 500);
  c1->SetFillColor(33);
  c1->SetFrameFillColor(41);
  c1->SetGrid();

  TH1F *histo = new TH1F("histo",
                         "Lorentzian Peak on Quadratic Background", 60, 0, 3);
  histo->SetMarkerStyle(21);
  histo->SetMarkerSize(0.8);
  histo->SetStats(0);

  for (int i = 0; i < nBins; i++)
  {
    histo->SetBinContent(i + 1, data[i]);
  }

  TF1 *fitFunc = new TF1("fitFunc", fitFunction, 0, 3, 6);
  fitFunc->SetNpx(500);
  fitFunc->SetLineWidth(4);
  fitFunc->SetLineColor(kMagenta);

  fitFunc->SetParameters(1, 1, 1, 1, 1, 1);
  histo->Fit("fitFunc");

  //histo->Draw();
  // fitFunc->Draw("Same");
}