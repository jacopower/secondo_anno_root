#include "TF1.h"

Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = 0;
  if (xx > 0 && xx < 1)
  {
    val = par[0] * xx * xx;
  }
  else
  {
    val = par[1];
  }
  return val;
}

void myMacro()
{
  TF1 *f = new TF1("f", myFunction, 0, 10, 2);
  f->SetParameters(1, 5);
  TCanvas *c = new TCanvas();
  f->Draw();
}