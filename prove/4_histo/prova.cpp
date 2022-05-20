#include "TF1.h"

void prova1()
{
  TF1 *f1 = new TF1("f1", "sin(x) / x", 0, 10);
  f1->Draw();

  TF1 *f2 = new TF1("f2", "f1 * 2", 0, 10);
  f2->Draw();

  TF1 *f3 = new TF1("f3", "[0] * sin([1] * x)", -3, 3);
  f3->SetParameter(0, 10);
  f3->SetParameter(1, 5);
  f3->Draw();
}

Double_t MyFunction(Double_t *x, Double_t *par)
{
  Float_t xx = x[0];
  Double_t val = TMath::Abs(par[0] * sin(par[1] * xx) / xx);
  return val;
}

void prova2()
{
  TF1 *f1 = new TF1("f1", MyFunction, 0, 10, 2); // 2 = numero parametri
  f1->SetParameters(2, 1);                       // inizializzo i parametri a 2 e 1
}