#include "TH1F.h"
#include "TRandom.h"

void func()
{
  TH1F *h = new TH1F("h", "prova", 100, 0, 10);

  Double_t x{};
  for (int i = 0; i < 1000; ++i)
  {
    x = gRandom->Uniform(0, 10);
    h->Fill(x);
  }

  h->GetXaxis()->SetTitle("Asse x");
  h->GetYaxis()->SetTitle("Asse y");
  h->GetXaxis()->ChangeLabel(1, 45, -1, -1, -1, -1, "x1");
  h->GetXaxis()->ChangeLabel(2, 45, -1, -1, -1, -1, "x2");
  h->GetXaxis()->ChangeLabel(3, 45, -1, -1, -1, -1, "x3");
  h->GetXaxis()->ChangeLabel(4, 45, -1, -1, -1, -1, "x4");
  h->GetXaxis()->ChangeLabel(5, 45, -1, -1, -1, -1, "x5");
  h->GetXaxis()->ChangeLabel(6, 45, -1, -1, -1, -1, "x6");
  h->GetXaxis()->ChangeLabel(7, 45, -1, -1, -1, -1, "x7");
  h->GetXaxis()->ChangeLabel(8, 45, -1, -1, -1, -1, "x8");
  h->GetXaxis()->ChangeLabel(9, 45, -1, -1, -1, -1, "SUCA");
  h->GetXaxis()->ChangeLabel(10, 45, -1, -1, -1, -1, "SUCA");
  h->GetXaxis()->ChangeLabel(11, 45, -1, -1, -1, -1, "SUCA");


  h->Draw();
}