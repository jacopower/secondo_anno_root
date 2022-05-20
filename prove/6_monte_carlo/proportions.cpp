#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"

void proportions(Int_t nGen)
{
  TH1F *h = new TH1F("h", "abudancies", 5, 0, 5);
  h->SetFillColor(kBlue);

  Double_t x = 0;

  for (Int_t i = 0; i < nGen; ++i)
  {
    x = gRandom->Rndm(); // uniform in [0, 1]

    if (x < 0.1) // 10%
    {
      h->Fill(0);
    }
    else if (x < 0.3) // 20%
    {
      h->Fill(1);
    }
    else if (x < 0.7) // 40%
    {
      h->Fill(2);
    }
    else if (x < 0.95) // 25%
    {
      h->Fill(3);
    }
    else // 5%
    {
      h->Fill(4);
    }
  }

  TCanvas *canvas = new TCanvas("c1", "Generate proportions Example", 200, 10, 600, 400);
  h->Draw("H");
  h->Draw("E, SAME");
}