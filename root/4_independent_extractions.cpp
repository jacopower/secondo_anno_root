#include "TH2F.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TMath.h"

void testCorrelation(Int_t nGen)
{
  gStyle->SetOptStat(0);

  TH2F *h = new TH2F("h", "", 1000, 0, TMath::Pi(), 1000, 0, 2 * TMath::Pi());
  TH3F *hSpace = new TH3F("hSpace", "", 100, -1, 1, 100, -1, 1, 100, -1, 1);

  h->SetLineColor(1);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitle("#theta (rad)");
  h->GetYaxis()->SetTitle("#phi (rad)");

  Double_t phi, theta = 0;

  for (Int_t i = 0; i < nGen; ++i)
  {
    theta = gRandom->Rndm() * TMath::Pi();
    phi = gRandom->Rndm() * 2 * TMath::Pi();
    h->Fill(theta, phi);
    hSpace->Fill(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  }
}