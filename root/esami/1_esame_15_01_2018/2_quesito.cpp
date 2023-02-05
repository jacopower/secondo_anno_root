// Scrivi Macro in cui
// 1) Si definiscono 2 istogrammi monodimensionali di 500 bin in un range da 0 a 5
// 2) Si riempe il primo istogramma con 10^6 occorrenze generate esplicitamente e singolarmente secondo gaussiana
// con media 2.5 e deviazione standard 0.25.
// 3) Si riempe il secondo istogramma con 10^4 occorrenze generate esplicitamente e singolarmente secondo una distribuzione uniforme del range [0,5]
// 4) Si fa la somma dei due istogrammi, e si effettua il Fit dell’istogramma somma secondo una forma funzionale consistente di una
// gaussiana (3 parametri: ampiezza,media e deviazione standard) e un polinomio di grado 0 (1 parametro), per un totale di 4 parametri liberi.
// 5) Si stampa a schermo il valore dei parametri dopo il fit, con relativo errore, e il ChiQuadro ridotto

#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"


void setStyle()
{
  gROOT->SetStyle("Plain");   // the Plain style gives"conventional" Postscript output (the only other style is Default)
  gStyle->SetPalette(57);     // sets the palette called kBird (kCool is num 109)
  gStyle->SetOptStat(2210);   // prints on screen (graph legend): kurtosis; skewness; integral; overflows; underflows; rms(stdev); mean; entries; name
  gStyle->SetOptFit(1111);    // prints on screen (graph legend): probability; chisquare; error; value prameters
}


Double_t myFunction(Double_t *x, Double_t *par)     // gauss distrib function + polinomio ( = par[3]) // x and par are used as arrays for variables and parameters
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]) + par[3];
  return val;
}     // xx=x[0] =variable; par[0] =ampiezza ; par[1] = media ; par[2] = dev st


void myMacro()
{
  constexpr Int_t nGen1 = 1E6;
  constexpr Int_t nGen2 = 1E4;

  TH1F *h1 = new TH1F("h1", "histo1", 500, 0., 5.);     // creating histo 1
  TH1F *h2 = new TH1F("h2", "histo2", 500, 0., 5.);     // creating histo 2

  Double_t x1, x2 = 0;
  for (Int_t i = 0; i < nGen1; ++i)     // generating x and filling h1 1E6 times
  {
    x1 = gRandom->Gaus(2.5, 0.25);      // generating x from a gauss distribution of mean=2.5 and stdev=0.25
    h1->Fill(x1);                       // filling h1 with x
  }

  for (Int_t i = 0; i < nGen2; ++i)     // filling h2 ...
  {
    x2 = gRandom->Uniform(0., 5.);      // generating x from a random uniform distribution from 0 to 5
    h2->Fill(x2);
  }

  TH1F *hSum = new TH1F(*h1);           // creating new histo
  hSum->SetTitle("Somma Istogrammi");
  hSum->SetName("hSum");
  hSum->Add(h1, h2, 1, 1);              // hSum = h1 + h2

  TF1 *f = new TF1("f", myFunction, 0., 5., 4);     // creating functional form with myFunction in [0,5]
  f->SetParameters(1, 1, 1, 1);                     // initializes the parameters of the function
  f->SetLineColor(kRed);                            // sets line color (1 to 5 or color code)
  f->SetLineWidth(1);                               // sets line with 

  TCanvas *canvas = new TCanvas("canvas");
  hSum->Fit("f");                     // fittin f on hSum
  hSum->Draw();                       // draw hSum (fitted) on canvas

  TCanvas *c = new TCanvas("c");      // (not requested...but ok)
  c->Divide(1, 2);                    // divides c in a 1x2 grid
  c->cd(1);                           // calls grid space n1
  h1->Draw();                         // draws h1 (in cd(1))
  c->cd(2);
  h2->Draw();

  TF1 *fitFunc = hSum->GetFunction("f");            // gets the fitted funcion extracted from the fitted histo hSum
  Double_t ampiezza = fitFunc->GetParameter(0);
  Double_t media = fitFunc->GetParameter(1);
  Double_t stdDev = fitFunc->GetParameter(2);
  Double_t A = fitFunc->GetParameter(3);
  Double_t ampiezzaErr = fitFunc->GetParError(0);
  Double_t mediaErr = fitFunc->GetParError(1);
  Double_t stdDevErr = fitFunc->GetParError(2);
  Double_t AErr = fitFunc->GetParError(3);
  Double_t chiSquare = fitFunc->GetChisquare();
  Double_t NDOF = fitFunc->GetNDF();

  std::cout << "***** PARAMETRI FIT *****" << '\n';
  std::cout << "Ampiezza: " << ampiezza << " +/- " << ampiezzaErr << '\n';
  std::cout << "Media: " << media << " +/- " << mediaErr << '\n';
  std::cout << "StdDEV: " << stdDev << " +/- " << stdDevErr << '\n';
  std::cout << "A: " << A << " +/- " << AErr << '\n';
  std::cout << "Chi/NDOF: " << chiSquare / NDOF << '\n';
  std::cout << "**********" << '\n';
}
