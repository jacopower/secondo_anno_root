// Si scriva la parte rilevante e autoconsistente di codice di una macro di ROOT in cui:
// Si apre un file ROOT (out.root) in modalità di scrittura
// Si scrive sul file ROOT un istogramma monodimensionale di Classe TH1D, il cui puntatore e nome dell’oggetto sono rispettivamente h (puntatore) e histo (nome),
// Si chiude il file root.
// Si apre il file root in modalità di lettura e si recupera dal file, utilizzando il metodo Get, l’istogramma precedentemente scritto sul file

#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
}

void myMacro()
{
  TFile *file = new TFile("out.root", "RECREATE");
  TH1D *h = new TH1D("histo", "titolo", 100, 0, 5);
  h->Write();
  file->Close();

  TFile *fileRead = new TFile("out.root");
  TH1D *histoRecup = (TH1D *)fileRead->Get("histo");
  fileRead->Close();
}
