#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"

#include "iostream"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

R__LOAD_LIBRARY(ParticleType_cpp.so)
R__LOAD_LIBRARY(ResonanceType_cpp.so)
R__LOAD_LIBRARY(Particle_cpp.so)

void test()
{
  // RECUPERO I DATI DEGLI ISTOGRAMMI
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  TFile *file = new TFile("mySimulation.root");

  if (file->IsOpen())
  {
    std::cout << "File opened successfully" << '\n';
  }

  TH1F *particlesHisto = (TH1F *)file->Get("particlesHisto");
  TH1F *phiHisto = (TH1F *)file->Get("phiHisto");
  TH1F *thetaHisto = (TH1F *)file->Get("thetaHisto");
  TH1F *impulseHisto = (TH1F *)file->Get("impulseHisto");
  TH1F *transverseImpulseHisto = (TH1F *)file->Get("transverseImpulseHisto");
  TH1F *energyHisto = (TH1F *)file->Get("energyHisto");

  TH1F *mass1Histo = (TH1F *)file->Get("mass1Histo");
  TH1F *mass2Histo = (TH1F *)file->Get("mass2Histo");
  TH1F *mass3Histo = (TH1F *)file->Get("mass3Histo");
  TH1F *mass4Histo = (TH1F *)file->Get("mass4Histo");
  TH1F *mass5Histo = (TH1F *)file->Get("mass5Histo");
  TH1F *mass6Histo = (TH1F *)file->Get("mass6Histo");

  // TESTI PARTICLE HISTO
  std::cout << '\n'
            << "-----------------------------" << '\n';
  std::cout << "TESTING PARTICLE HISTO" << '\n';
  std::cout << "Underflows " << particlesHisto->GetBinContent(0) << " +/- " << particlesHisto->GetBinError(0) << '\n';
  std::cout << "Pi+: " << particlesHisto->GetBinContent(1) << " +/- " << particlesHisto->GetBinError(1) << '\n';
  std::cout << "Pi-: " << particlesHisto->GetBinContent(2) << " +/- " << particlesHisto->GetBinError(2) << '\n';
  std::cout << "K+: " << particlesHisto->GetBinContent(3) << " +/- " << particlesHisto->GetBinError(3) << '\n';
  std::cout << "K-: " << particlesHisto->GetBinContent(4) << " +/- " << particlesHisto->GetBinError(4) << '\n';
  std::cout << "P+: " << particlesHisto->GetBinContent(5) << " +/- " << particlesHisto->GetBinError(5) << '\n';
  std::cout << "P-: " << particlesHisto->GetBinContent(6) << " +/- " << particlesHisto->GetBinError(6) << '\n';
  std::cout << "K*: " << particlesHisto->GetBinContent(7) << " +/- " << particlesHisto->GetBinError(7) << '\n';
  std::cout << "Overflows " << particlesHisto->GetBinContent(8) << " +/- " << particlesHisto->GetBinError(8) << '\n';
  std::cout << "---------------------------------" << '\n'
            << '\n';

  // TESTING PHI AND THETA (UNIFORM FIT)
  TCanvas *canvasAngle = new TCanvas("canvasAngle");
  canvasAngle->Divide(2, 0);

  canvasAngle->cd(1);
  phiHisto->Fit("pol0");
  phiHisto->Draw();

  canvasAngle->cd(2);
  thetaHisto->Fit("pol0");
  thetaHisto->Draw();

  // CHE COSA DEVO GUARDARE SULLA MEDIA????????????????????????????????????''
  TCanvas *canvasImpulse = new TCanvas("canvasImpulse");
  impulseHisto->Fit("expo");
  impulseHisto->Draw();

}