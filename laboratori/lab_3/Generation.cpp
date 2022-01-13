#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include <iostream>

void compileGeneration()
{
  gROOT->LoadMacro("ParticleType.cpp+");
  gROOT->LoadMacro("ResonanceType.cpp+");
  gROOT->LoadMacro("Particle.cpp+");
  gROOT->LoadMacro("main.cpp+");
}

void analyseGeneration()
{
  // ***** ROOT STYLE *****
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  // ***** READ HISTOGRAMS FROM FILE *****
  TFile *file = new TFile("mySimulation.root");

  if (file->IsOpen())
  {
    std::cout << "File opened successfully" << '\n';
  }

  TH1F *particlesHisto = (TH1F *)file->Get("particlesHisto");
  TH1F *phiHisto = (TH1F *)file->Get("phiHisto");
  TH1F *thetaHisto = (TH1F *)file->Get("thetaHisto");
  TH1F *impulseHisto = (TH1F *)file->Get("impulseHisto");
  TH1F *mass1Histo = (TH1F *)file->Get("mass1Histo");
  TH1F *mass2Histo = (TH1F *)file->Get("mass2Histo");
  TH1F *mass3Histo = (TH1F *)file->Get("mass3Histo");
  TH1F *mass4Histo = (TH1F *)file->Get("mass4Histo");
  TH1F *mass5Histo = (TH1F *)file->Get("mass5Histo");

  // ***** SUBTRACTIONS *****
  TH1F *subtractPK = new TH1F("subtractPK", "subtraction_invMass_PK", 80, 0, 2);
  TH1F *subtractAll = new TH1F("subtractAll", "subtraction_invMass", 80, 0, 2);
  subtractPK->Sumw2();
  subtractPK->Add(mass3Histo, mass4Histo, 1, -1);
  subtractAll->Sumw2();
  subtractAll->Add(mass1Histo, mass2Histo, 1, -1);

  // ***** FITTING *****
  TF1 *uniformFit = new TF1("uniformFit", "pol0");
  TF1 *expoFit = new TF1("expoFit", "expo");
  TF1 *gausFit = new TF1("gausFit", "gaus");

  uniformFit->SetLineColor(kRed);
  uniformFit->SetLineWidth(2);
  expoFit->SetLineColor(kRed);
  expoFit->SetLineWidth(2);
  gausFit->SetLineColor(kRed);
  gausFit->SetLineWidth(2);

  // ***** ADJUST ENTRIES *****
  subtractPK->SetEntries(subtractPK->Integral());
  subtractAll->SetEntries(subtractAll->Integral());

  // ***** SAVE IN A ROOT FILE *****
  TFile *finalHisto = new TFile("finalHisto.root", "RECREATE");
  thetaHisto->Write();
  phiHisto->Write();
  impulseHisto->Write();
  subtractPK->Write();
  subtractAll->Write();
  mass5Histo->Write();
  finalHisto->Close();

  // ***** FIRST CANVAS WITH INV_MASS RESULTS *****
  TCanvas *cResults = new TCanvas("cResults", "Risultati Finali", 200, 10, 1400, 900);
  cResults->Divide(2, 2);

  // _____ pad 1 _____
  cResults->cd(1);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  particlesHisto->GetXaxis()->SetTitle("Particles Generated");
  particlesHisto->GetYaxis()->SetTitle("Number of Particles");
  particlesHisto->GetXaxis()->SetLabelSize(0.05);
  particlesHisto->GetXaxis()->SetBinLabel(1, "Pi+");
  particlesHisto->GetXaxis()->SetBinLabel(2, "Pi-");
  particlesHisto->GetXaxis()->SetBinLabel(3, "K+");
  particlesHisto->GetXaxis()->SetBinLabel(4, "K-");
  particlesHisto->GetXaxis()->SetBinLabel(5, "P+");
  particlesHisto->GetXaxis()->SetBinLabel(6, "P-");
  particlesHisto->GetXaxis()->SetBinLabel(7, "K*");
  particlesHisto->GetXaxis()->SetTitleOffset(1.3);
  particlesHisto->GetYaxis()->SetTitleOffset(1.3);
  particlesHisto->GetXaxis()->SetTitleSize(0.04);
  particlesHisto->GetYaxis()->SetTitleSize(0.04);
  particlesHisto->SetLineColor(kBlue - 3);
  particlesHisto->SetLineWidth(2);
  particlesHisto->SetFillColor(kAzure + 1);

  particlesHisto->DrawCopy();

  // _____ pad 2 _____
  cResults->cd(2);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  impulseHisto->GetXaxis()->SetTitle("Impulse (GeV)");
  impulseHisto->GetYaxis()->SetTitle("Occurrences");
  impulseHisto->GetXaxis()->SetTitleOffset(1.1);
  impulseHisto->GetYaxis()->SetTitleOffset(1.3);
  impulseHisto->GetXaxis()->SetTitleSize(0.04);
  impulseHisto->GetYaxis()->SetTitleSize(0.04);
  impulseHisto->SetLineColor(kBlue - 3);
  impulseHisto->SetLineWidth(2);
  impulseHisto->SetFillColor(kAzure + 1);

  impulseHisto->Fit("expoFit");
  impulseHisto->DrawCopy();

  // _____ pad 3 _____
  cResults->cd(3);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  thetaHisto->GetXaxis()->SetTitle("Angle (rad)");
  thetaHisto->GetYaxis()->SetTitle("Occurrences");
  thetaHisto->GetXaxis()->SetTitleOffset(1.1);
  thetaHisto->GetYaxis()->SetTitleOffset(1.3);
  thetaHisto->GetXaxis()->SetTitleSize(0.04);
  thetaHisto->GetYaxis()->SetTitleSize(0.04);
  thetaHisto->SetLineColor(kBlue - 3);
  thetaHisto->SetLineWidth(2);

  thetaHisto->Fit("uniformFit", "Q");
  thetaHisto->DrawCopy();

  // _____ pad 4 _____
  cResults->cd(4);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  phiHisto->GetXaxis()->SetTitle("Angle (rad)");
  phiHisto->GetYaxis()->SetTitle("Occurrences");
  phiHisto->GetXaxis()->SetTitleOffset(1.1);
  phiHisto->GetYaxis()->SetTitleOffset(1.3);
  phiHisto->GetXaxis()->SetTitleSize(0.04);
  phiHisto->GetYaxis()->SetTitleSize(0.04);
  phiHisto->SetLineColor(kBlue - 3);
  phiHisto->SetLineWidth(2);

  phiHisto->Fit("uniformFit", "Q");
  phiHisto->DrawCopy();

  // ***** SECOND CANVAS *****
  TCanvas *cInvMass = new TCanvas("cInvMass", "Masse Invarianti", 200, 10, 1400, 900);
  cInvMass->Divide(2, 2);

  // _____ pad 1 _____
  cInvMass->cd(1);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  subtractPK->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  subtractPK->GetYaxis()->SetTitle("Number of Particles");
  subtractPK->GetXaxis()->SetTitleOffset(1.1);
  subtractPK->GetYaxis()->SetTitleOffset(1.3);
  subtractPK->GetXaxis()->SetTitleSize(0.04);
  subtractPK->GetYaxis()->SetTitleSize(0.04);
  subtractPK->SetLineColor(kBlue - 3);
  subtractPK->SetLineWidth(2);

  subtractPK->Fit("gausFit", "Q", "", 0.6, 1.2);
  subtractPK->DrawCopy();

  // _____ pad 2 _____
  cInvMass->cd(2);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  mass5Histo->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  mass5Histo->GetYaxis()->SetTitle("Number of Particles");
  mass5Histo->GetXaxis()->SetTitleOffset(1.1);
  mass5Histo->GetYaxis()->SetTitleOffset(1.3);
  mass5Histo->GetXaxis()->SetTitleSize(0.04);
  mass5Histo->GetYaxis()->SetTitleSize(0.04);
  mass5Histo->SetLineColor(kBlue - 3);
  mass5Histo->SetLineWidth(2);
  mass5Histo->SetFillColor(kAzure + 1);

  mass5Histo->Fit("gausFit", "Q", "", 0.6, 1.2);
  mass5Histo->DrawCopy();

  // _____ pad 3 _____
  cInvMass->cd(3);

  gPad->SetGrid(1);
  gPad->SetFrameFillColor(18);
  subtractAll->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  subtractAll->GetYaxis()->SetTitle("Number of Particles");
  subtractAll->GetXaxis()->SetTitleOffset(1.1);
  subtractAll->GetYaxis()->SetTitleOffset(1.3);
  subtractAll->GetXaxis()->SetTitleSize(0.04);
  subtractAll->GetYaxis()->SetTitleSize(0.04);
  subtractAll->SetLineColor(kBlue - 3);
  subtractAll->SetLineWidth(2);

  subtractAll->Fit("gausFit", "Q", "", 0.6, 1.2);
  subtractAll->DrawCopy();

  // _____ pad 4 _____
  cInvMass->cd(4);

  // ***** STANDARD OUTPUT *****

  // _____ particles generated _____
  std::cout << "***** PARTICLES GENERATED *****" << '\n';
  std::cout << "Pi+: " << particlesHisto->GetBinContent(1) << " +/- " << particlesHisto->GetBinError(1) << '\n';
  std::cout << "Pi-: " << particlesHisto->GetBinContent(2) << " +/- " << particlesHisto->GetBinError(2) << '\n';
  std::cout << "K+: " << particlesHisto->GetBinContent(3) << " +/- " << particlesHisto->GetBinError(3) << '\n';
  std::cout << "K-: " << particlesHisto->GetBinContent(4) << " +/- " << particlesHisto->GetBinError(4) << '\n';
  std::cout << "P+: " << particlesHisto->GetBinContent(5) << " +/- " << particlesHisto->GetBinError(5) << '\n';
  std::cout << "P-: " << particlesHisto->GetBinContent(6) << " +/- " << particlesHisto->GetBinError(6) << '\n';
  std::cout << "K*: " << particlesHisto->GetBinContent(7) << " +/- " << particlesHisto->GetBinError(7) << '\n';
  std::cout << "**********" << '\n';

  // ***** PRINTING & CLOSING FILE *****
  cResults->Print("Risultati Finali.pdf");
  cInvMass->Print("Masse Invarianti.pdf");
  file->Close();
}