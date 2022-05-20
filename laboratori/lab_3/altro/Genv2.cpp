#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
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
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  // open a ROOT file
  TFile *file = new TFile("mySimulation.root");

  if (file->IsOpen())
  {
    std::cout << "File opened successfully" << '\n';
  }

  // retrieving histograms
  TH1F *typesGenerated = (TH1F *)file->Get("particlesHisto");
  TH1F *hPhi = (TH1F *)file->Get("phiHisto");
  TH1F *hTheta = (TH1F *)file->Get("thetaHisto");
  TH1F *hImpulse = (TH1F *)file->Get("impulseHisto");
  // TH1F *transverseImpulseHisto = (TH1F *)file->Get("transverseImpulseHisto");
  // TH1F *energyHisto = (TH1F *)file->Get("energyHisto");

  TH1F *invMass_discord = (TH1F *)file->Get("mass1Histo");
  TH1F *invMass_concord = (TH1F *)file->Get("mass2Histo");
  TH1F *invMass_PK_discord = (TH1F *)file->Get("mass3Histo");
  TH1F *invMass_PK_concord = (TH1F *)file->Get("mass4Histo");
  TH1F *invMass_decayment = (TH1F *)file->Get("mass5Histo");
  // TH1F *mass6Histo = (TH1F *)file->Get("mass6Histo");

  // ANALYSIS: SUBTRACTION
  TH1F *subtraction_invMass_PK = new TH1F("subtraction_invMass_PK", "subtraction_invMass_PK", 80, 0, 2);
  TH1F *subtraction_invMass = new TH1F("subtraction_invMass", "subtraction_invMass", 80, 0, 2);
  subtraction_invMass_PK->Sumw2();
  subtraction_invMass_PK->Add(invMass_PK_discord, invMass_PK_concord, 1, -1);
  subtraction_invMass->Sumw2();
  subtraction_invMass->Add(invMass_discord, invMass_concord, 1, -1);

  // fitting
  hTheta->Fit("pol0", "Q");
  hPhi->Fit("pol0", "Q");
  hImpulse->Fit("expo", "Q");
  subtraction_invMass_PK->Fit("gaus", "Q", "", 0.6, 1.2);
  subtraction_invMass->Fit("gaus", "Q", "", 0.6, 1.2);
  invMass_decayment->Fit("gaus", "Q", "", 0.6, 1.2);

  // adjusting entries
  subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral());
  subtraction_invMass->SetEntries(subtraction_invMass->Integral());

  // saving in a new root file
  TFile *file2 = new TFile("testrootfit.root", "RECREATE");
  hTheta->Write();
  hPhi->Write();
  hImpulse->Write();
  subtraction_invMass_PK->Write();
  subtraction_invMass->Write();
  invMass_decayment->Write();
  file2->Close();

  // Create new Canvas with inv mass results
  TCanvas *first = new TCanvas("c1", "Generation results", 200, 10, 1400, 900);
  first->Divide(2, 2);

  // working on pad 1 ------ types particle generated
  first->cd(1);
  gPad->SetGrid(1);
  typesGenerated->GetXaxis()->SetTitle("Types of Particles");
  typesGenerated->GetYaxis()->SetTitle("NUmber of particles");
  typesGenerated->GetYaxis()->SetTitleOffset(1.3);
  typesGenerated->SetLineColor(65);
  typesGenerated->SetLineWidth(2);
  typesGenerated->SetFillColorAlpha(68, 0.7);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  typesGenerated->DrawCopy();

  // working on pad 2 --------- impulse distribution
  first->cd(1);
  gPad->SetGrid(1);
  hImpulse->GetXaxis()->SetTitle("Impulse (GeV)");
  hImpulse->GetYaxis()->SetTitle("Occurrences");
  hImpulse->GetYaxis()->SetTitleOffset(1.3);
  hImpulse->SetLineColor(65);
  hImpulse->SetLineWidth(2);
  hImpulse->SetFillColorAlpha(68, 0.7);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  hImpulse->DrawCopy();

  // working on pad 3 -------- Theta distribution
  first->cd(3);
  gPad->SetGrid(1);
  hTheta->GetXaxis()->SetTitle("Angle Values (rad)");
  hTheta->GetYaxis()->SetTitle("Occurrences");
  hTheta->GetYaxis()->SetTitleOffset(1.3);
  hTheta->SetLineColor(65);
  hTheta->SetLineWidth(2);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  hTheta->DrawCopy();

  // working on pad 4 ----------------- Phi distribution
  first->cd(4);
  gPad->SetGrid(1);
  hPhi->GetXaxis()->SetTitle("Angle Values (rad)");
  hPhi->GetYaxis()->SetTitle("Occurrences");
  hPhi->GetYaxis()->SetTitleOffset(1.3);
  hPhi->SetLineColor(65);
  hPhi->SetLineWidth(2);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  hPhi->DrawCopy();

  // CREATING SECOND CANVAS WITH INVARIANT MASSES RESULTS
  TCanvas *second = new TCanvas("second", "Comparison of invariant masses", 200, 10, 1400, 900);
  second->Divide(2, 2);

  // Working on pad 5 ------------------ invMassSubtraction 3 - 4
  second->cd(1);
  gPad->SetGrid(1);
  subtraction_invMass_PK->GetXaxis()->SetTitle("Invariant Masses (GeV/c^2)");
  subtraction_invMass_PK->GetYaxis()->SetTitle("Number of particles");
  subtraction_invMass_PK->GetYaxis()->SetTitleOffset(1.3);
  subtraction_invMass_PK->SetLineColor(65);
  subtraction_invMass_PK->SetLineWidth(2);
  subtraction_invMass_PK->SetFillColorAlpha(68, 0.7);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  subtraction_invMass_PK->DrawCopy();

  // Working on pad 6 ------------ invMassSubtraction 1 - 2
  second->cd(3);
  gPad->SetGrid(1);
  subtraction_invMass->GetXaxis()->SetTitle("Invariant Masses (GeV/c^2)");
  subtraction_invMass->GetYaxis()->SetTitle("Number of particles");
  subtraction_invMass->GetYaxis()->SetTitleOffset(1.3);
  subtraction_invMass->SetLineColor(65);
  subtraction_invMass->SetLineWidth(2);
  subtraction_invMass->SetFillColorAlpha(68, 0.7);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  subtraction_invMass->DrawCopy();

  // Working on pad 7 ---------- inMass Decay
  second->cd(2);
  gPad->SetGrid(1);
  invMass_decayment->GetXaxis()->SetTitle("Invariant Masses (GeV/c^2)");
  invMass_decayment->GetYaxis()->SetTitle("Number of particles");
  invMass_decayment->GetYaxis()->SetTitleOffset(1.3);
  invMass_decayment->SetLineColor(65);
  invMass_decayment->SetLineWidth(2);
  invMass_decayment->SetFillColorAlpha(68, 0.7);
  gPad->SetFrameFillColor(19);
  gPad->SetFrameBorderSize(5);
  gPad->SetFrameBorderMode(-1);
  invMass_decayment->DrawCopy();

  first->Print("Generation Results.pdf");
  second->Print("Comparison of Invariant Masses.pdf");
  file->Close();
}