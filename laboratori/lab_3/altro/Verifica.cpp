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

// DA LEVARE

void verifica()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(222222);
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

  // DISEGNO ISTOGRAMMI E VEDO LE ENTRIES
  std::cout << "*******************************************" << '\n'
            << "ENTRIES:" << '\n';

  // TCanvas *cparticlesHisto = new TCanvas("cparticlesHisto", "cparticlesHisto");
  // particlesHisto->Draw();
  std::cout << "particlesHisto " << particlesHisto->GetEntries() << '\n';

  // TCanvas *cphiHisto = new TCanvas("cphiHisto", "cphiHisto");
  //  phiHisto->Fit(" pol0 ", "Q ");
  // phiHisto->Draw();
  std::cout << "phiHisto " << phiHisto->GetEntries() << '\n';

  // TCanvas *cthetaHisto = new TCanvas("cthetaHisto", "cthetaHisto");
  //  thetaHisto->Fit(" pol0 ", "Q ");
  // thetaHisto->Draw();
  std::cout << "thetaHisto " << thetaHisto->GetEntries() << '\n';

  // TCanvas *cimpulseHisto = new TCanvas("cimpulseHisto", "cimpulseHisto");
  // impulseHisto->Draw();
  std::cout << "impulseHisto " << impulseHisto->GetEntries() << '\n';

  // TCanvas *ctransverseImpulseHisto = new TCanvas("ctransverseImpulseHisto", "ctransverseImpulseHisto");
  // transverseImpulseHisto->Draw();
  std::cout << "transverseImpulseHisto " << transverseImpulseHisto->GetEntries() << '\n';

  // TCanvas *cenergyHisto = new TCanvas("cenergyHisto", "cenergyHisto");
  // energyHisto->Draw();
  std::cout << "energyHisto " << energyHisto->GetEntries() << '\n';

  // TCanvas *cmass1Histo = new TCanvas("cmass1Histo", "cmass1Histo");
  // mass1Histo->Draw();
  std::cout << "mass1Histo " << mass1Histo->GetEntries() << '\t' << "Invariant Mass - Discord Charges" << '\n';

  // TCanvas *cmass2Histo = new TCanvas("cmass2Histo", "cmass2Histo");
  // mass2Histo->Draw();
  std::cout << "mass2Histo " << mass2Histo->GetEntries() << '\t' << "Invariant Mass - Concord Charges" << '\n';

  // TCanvas *cmass3Histo = new TCanvas("cmass3Histo", "cmass3Histo");
  // mass3Histo->Draw();
  std::cout << "mass3Histo " << mass3Histo->GetEntries() << '\t' << "Invariant Mass - P & K with Discord Charges" << '\n';

  // TCanvas *cmass4Histo = new TCanvas("cmass4Histo", "cmass4Histo");
  // mass4Histo->Draw();
  std::cout << "mass4Histo " << mass4Histo->GetEntries() << '\t' << "Invariant Mass - P & K with Concord Charges" << '\n';

  // TCanvas *cmass5Histo = new TCanvas("cmass5Histo", "cmass5Histo");
  // mass5Histo->Draw();
  std::cout << "mass5Histo " << mass5Histo->GetEntries() << '\t' << "Invariant Mass - K* Decays" << '\n';

  // TCanvas *cmass6Histo = new TCanvas("cmass6Histo", "cmass6Histo");
  // mass6Histo->Draw();
  std::cout << "mass6Histo " << mass6Histo->GetEntries() << '\t' << "Invariant Mass - All Particles" << '\n';

  std::cout << "***************************************" << '\n'
            << '\n';

  // CONTROLLO LE PROPORZIONI DI GENERAZIONE DELLE PARTICELLE
  std::cout << "***************************************" << '\n'
            << "PARTICELLE GENERATE" << '\n';

  std::cout << "Underflows: " << particlesHisto->GetBinContent(0) << " +/- " << particlesHisto->GetBinError(0) << '\n';
  std::cout << "Pi+: " << particlesHisto->GetBinContent(1) << " +/- " << particlesHisto->GetBinError(1) << '\n';
  std::cout << "Pi-: " << particlesHisto->GetBinContent(2) << " +/- " << particlesHisto->GetBinError(2) << '\n';
  std::cout << "K+: " << particlesHisto->GetBinContent(3) << " +/- " << particlesHisto->GetBinError(3) << '\n';
  std::cout << "K-: " << particlesHisto->GetBinContent(4) << " +/- " << particlesHisto->GetBinError(4) << '\n';
  std::cout << "P+: " << particlesHisto->GetBinContent(5) << " +/- " << particlesHisto->GetBinError(5) << '\n';
  std::cout << "P-: " << particlesHisto->GetBinContent(6) << " +/- " << particlesHisto->GetBinError(6) << '\n';
  std::cout << "K*: " << particlesHisto->GetBinContent(7) << " +/- " << particlesHisto->GetBinError(7) << '\n';
  std::cout << "Overflows: " << particlesHisto->GetBinContent(8) << " +/- " << particlesHisto->GetBinError(8) << '\n';

  std::cout << '\n'
            << "PERCENTUALI: " << '\n';
  std::cout << "Pi+: " << (particlesHisto->GetBinContent(1) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "Pi-: " << (particlesHisto->GetBinContent(2) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K+: " << (particlesHisto->GetBinContent(3) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K-: " << (particlesHisto->GetBinContent(4) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "P+: " << (particlesHisto->GetBinContent(5) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "P-: " << (particlesHisto->GetBinContent(6) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K*: " << (particlesHisto->GetBinContent(7) / particlesHisto->GetEntries()) * 100 << '\n';

  std::cout << "**********************************************************" << '\n'
            << '\n';

  // FIT DISTRIBUZIONI ANGOLI E MODULO IMPULSO
  TCanvas *fitCanvas = new TCanvas("fitCanvas", "Fit Istogrammi");
  fitCanvas->Divide(2, 2);

  fitCanvas->cd(1);
  phiHisto->Fit("pol0", "Q");
  phiHisto->Draw("E, SAME");

  fitCanvas->cd(2);
  thetaHisto->Fit("pol0", "Q");
  thetaHisto->Draw("E, SAME");

  fitCanvas->cd(3);
  impulseHisto->Fit("expo", "Q");
  impulseHisto->Draw("E, SAME");

  // FACCIO LE SOTTRAZIONI DEGLI ISTOGRAMMI
  TCanvas *diffCanvas = new TCanvas("diffCanvas", "Differenze Histogrammi");
  diffCanvas->Divide(2, 2);

  // 3 - 4
  diffCanvas->cd(1);
  TH1F *diffKP(mass3Histo);
  diffKP->SetTitle("SOTTRAZIONE P & K CONCORDI E DISCORDI");
  diffKP->SetName("DiffKP");
  // std::cout << "ERORRE MEDIA DIFFKP" << diffKP->GetMeanError() << '\n';
  diffKP->Add(mass3Histo, mass4Histo, 1, -1);
  diffKP->Fit("gaus", "Q");
  diffKP->Draw();

  // histo 5
  diffCanvas->cd(2);
  // mass5Histo->Fit("gaus", "Q");
  mass5Histo->Draw();

  // 1 - 2
  diffCanvas->cd(3);
  TH1F *diffAll(mass2Histo);
  diffAll->SetTitle("SOTTRAZIONE TUTTE PARTICELLE CONCORDI E DISCORDI");
  diffAll->SetName("DiffAll");
  diffAll->Add(mass1Histo, mass2Histo, 1, -1);
  diffAll->Draw();

  std::cout << '\n'
            << "***********************************" << '\n'
            << "ANALISI ISTOGRAMMI DIFFERENZA" << '\n';
  std::cout << "Massima x in DiffPK: " << diffKP->GetXaxis()->GetBinCenter(diffKP->GetMaximumBin()) << '\n';
  std::cout << "Massima x in mass5Histo: " << mass5Histo->GetXaxis()->GetBinCenter(mass5Histo->GetMaximumBin()) << '\n';
  std::cout << "Massima x in diffAll: " << diffAll->GetXaxis()->GetBinCenter(diffAll->GetMaximumBin()) << '\n';

  // FACCIO COME IL BRODISMAN
  TCanvas *provaCanvas = new TCanvas("provaCanvas");
  TH1F *subtraction_invMass_PK = new TH1F("subtraction_invMass_PK", "subtraction_invMass_PK", 80, 0, 2);
  TH1F *subtraction_invMass = new TH1F("subtraction_invMass", "subtraction_invMass", 80, 0, 2);
  subtraction_invMass_PK->Sumw2();
  subtraction_invMass_PK->Add(mass3Histo, mass4Histo, 1, -1);
  subtraction_invMass->Sumw2();
  subtraction_invMass->Add(mass1Histo, mass2Histo, 1, -1);

  thetaHisto->Fit("pol0", "Q");
  phiHisto->Fit("pol0", "Q");
  impulseHisto->Fit("expo", "Q");
  subtraction_invMass_PK->Fit("gaus", "Q", "", 0.6, 1.2);
  subtraction_invMass->Fit("gaus", "Q", "", 0.6, 1.2);
  mass5Histo->Fit("gaus", "Q", "", 0.6, 1.2);

  subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral());
  subtraction_invMass->SetEntries(subtraction_invMass->Integral());
}

void verifica()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(222222);
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

  // DISEGNO ISTOGRAMMI E VEDO LE ENTRIES
  std::cout << "*******************************************" << '\n'
            << "ENTRIES:" << '\n';

  // TCanvas *cparticlesHisto = new TCanvas("cparticlesHisto", "cparticlesHisto");
  // particlesHisto->Draw();
  std::cout << "particlesHisto " << particlesHisto->GetEntries() << '\n';

  // TCanvas *cphiHisto = new TCanvas("cphiHisto", "cphiHisto");
  //  phiHisto->Fit(" pol0 ", "Q ");
  // phiHisto->Draw();
  std::cout << "phiHisto " << phiHisto->GetEntries() << '\n';

  // TCanvas *cthetaHisto = new TCanvas("cthetaHisto", "cthetaHisto");
  //  thetaHisto->Fit(" pol0 ", "Q ");
  // thetaHisto->Draw();
  std::cout << "thetaHisto " << thetaHisto->GetEntries() << '\n';

  // TCanvas *cimpulseHisto = new TCanvas("cimpulseHisto", "cimpulseHisto");
  // impulseHisto->Draw();
  std::cout << "impulseHisto " << impulseHisto->GetEntries() << '\n';

  // TCanvas *ctransverseImpulseHisto = new TCanvas("ctransverseImpulseHisto", "ctransverseImpulseHisto");
  // transverseImpulseHisto->Draw();
  std::cout << "transverseImpulseHisto " << transverseImpulseHisto->GetEntries() << '\n';

  // TCanvas *cenergyHisto = new TCanvas("cenergyHisto", "cenergyHisto");
  // energyHisto->Draw();
  std::cout << "energyHisto " << energyHisto->GetEntries() << '\n';

  // TCanvas *cmass1Histo = new TCanvas("cmass1Histo", "cmass1Histo");
  // mass1Histo->Draw();
  std::cout << "mass1Histo " << mass1Histo->GetEntries() << '\t' << "Invariant Mass - Discord Charges" << '\n';

  // TCanvas *cmass2Histo = new TCanvas("cmass2Histo", "cmass2Histo");
  // mass2Histo->Draw();
  std::cout << "mass2Histo " << mass2Histo->GetEntries() << '\t' << "Invariant Mass - Concord Charges" << '\n';

  // TCanvas *cmass3Histo = new TCanvas("cmass3Histo", "cmass3Histo");
  // mass3Histo->Draw();
  std::cout << "mass3Histo " << mass3Histo->GetEntries() << '\t' << "Invariant Mass - P & K with Discord Charges" << '\n';

  // TCanvas *cmass4Histo = new TCanvas("cmass4Histo", "cmass4Histo");
  // mass4Histo->Draw();
  std::cout << "mass4Histo " << mass4Histo->GetEntries() << '\t' << "Invariant Mass - P & K with Concord Charges" << '\n';

  // TCanvas *cmass5Histo = new TCanvas("cmass5Histo", "cmass5Histo");
  // mass5Histo->Draw();
  std::cout << "mass5Histo " << mass5Histo->GetEntries() << '\t' << "Invariant Mass - K* Decays" << '\n';

  // TCanvas *cmass6Histo = new TCanvas("cmass6Histo", "cmass6Histo");
  // mass6Histo->Draw();
  std::cout << "mass6Histo " << mass6Histo->GetEntries() << '\t' << "Invariant Mass - All Particles" << '\n';

  std::cout << "***************************************" << '\n'
            << '\n';

  // CONTROLLO LE PROPORZIONI DI GENERAZIONE DELLE PARTICELLE
  std::cout << "***************************************" << '\n'
            << "PARTICELLE GENERATE" << '\n';

  std::cout << "Underflows: " << particlesHisto->GetBinContent(0) << " +/- " << particlesHisto->GetBinError(0) << '\n';
  std::cout << "Pi+: " << particlesHisto->GetBinContent(1) << " +/- " << particlesHisto->GetBinError(1) << '\n';
  std::cout << "Pi-: " << particlesHisto->GetBinContent(2) << " +/- " << particlesHisto->GetBinError(2) << '\n';
  std::cout << "K+: " << particlesHisto->GetBinContent(3) << " +/- " << particlesHisto->GetBinError(3) << '\n';
  std::cout << "K-: " << particlesHisto->GetBinContent(4) << " +/- " << particlesHisto->GetBinError(4) << '\n';
  std::cout << "P+: " << particlesHisto->GetBinContent(5) << " +/- " << particlesHisto->GetBinError(5) << '\n';
  std::cout << "P-: " << particlesHisto->GetBinContent(6) << " +/- " << particlesHisto->GetBinError(6) << '\n';
  std::cout << "K*: " << particlesHisto->GetBinContent(7) << " +/- " << particlesHisto->GetBinError(7) << '\n';
  std::cout << "Overflows: " << particlesHisto->GetBinContent(8) << " +/- " << particlesHisto->GetBinError(8) << '\n';

  std::cout << '\n'
            << "PERCENTUALI: " << '\n';
  std::cout << "Pi+: " << (particlesHisto->GetBinContent(1) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "Pi-: " << (particlesHisto->GetBinContent(2) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K+: " << (particlesHisto->GetBinContent(3) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K-: " << (particlesHisto->GetBinContent(4) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "P+: " << (particlesHisto->GetBinContent(5) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "P-: " << (particlesHisto->GetBinContent(6) / particlesHisto->GetEntries()) * 100 << '\n';
  std::cout << "K*: " << (particlesHisto->GetBinContent(7) / particlesHisto->GetEntries()) * 100 << '\n';

  std::cout << "**********************************************************" << '\n'
            << '\n';

  // FIT DISTRIBUZIONI ANGOLI E MODULO IMPULSO
  TCanvas *fitCanvas = new TCanvas("fitCanvas", "Fit Istogrammi");
  fitCanvas->Divide(2, 2);

  fitCanvas->cd(1);
  phiHisto->Fit("pol0", "Q");
  phiHisto->Draw("E, SAME");

  fitCanvas->cd(2);
  thetaHisto->Fit("pol0", "Q");
  thetaHisto->Draw("E, SAME");

  fitCanvas->cd(3);
  impulseHisto->Fit("expo", "Q");
  impulseHisto->Draw("E, SAME");

  // FACCIO LE SOTTRAZIONI DEGLI ISTOGRAMMI
  TCanvas *diffCanvas = new TCanvas("diffCanvas", "Differenze Histogrammi");
  diffCanvas->Divide(2, 2);

  // 3 - 4
  diffCanvas->cd(1);
  TH1F *diffKP(mass3Histo);
  diffKP->SetTitle("SOTTRAZIONE P & K CONCORDI E DISCORDI");
  diffKP->SetName("DiffKP");
  // std::cout << "ERORRE MEDIA DIFFKP" << diffKP->GetMeanError() << '\n';
  diffKP->Add(mass3Histo, mass4Histo, 1, -1);
  diffKP->Fit("gaus", "Q");
  diffKP->Draw();

  // histo 5
  diffCanvas->cd(2);
  // mass5Histo->Fit("gaus", "Q");
  mass5Histo->Draw();

  // 1 - 2
  diffCanvas->cd(3);
  TH1F *diffAll(mass2Histo);
  diffAll->SetTitle("SOTTRAZIONE TUTTE PARTICELLE CONCORDI E DISCORDI");
  diffAll->SetName("DiffAll");
  diffAll->Add(mass1Histo, mass2Histo, 1, -1);
  diffAll->Draw();

  std::cout << '\n'
            << "***********************************" << '\n'
            << "ANALISI ISTOGRAMMI DIFFERENZA" << '\n';
  std::cout << "Massima x in DiffPK: " << diffKP->GetXaxis()->GetBinCenter(diffKP->GetMaximumBin()) << '\n';
  std::cout << "Massima x in mass5Histo: " << mass5Histo->GetXaxis()->GetBinCenter(mass5Histo->GetMaximumBin()) << '\n';
  std::cout << "Massima x in diffAll: " << diffAll->GetXaxis()->GetBinCenter(diffAll->GetMaximumBin()) << '\n';

  // FACCIO COME IL BRODISMAN
  TCanvas *provaCanvas = new TCanvas("provaCanvas");
  TH1F *subtraction_invMass_PK = new TH1F("subtraction_invMass_PK", "subtraction_invMass_PK", 80, 0, 2);
  TH1F *subtraction_invMass = new TH1F("subtraction_invMass", "subtraction_invMass", 80, 0, 2);
  subtraction_invMass_PK->Sumw2();
  subtraction_invMass_PK->Add(mass3Histo, mass4Histo, 1, -1);
  subtraction_invMass->Sumw2();
  subtraction_invMass->Add(mass1Histo, mass2Histo, 1, -1);

  thetaHisto->Fit("pol0", "Q");
  phiHisto->Fit("pol0", "Q");
  impulseHisto->Fit("expo", "Q");
  subtraction_invMass_PK->Fit("gaus", "Q", "", 0.6, 1.2);
  subtraction_invMass->Fit("gaus", "Q", "", 0.6, 1.2);
  mass5Histo->Fit("gaus", "Q", "", 0.6, 1.2);

  subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral());
  subtraction_invMass->SetEntries(subtraction_invMass->Integral());
}

  /*
    // ***** FITTING *****
    thetaHisto->Fit("pol0", "Q");
    phiHisto->Fit("pol0", "Q");
    impulseHisto->Fit("expo", "Q");
    subtractPK->Fit("gaus", "Q", "", 0.6, 1.2);
    subtractAll->Fit("gaus", "Q", "", 0.6, 1.2);
    mass5Histo->Fit("gaus", "Q", "", 0.6, 1.2);

    // CONTROLLA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TF1 *thetaFit = thetaHisto->GetFunction("thetaFit");
    TF1 *phiFit = phiHisto->GetFunction("phiFit");
    TF1 *impulseFit = impulseHisto->GetFunction("impulseFit");
    TF1 *subtractPKFit = subtractPK->GetFunction("subtractPKFit");
    TF1 *subtractAllFit = subtractAll->GetFunction("subtractAllFit");
    TF1 *mass5Fit = mass5Histo->GetFunction("mass5Fit");
  */
