#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"
#include "TROOT.h" // for R__LOAD_LIBRARY (LEVAAAA)
#include "TH1F.h"
#include "TFile.h"

R__LOAD_LIBRARY(ParticleType_cpp.so)
R__LOAD_LIBRARY(ResonanceType_cpp.so)
R__LOAD_LIBRARY(Particle_cpp.so)

#include "TRandom.h"
#include "TMath.h"
#include <iostream>

int main()
{
  gRandom->SetSeed();
  constexpr Int_t nGen = 1E3;
  constexpr Int_t nParticles = 100;

  // CREO I TIPI DI PARTICELLE
  Particle::AddParticleType("Pi+", 0.13957, 1);
  Particle::AddParticleType("Pi-", 0.13957, -1);
  Particle::AddParticleType("K+", 0.49367, 1);
  Particle::AddParticleType("K-", 0.49367, -1);
  Particle::AddParticleType("P+", 0.93827, 1);
  Particle::AddParticleType("P-", 0.93827, -1);
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);

  TFile *file = new TFile("mySimulation.root", "RECREATE");

  //CREO GLI ISTOGRAMMI ROOT
  TH1F *particlesHisto = new TH1F("particlesHisto", "Types of Particles Generated", 7, 0, 7);
  TH1F *phiHisto = new TH1F("phiHisto", "Distribution of Azimutal Angle", 200, 0, 2 * TMath::Pi());                // 100
  TH1F *thetaHisto = new TH1F("thetaHisto", "Distribution of Polar Angle", 200, 0, TMath::Pi());                   // 100
  TH1F *impulseHisto = new TH1F("impulseHisto", "Impulse Distribution", 200, 0, 5);                                // 100
  TH1F *transverseImpulseHisto = new TH1F("transverseImpulseHisto", "Transverse Impulse Distribution", 200, 0, 5); // 100
  TH1F *energyHisto = new TH1F("energyHisto", "Distribution of Particles' Energy", 200, 0, 5);                     // 100

  TH1F *mass1Histo = new TH1F("mass1Histo", "Invariant Mass - All Particles", 160, 0, 2);              // 160
  TH1F *mass2Histo = new TH1F("mass2Histo", "Invariant Mass - Discord Charges", 160, 0, 2);            // 160
  TH1F *mass3Histo = new TH1F("mass3Histo", "Invariant Mass - Concord Charges", 160, 0, 2);            // 160
  TH1F *mass4Histo = new TH1F("mass4Histo", "Invariant Mass - P & K with Discord Charges", 160, 0, 2); // 160
  TH1F *mass5Histo = new TH1F("mass5Histo", "Invariant Mass - P & K with Concord Charges", 160, 0, 2); // 160
  TH1F *mass6Histo = new TH1F("mass6Histo", "Invariant Mass - K* Decays", 320, 0, 2);                  // 160

  // FACCIO I MIEI CONTICINI
  Particle particles[120];

  Double_t phi, theta = 0.;
  Double_t px, py, pz = 0.;
  Double_t impulse = 0.;

  Double_t randomParticle, randomDecay = 0.;

  for (Int_t i = 0; i < nGen; ++i) // 1E5 events
  {
    Int_t decaysCounter = 0.;

    for (Int_t j = 0; j < nParticles; ++j) // 100 particles
    {
      phi = gRandom->Uniform(0, 2 * TMath::Pi()); // uniform [0, 2pi]
      phiHisto->Fill(phi);

      theta = gRandom->Uniform(0, TMath::Pi()); // uniform [0, pi]
      thetaHisto->Fill(theta);

      impulse = gRandom->Exp(1); // exponential e^-x
      impulseHisto->Fill(impulse);

      px = impulse * sin(theta) * cos(phi);
      py = impulse * sin(theta) * sin(phi);
      pz = impulse * cos(theta);
      particles[j].SetP(px, py, pz);
      transverseImpulseHisto->Fill(sqrt(px * px + py * py));

      randomParticle = gRandom->Rndm(); // uniform [0,1]

      if (randomParticle < 0.4) // 40% -> Pi+
      {
        particles[j].SetIndex("Pi+");
      }
      else if (randomParticle < 0.8) // 40% -> Pi-
      {
        particles[j].SetIndex("Pi-");
      }
      else if (randomParticle < 0.85) // 5% -> K+
      {
        particles[j].SetIndex("K+");
      }
      else if (randomParticle < 0.90) // 5% -> K-
      {
        particles[j].SetIndex("K-");
      }
      else if (randomParticle < 0.945) // 4.5% -> P+
      {
        particles[j].SetIndex("P+");
      }
      else if (randomParticle < 0.99) // 4.5% -> P-
      {
        particles[j].SetIndex("P-");
      }
      else // 1% -> K*
      {
        particles[j].SetIndex("K*");

        randomDecay = gRandom->Rndm();
        if (randomDecay < 0.5) // K* -> Pi+ & K-
        {
          particles[100 + decaysCounter].SetIndex("Pi+");
          particles[100 + decaysCounter + 1].SetIndex("K-");
        }
        else // K* -> Pi- & K+
        {
          particles[100 + decaysCounter].SetIndex("Pi-");
          particles[100 + decaysCounter + 1].SetIndex("K+");
        }

        // Assegno le giuste quantità di moto dopo il decadimento:
        particles[j].Decay2body(particles[100 + decaysCounter], particles[100 + decaysCounter + 1]);
        mass6Histo->Fill(particles[100 + decaysCounter].InvMass(particles[100 + decaysCounter + 1]));
        decaysCounter += 2;
      } // fine else k*

      // RIEMPIO ISTOGRAMMA PARTICELLE GENERATE E ENERGIA
      particlesHisto->Fill(particles[j].GetIndex());
      energyHisto->Fill(particles[j].ComputeEnergy());

    } // fine primo for (j)

    // RIEMPIO GLI ISTROGRAMMI DELLE MASSE INVARIANTI
    for (Int_t j = 0; j < (100 + decaysCounter); ++j)
    {
      for (Int_t k = j + 1; k < (100 + decaysCounter); ++k)
      {
        Double_t relativeCharge = particles[j].GetCharge() * particles[k].GetCharge();
        Bool_t discordPK = (particles[j].GetIndex() == 0 && particles[k].GetIndex() == 3) || (particles[j].GetIndex() == 1 && particles[k].GetIndex() == 2);
        Bool_t concordPK = (particles[j].GetIndex() == 0 && particles[k].GetIndex() == 2) || (particles[j].GetIndex() == 1 && particles[k].GetIndex() == 3);

        mass1Histo->Fill(particles[j].InvMass(particles[k]));

        if (relativeCharge == -1)
        {
          mass2Histo->Fill(particles[j].InvMass(particles[k]));
        }
        if (relativeCharge == 1)
        {
          mass3Histo->Fill(particles[j].InvMass(particles[k]));
        }
        if (discordPK)
        {
          mass4Histo->Fill(particles[j].InvMass(particles[k]));
        }
        if (concordPK)
        {
          mass5Histo->Fill(particles[j].InvMass(particles[k]));
        }
      } // fine primo for
    }   // fine secondo for
    // FINE RIEMPIMENTO MASSE INVARIANTI

  } // fine secondo for (i)

  file->Write();
  file->Close();
}