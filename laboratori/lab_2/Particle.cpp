#include "Particle.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>

Particle::Particle() : fIndex{-1}, fPx{0.}, fPy{0.}, fPz{0.} {}

Particle::Particle(std::string name, double Px /*= 0.*/, double Py /*= 0.*/, double Pz /*= 0.*/)
{
  int const index = FindParticle(name);
  if (index == -1)
  {
    std::cerr << "Can't find correspondence" << '\n';
  }

  fIndex = index;
  fPx = Px;
  fPy = Py;
  fPz = Pz;
}

double Particle::GetPx() const { return fPx; }
double Particle::GetPy() const { return fPy; }
double Particle::GetPz() const { return fPz; }

void Particle::SetP(double Px, double Py, double Pz)
{
  fPx = Px;
  fPy = Py;
  fPz = Pz;
}

double Particle::GetMass() const
{
  return fParticleType[fIndex]->GetMass(); // ParticleType::GetMass()
}

double Particle::GetCharge() const
{
  return fParticleType[fIndex]->GetCharge(); // ParticleType::GetMass()
}

int Particle::GetIndex() const { return fIndex; }

void Particle::SetIndex(int index)
{
  if (index < 0 || index >= fNParticleType)
  {
    std::cerr << "This kind of particle doesn't exist" << '\n';
  }
  else
  {
    fIndex = index;
  }
}

void Particle::SetIndex(std::string name)
{
  int const index = FindParticle(name);

  if (index == -1)
  {
    std::cerr << "This kind of particle doesn't exists!" << '\n';
  }
  else
  {
    fIndex = index;
  }
}

double Particle::ComputeEnergy() const
{
  double const pSquared = fPx * fPx + fPy * fPy + fPz * fPz;
  double const mass = GetMass(); // Particle::GetMass()
  return std::sqrt(mass * mass + pSquared);
}

double Particle::InvMass(Particle &p) const // FUNZIONA SENZA THIS? OVVIAMENTE SI
{
  double E1 = ComputeEnergy(); // this->ComputeEnergy()
  double E2 = p.ComputeEnergy();
  double PxSum = fPx + p.fPx;
  double PySum = fPy + p.fPy;
  double PzSum = fPz + p.fPz;

  double PSumSquared = PxSum * PxSum + PySum * PySum + PzSum * PzSum;
  double result = sqrt((E1 + E2) * (E1 + E2) - PSumSquared);
  return result;
}

void Particle::AddParticleType(std::string name, double mass, int charge, double width /*= 0.*/)
{
  if (fNParticleType == fMaxNumParticleType)
  {
    std::cerr << "Can't add another particle: maximum size reached up!" << '\n';
    return;
  }

  if (FindParticle(name) == -1)
  {
    if (width == 0) // CHIEDI PER QUESTO
    {
      fParticleType[fNParticleType] = new ParticleType(name, mass, charge);
      fNParticleType++;
    }
    else
    {
      fParticleType[fNParticleType] = new ResonanceType(name, mass, charge, width);
      fNParticleType++;
    }
  }
  else
  {
    std::cerr << "Particle " << name << " already exists!" << '\n';
  }
}

void Particle::PrintList()
{
  for (int i = 0; i < fNParticleType; ++i)
  {
    std::cout << "Particle " << i << '\n';
    fParticleType[i]->Print();
    std::cout << '\n';
  }
}

void Particle::PrintData() const
{
  std::cout << "Index: " << fIndex << '\n'
            << "Name: " << fParticleType[fIndex]->GetName() << '\n'
            << "Px: " << fPx << '\n'
            << "Py: " << fPy << '\n'
            << "Pz: " << fPz << '\n';
}

// CHIEDI SPIEGAZIONI E SCOMMENTA Y2 -> E' UNUSED PARAMETER
int Particle::Decay2body(Particle &dau1, Particle &dau2) const
{
  if (GetMass() == 0.0)
  {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex > -1)
  {
    float x1, x2, w, y1 /*, y2*/;

    double invnum = 1. / RAND_MAX;
    do
    {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    //y2 = x2 * w;

    massMot += fParticleType[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2)
  {
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }

  double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

ParticleType *Particle::fParticleType[fMaxNumParticleType] = {};
int Particle::fNParticleType = 0;

int Particle::FindParticle(std::string name) // CHIEDI PER CERR + posso partire da zero?
{
  for (int index = 0; index < fNParticleType; ++index)
  {
    if (fParticleType[index]->GetName() == name)
    {
      return index;
    }
  }

  //std::cerr << "Can't find correspondence" << '\n'; // l'ho levato
  return -1;
}

void Particle::Boost(double bx, double by, double bz)
{
  double energy = ComputeEnergy();

  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}