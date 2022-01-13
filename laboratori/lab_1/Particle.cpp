#include "Particle.hpp"

#include <iostream>
#include <cmath>

Particle::Particle(std::string name, double Px, double Py, double Pz) // rivedi
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
  return fParticleType[fIndex]->ParticleType::GetMass(); // ridondante?
}

int Particle::GetIndex() const { return fIndex; }

void Particle::SetIndex(int index)
{
  if (index < 0 || index > (fNParticleType - 1)) // Controlla sia giusto
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
  auto const index = FindParticle(name);

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
  double pSquared = fPx * fPx + fPy * fPy + fPz * fPz;
  //double mass = this->GetMass(); // giusto?
  double const mass = GetMass();
  return std::sqrt(mass * mass * pSquared);
}

double Particle::InvMass(Particle &p) const
{
  double const E1 = this->ComputeEnergy();
  double const E2 = p.ComputeEnergy();
  double ESumSquared = (E1 + E2) * (E1 + E2);
  double PxSum = this->fPx + p.fPx;
  double PySum = this->fPy + p.fPy;
  double PzSum = this->fPz + p.fPz;
  double PSumSquared = PxSum * PxSum + PySum * PySum + PzSum * PzSum;
  return ESumSquared * PSumSquared;
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
    if (width == 0)
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
    std::cout << "PARTICLE " << i << '\n';
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

ParticleType *Particle::fParticleType[fMaxNumParticleType] = {}; // controlla inizializzazione
int Particle::fNParticleType = 0;

int Particle::FindParticle(std::string name)
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