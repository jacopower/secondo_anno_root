#include "ParticleType.hpp"

#include <iostream>
//#include <string>
#include <stdexcept>

ParticleType::ParticleType(std::string name, double mass, int charge) : fName{name}, fMass{mass}, fCharge{charge}
{
  if (mass < 0)
  {
    throw std::runtime_error("Mass can't be negative!");
  }
  // vedi se ci sono altre resitrizioni e se mettere un cerr magari
}

std::string ParticleType::GetName() const { return fName; }
double ParticleType::GetMass() const { return fMass; }
int ParticleType::GetCharge() const { return fCharge; }

void ParticleType::Print() const
{
  std::cout << "Name: " << fName << '\n'
            << "Mass: " << fMass << '\n'
            << "Charge: " << fCharge << '\n';
}