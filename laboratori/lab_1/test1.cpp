#include "ParticleType.hpp"
#include "ResonanceType.hpp"

#include <iostream>

int main()
{
  ParticleType particle{"Nome Particle", 124.4, 122};
  ResonanceType resonance{"Nome Resonance", 1222.001, 10, 123.2};

  particle.Print();
  resonance.Print();

  std::cout << "Particle Name: " << particle.GetName() << '\n'
            << "Particle Mass: " << particle.GetMass() << '\n'
            << "Particle Charge: " << particle.GetCharge() << '\n';

  std::cout << "Resonance Name: " << resonance.GetName() << '\n'
            << "Resonance Mass: " << resonance.GetMass() << '\n'
            << "Resonance Charge: " << resonance.GetCharge() << '\n'
            << "Resonance Width: " << resonance.GetWidth() << '\n';

  std::cout << '\n'
            << "-----------------------------" << '\n'
            << '\n';

  ParticleType *provaArray[2];
  provaArray[0] = new ParticleType{"Particle Array", 678.2, 87};
  provaArray[1] = new ResonanceType{"Resonance Array", 455.342, 46, 234.2};

  for (int i = 0; i < 2; ++i)
  {
    provaArray[i]->Print();
  }

  //delete[] provaArray;
}