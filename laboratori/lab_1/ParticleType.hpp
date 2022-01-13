#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP

#include <string>

class ParticleType
{
public:
  ParticleType(std::string name, double mass, int charge);

  std::string GetName() const;
  double GetMass() const;
  int GetCharge() const;
  virtual void Print() const;

private:
  const std::string fName;
  const double fMass;
  const int fCharge;
};

#endif