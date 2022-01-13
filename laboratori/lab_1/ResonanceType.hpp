#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP

#include "ParticleType.hpp"
//#include <string>

class ResonanceType : public ParticleType
{
public:
  ResonanceType(std::string name, double mass, int charge, double width);

  double GetWidth() const;
  void Print() const override; // rivedi override

private:
  const double fWidth;
};

#endif