#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <stdexcept>
#include <iostream>

class Particle
{
public:
  //Particle() = delete;
  Particle(std::string name, double Px = 0., double Py = 0., double Pz = 0.);

  double GetPx() const;
  double GetPy() const;
  double GetPz() const;

  void SetP(double Px, double Py, double Pz);

  double GetMass() const;

  int GetIndex() const;

  void SetIndex(int index);
  void SetIndex(std::string name);

  double ComputeEnergy() const;
  double InvMass(Particle &p) const;

  static void AddParticleType(std::string name, double mass, int charge, double width = 0.);

  static void PrintList();
  void PrintData() const;

private:
  int fIndex;
  double fPx, fPy, fPz;

  static const int fMaxNumParticleType = 10; // max array dimension
  static ParticleType *fParticleType[fMaxNumParticleType];
  static int fNParticleType; // counter of array's elements

  static int FindParticle(std::string name);
};

#endif