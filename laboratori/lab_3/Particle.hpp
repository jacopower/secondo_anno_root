#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

class Particle
{
public:
  Particle();
  Particle(std::string name, double Px = 0., double Py = 0., double Pz = 0.);

  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  void SetP(double Px, double Py, double Pz);

  double GetMass() const;
  double GetCharge() const;

  int GetIndex() const;
  void SetIndex(int index);
  void SetIndex(std::string name);

  double ComputeEnergy() const;
  double InvMass(Particle &p) const;

  static void AddParticleType(std::string name, double mass, int charge, double width = 0.);

  static void PrintList();
  void PrintData() const;

  int Decay2body(Particle &dau1, Particle &dau2) const;

private:
  int fIndex;
  double fPx, fPy, fPz;

  static const int fMaxNumParticleType = 10;
  static ParticleType *fParticleType[fMaxNumParticleType];
  static int fNParticleType;

  static int FindParticle(std::string name);
  void Boost(double bx, double by, double bz);
};

#endif