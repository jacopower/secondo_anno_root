#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "Particle.hpp"

int main()
{
  // Add Particle Type
  Particle::AddParticleType("protone", 1.013, 2);
  Particle::AddParticleType("elettrone", 3, 4, 5);
  Particle::AddParticleType("neutrone", 6, 7);

  // Particle List (static array)
  std::cout << "-------------------------------------------------------------" << '\n'
            << "LISTA PARTICELLE" << '\n';
  Particle::PrintList();
  std::cout << "-------------------------------------------" << '\n'
            << '\n';

  // Constructor errors
  std::cout << "ERRORI COSTRUTTORE" << '\n';
  Particle particella_errore("ciao");
  std::cout << "particella_errore index: " << particella_errore.GetIndex() << '\n'
            << "---------------------------------" << '\n'
            << '\n';

  // create Particle
  Particle particella_p("protone", 12, 23, 43);
  Particle particella_e("elettrone", 12, 32, 45);
  Particle particella_n("neutrone", 14, 35, 75);

  // Index
  std::cout << "INDICI" << '\n'
            << "particella_p: " << particella_p.GetIndex() << '\n'
            << "particella_e: " << particella_e.GetIndex() << '\n'
            << "particella_n: " << particella_n.GetIndex() << '\n'
            << "---------------------------------------------" << '\n'
            << '\n';

  // Testing Getters
  std::cout << "PROVA GETTER PARTICELLA_P" << '\n'
            << "Px: " << particella_p.GetPx() << '\n'
            << "Py: " << particella_p.GetPy() << '\n'
            << "Pz: " << particella_p.GetPz() << '\n'
            << "------------------------------------------------------" << '\n'
            << '\n';

  // Testing Set/Get P
  particella_p.SetP(7, 8, 9);
  std::cout << "PROVA SETTER PARTICELLA_P" << '\n'
            << "Px: " << particella_p.GetPx() << '\n'
            << "Py: " << particella_p.GetPy() << '\n'
            << "Pz: " << particella_p.GetPz() << '\n'
            << "------------------------------------------------------" << '\n'
            << '\n';

  // Testing Set Mass / Index
  std::cout << "PROVA GETTERS" << '\n'
            << "particella_p mass: " << particella_p.GetMass() << '\n'
            << "particella_p index: " << particella_p.GetIndex() << '\n'
            << "------------------------------" << '\n'
            << '\n';

  // Testing Setter Index
  std::cout << "PROVE SETTER INDEX" << '\n';
  particella_p.SetIndex(2);
  std::cout << "particella_p index: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("elettrone");
  std::cout << "particella_p index: " << particella_p.GetIndex() << '\n';
  std::cout << "-------------------------------------" << '\n'
            << '\n';

  // Testing Setter Index Errors

  // Test (1)
  std::cout << "PROVA ERRORI SETTER INDEX" << '\n'
            << "particella_p index: " << particella_p.GetIndex() << '\n'
            << "Provo a settare index a 3: " << '\n';
  particella_p.SetIndex(3);
  std::cout << "particella_p index: " << particella_p.GetIndex() << '\n'
            << '\n';

  // Test(2)
  std::cout << "Indice particella_p: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(-1);
  std::cout << "Settato indice a -1: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(0);
  std::cout << "Settato indice a 0: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(1);
  std::cout << "Settato indice a 1: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(2);
  std::cout << "Settato indice a 2: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(3);
  std::cout << "Settato indice a 3: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex(4);
  std::cout << "Settato indice a 4: " << particella_p.GetIndex() << '\n'
            << '\n';

  // test(3)
  std::cout << "Indice particella_p: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("ciao");
  std::cout << "Settato indice a ciao: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("protone");
  std::cout << "Settato indice a protone: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("elettrone");
  std::cout << "Settato indice a elettrone: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("neutrone");
  std::cout << "Settato indice a neutrone: " << particella_p.GetIndex() << '\n';
  particella_p.SetIndex("protone");
  std::cout << "Settato indice a protone: " << particella_p.GetIndex() << '\n';
  std::cout << "------------------------------------" << '\n'
            << '\n';

  // Testing Energy and InvMass
  std::cout << "PROVA ENERGIA E MASSA INVARIANTE" << '\n'
            << "particella_p energy: " << particella_p.ComputeEnergy() << '\n'
            << "particella_p invariant mass: " << particella_p.InvMass(particella_n) << '\n'
            << "------------------------------------" << '\n'
            << '\n';

  // Testing AddParticleType Errors
  std::cout << "ERRORI ADD PARTICLE TYPE (VEDI MAIN)" << '\n';
  Particle::AddParticleType("neutrone", 6, 7);
  Particle::AddParticleType("4", 12, 23);
  Particle::AddParticleType("5", 12, 23);
  Particle::AddParticleType("6", 12, 23);
  Particle::AddParticleType("7", 12, 23);
  Particle::AddParticleType("8", 12, 23);
  Particle::AddParticleType("9", 12, 23);
  Particle::AddParticleType("9", 12, 23);
  Particle::AddParticleType("10", 12, 23);
  Particle::AddParticleType("11", 123, 3, 12);
  std::cout << "LIST: " << '\n';
  Particle::PrintList();
  std::cout
      << "--------------------------------------" << '\n'
      << '\n';
}