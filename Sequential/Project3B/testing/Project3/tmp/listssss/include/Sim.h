#include "ParticleList.h"
#include "vec.h"

class Sim {
private:
  double (*PotentialEnergyFunction)(Vector3D rij);
  Vector3D (*ForceFunction)(Vector3D rij);
  Vector3D SimBox;
  double dt, KE, PE, T;
  int nAtoms;
  ParticleList *particles;

public:
  Sim();
  Sim(double _dt, Vector3D (*_Force)(Vector3D rij),
      double (*_Potential)(Vector3D rij), double _r_verlet);
  ~Sim() {}

  void InitSim(Vector3D _SimBox, Vector3D (*InitPositions)(int),
               Vector3D (*InitVelocities)(int));
};

void Sim::InitSim(Vector3D _SimBox, Vector3D (*InitPositions)(int) = 0,
                  Vector3D (*InitVelocities)(int) = 0) {
  int i = 0;
  ParticleNode *entry;
  for (entry = this->particles->GetFirst(); entry; entry = entry->GetNext())
    entry->GetParticle->r = InitPositions(i++);
}
