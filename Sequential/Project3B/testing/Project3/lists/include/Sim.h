#include "ParticleList.h"
#include "vec.h"

class Sim {
private:
  double (*U_ij)(const VectorD &);
  VectorD (*F_ij)(const VectorD &);
  VectorD (*r0)(int);
  VectorD (*v0)(int);

  int fDIM;
  int fNumParticles;
  VectorD SimBox{fDIM};
  ParticleList *particles;
  double h;

public:
  Sim() : fDIM{3}, fNumParticles{108}, h{0.01}, particles{new ParticleList()} {
    for (int i = 0; i < fNumParticles; ++i) {
      this->particles->AddParticle(new Particle(fDIM));
    }
  }
  void Setup(VectorD (*_r0)(int), VectorD (*_v0)(int), VectorD _SimBox,
             VectorD (*_F)(const VectorD &), double (*_U)(const VectorD &)) {
    U_ij = _U;
    F_ij = _F;
    SimBox = _SimBox;
    ParticleNode *entry;
    if (_r0 && _v0) {
      int i = 0;
      for (entry = this->particles->GetFirst(); entry; entry->GetNext()) {
        entry->GetParticle()->r = _r0(i);
        entry->GetParticle()->v = _v0(i++);
      }
    }
  }
  Sim(double _dt, VectorD (*_Force)(VectorD rij),
      double (*_Potential)(VectorD rij), double _r_verlet);
  ~Sim() {}
};
