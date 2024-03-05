#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "ParticleList.h"
#include "vec.h"

class Particle {
public:
  VectorD r, v, a;
  ParticleList *VerletList;

  Particle() : r{}, v{}, a{}, VerletList{NULL} {}
  Particle(int _DIM) : r{_DIM}, v{_DIM}, a{_DIM}, VerletList{NULL} {}
  Particle(VectorD _r, VectorD _v) : r{_r}, v{_v}, a{}, VerletList{NULL} {}
}

#endif