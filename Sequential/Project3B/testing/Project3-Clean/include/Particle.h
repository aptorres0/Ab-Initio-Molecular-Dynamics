#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "ParticleList.h"
#include "vec.h"

class Particle {
public:
  VectorD r, v, a;

  Particle() : r{}, v{}, a{} {}
  Particle(int _DIM) : r{_DIM}, v{_DIM}, a{_DIM} {}
  Particle(VectorD _r, VectorD _v) : r{_r}, v{_v}, a{} {}
}

#endif