#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "ParticleList.h"
#include "vec.h"

class Particle {
public:
  Vector3D r, v, a;
  ParticleList *VerletList;

  Particle() : r{}, v{}, a{}, VerletList{NULL} {}
  Particle(Vector3D _r, Vector3D _v) : r{_r}, v{_v}, a{}, VerletList{NULL} {}
}

#endif