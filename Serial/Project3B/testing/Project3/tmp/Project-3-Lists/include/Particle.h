#ifndef _PARTICLE_H
#define _PARTICLE_H

#ifndef _VEC_H
#include "vec.h"
#endif

class ParticleList;

class Particle {
public:
  // Constructors
  Particle(int dim)
      : r(0, dim), v(0, dim), a(0, dim), verletList(new ParticleList()) {}

  Particle(vec r_, vec v_)
      : dim(r_.get_dim()), r(r_), v(v_), verletList(new ParticleList()) {}

  // Destructor
  ~Particle();

  int dim;
  vec r, v, a;
  ParticleList *verletList;
};

#endif