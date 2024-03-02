#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "ParticleList.h"
#include "Vector3D.h"

// NOTE: Using `new` allocated memory on the heap

// Data element of node (ParticleListEntry)
class Particle {
public:
  // Constructors
  Particle() : r{}, v{}, a{} {}

  Particle(Vector3D _r, Vector3D _v) : r{_r}, v{_v} {}

  // Destructor
  ~Particle();

  Vector3D r, v, a;
};

#endif