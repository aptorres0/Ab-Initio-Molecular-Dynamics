#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "Vector3.h"

class Particle {
public:
  Vector3D r, v, a;
  Particle *next; // pointer to the next particle in the list
  int ID;

  Particle() : r{}, v{}, a{}, next{nullptr} {}
  Particle(Vector3D _r, Vector3D _v, int _ID = 0)
      : r{_r}, v{_v}, a{}, next{nullptr}, ID{_ID} {}
  Particle(Vector3D _r, Vector3D _v) : r{_r}, v{_v}, a{}, next{nullptr} {}
};

#endif