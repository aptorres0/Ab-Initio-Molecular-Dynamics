#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "Vector3.h"

float LennardJonesPotential(const Vector3F &rij) {
  return (4.0f / pow6(rij)) * (1.0f / pow6(rij) - 1.0f);
}

// function called a lot, so it's worth optimizing
Vector3F LennardJonesForce(const Vector3F &rij) {
  // return (48.0f / pow8(rij)) * (1.0f / pow6(rij) - 0.5f) * rij;
  float factor = (48.0f / pow8(rij)) * (1.0f / pow6(rij) - 0.5f);
  Vector3F result{};
  result[0] = factor * rij[0];
  result[1] = factor * rij[1];
  result[2] = factor * rij[2];
  return result;
}

#endif
