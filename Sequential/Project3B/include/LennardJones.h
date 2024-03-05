#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "Vector3.h"

double LennardJonesPotential(const Vector3D &rij) {
  return (4.0 / pow6(rij)) * (1.0 / pow6(rij) - 1.0);
}

// function called a lot, so it's worth optimizing
Vector3D LennardJonesForce(const Vector3D &rij) {
  // return (48.0f / pow8(rij)) * (1.0f / pow6(rij) - 0.5f) * rij;
  double factor = (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5);
  Vector3D result{};
  result[0] = factor * rij[0];
  result[1] = factor * rij[1];
  result[2] = factor * rij[2];
  return result;
}

#endif
