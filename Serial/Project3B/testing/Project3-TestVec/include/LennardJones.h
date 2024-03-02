#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "vec.h"

double LennardJonesPotential(const Vector3D &rij) {
  return (4.0 / pow6(rij)) * (1.0 / pow6(rij) - 1.0);
}

Vector3D LennardJonesForce(const Vector3D &rij) {
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}

#endif
