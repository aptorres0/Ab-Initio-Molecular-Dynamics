#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "vec.h"

double LennardJonesPotential(const VectorD &rij) {
  return (4.0 / pow6(rij)) * (1.0 / pow6(rij) - 1.0);
}

VectorD LennardJonesForce(const VectorD &rij) {
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}

#endif
