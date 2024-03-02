#include "LennardJones.h"

double LennardJonesPotential(const VectorD &rij) {
  return (4.0 / pow6(rij)) * (1.0 / pow6(rij) - 1.0);
}

// Return by reference
/*void LennardJonesForce(VectorD &_force, const VectorD &rij) {
  _force.assign_to((48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij,
                   rij.length()); // avoid a copy
}*/

VectorD LennardJonesForce(const VectorD &rij) {
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}