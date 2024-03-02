#include "Vector3D.h"

double Potential(Vector3D rij) {
  return (4.0 / pow6(rij)) * (1.0 / pow6(rij) - 1.0);
}

Vector3D Force(Vector3D rij) {
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}

Vector3D gFCC[108];

Vector3D FCC(int i) { return gFCC[i]; }

Vector3D gV0[108];

Vector3D v0(int i) { return gV0[i]; }

void GenerateFCC() {
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * l;
        gFCC[i * 36 + j * 12 + k * 4] = Vector3D(xmin, ymin, zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 1] =
            Vector3D(xmin, 0.5 * l + ymin, 0.5 * l + zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 2] =
            Vector3D(0.5 * l + xmin, ymin, 0.5 * l + zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 3] =
            Vector3D(0.5 * l + xmin, 0.5 * l + ymin, zmin);
      }
    }
  }
}

void System::InitVelocities() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless variance
   *     of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // initialize the rng
  init_random();

  // 1. set velocities from gaussian and compute the resulting momentum
  Vector3D lNetMomentum(0.0, 0.0, 0.0);
  for (int i = 0; i < nAtoms; ++i) {
    /*atom[i].velocity = Vector3D(rnd(), rnd(), rnd());
    lNetMomentum += atom[i].velocity;*/
    p(i) = Vector3D(rnd(), rnd(), rnd());
    lNetMomentum += p(i);
  }
  // save that random state
  save_random();

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  double Temp = 0;
  for (int i = 0; i < nAtoms; ++i) {
    /*atom[i].velocity -= lNetMomentum / ((double)nAtoms);
    Temp += atom[i].velocity.sq();*/
    p(i) -= lNetMomentum / ((double)nAtoms);
    Temp += p(i).sq();
  }
  Temp /= 3.0 * (nAtoms - 1.0);

  // 3. set temperature of the system by multiplying by a correction factor
  double T_corr = sqrt(1.0 / Temp);
  for (int i = 0; i < nAtoms; ++i) {
    // atom[i].velocity *= T_corr;
    p(i) *= T_corr;
  }
}