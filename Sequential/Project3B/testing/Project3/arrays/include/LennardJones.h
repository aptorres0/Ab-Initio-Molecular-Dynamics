#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "vec.h"

double LennardJonesPotential(const VectorD &rij);

VectorD LennardJonesForce(const VectorD &rij);
// void LennardJonesForce(VectorD &_force, const VectorD &rij);

/*VectorD gFCC[108];

VectorD FCC(int i) { return gFCC[i]; }

VectorD gV0[108];

VectorD v0(int i) { return gV0[i]; }

void GenerateFCC() {
  double l = cbrt(108.0 / 0.8) / 3.0;
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * l;
        gFCC[i * 36 + j * 12 + k * 4] = VectorD(xmin, ymin, zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 1] =
            VectorD(xmin, 0.5 * l + ymin, 0.5 * l + zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 2] =
            VectorD(0.5 * l + xmin, ymin, 0.5 * l + zmin);
        gFCC[i * 36 + j * 12 + k * 4 + 3] =
            VectorD(0.5 * l + xmin, 0.5 * l + ymin, zmin);
      }
    }
  }
}*/

#endif
