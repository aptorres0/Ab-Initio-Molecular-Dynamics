#include "Ar/system.h"
#include "Vector3D.h"
#include <cstdio>
#include <iostream>

int main() {
  /*FILE *fsteps;
  fsteps = fopen("steps.txt", "w");
  FILE *fEn;
  fEn = fopen("energy.txt", "w");*/
  Ar::System sys;
  sys.InitPositions();
  sys.InitVelocities();
  for (int i = 0; i < 2'000'000; ++i) {
    /*  for (int j = 0; j < 108; ++j) {
        fprintf(fsteps, "%f %f %f ", sys.q(j).GetX(), sys.q(j).GetY(),
                sys.q(j).GetZ());
      }
      fprintf(fsteps, "\n");*/
    sys.Step();
    // fprintf(fEn, "%d %f %f %f\n", i, sys.GetE(), sys.GetT(), sys.GetU());
  }

  /*fclose(fsteps);
  fclose(fEn);*/

  return 0;
}