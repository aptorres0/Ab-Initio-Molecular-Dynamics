#include "LennardJones.h"
#include "Sim.h"
#include "vec.h"
#include <cmath>
#include <cstdio>

int main() {
  FILE *fptr;
  fptr = fopen("energies.txt", "w");
  double L = cbrt(108.0 / 0.8);
  VectorD SimBox = {L, L, L};
  Sim sim;
  sim.PrintSample();
  sim.InitSim(3, 108, LennardJonesPotential, LennardJonesForce, SimBox, 0.01);
  sim.PrintSample();
  sim.InitPositions();
  sim.PrintSample();
  sim.InitVelocities();
  sim.PrintSample();

  for (int i = 0; i < 2'000; ++i) {
    sim.VelocityVerletStep();
    /*fprintf(fptr, "%d %f %f %f\n", i, sim.GetInternalEnergy(),
            sim.GetKineticEnergy(), sim.GetPotentialEnergy());*/
    std::cout << sim.GetInternalEnergy() << " " << sim.GetKineticEnergy() << " "
              << sim.GetPotentialEnergy() << '\n';
    // sim.PrintSample();
  }
  fclose(fptr);

  return 0;
}