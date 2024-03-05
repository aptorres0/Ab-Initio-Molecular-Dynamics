#include "Ensemble.h"
#include "Lennard-Jones.h"
#include "Sim.h"
#include "vec.h"
#include <cmath>

int main() {
  double L = cbrt(108.0 / 0.8);
  VectorD SimBox = {L, L, L};
  Sim sim;
  sim.InitSim(3, 108, LennardJonesPotential, LennardJonesForce, SimBox, 0.01);
  sim.InitPositions();
  sim.InitVelocities();

  sim.PrintSample();
  for (int i = 0; i < 2'000; ++i) {
    sim.VelocityVerletStep();
    sim.PrintSample();
  }
}