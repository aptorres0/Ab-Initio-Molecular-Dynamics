#include "../include/Sim.h"
#include "LennardJones.h"
#include "vec.h"
#include <cmath>
#include <cstdio>

int main() {
  // Create output text file to write out energies
  /*FILE *fptr;
  fptr = fopen("energies.txt", "w");
  fprintf(fptr, "#totalenergy kineticenergy potentialenergy\n");*/

  // Define the simulation cell
  double L = cbrt(108.0 / 0.8);
  VectorD SimBox = {L, L, L};

  // Define arguments to simulation class
  SimArgs args = {.NumberOfDimensions = 3,
                  .NumberOfAtoms = 108,
                  .SimulationCellBounds = SimBox,
                  .TimeStepSize = 0.01,
                  .PotentialEnergyFunction = LennardJonesPotential,
                  .ForceFunction = LennardJonesForce};

  // Declare an instance of the simulation class
  Sim sim;
  // initialize arguments
  sim.InitSim(args);
  // initialize positions of atoms to FCC crystal lattice
  sim.InitPositions();
  // initialize the random velocities for a system with 0 net momentum and at
  // 120 K
  sim.InitVelocities();
  std::cout << "Initial Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nInitial Temperature [reduced units] = "
            << sim.GetTemperature() << '\n';

  for (int i = 0; i < 20'000; ++i) {
    sim.VelocityVerletStep();
    /*fprintf(fptr, "%d %f %f %f\n", i, sim.GetInternalEnergy(),
            sim.GetKineticEnergy(), sim.GetPotentialEnergy());*/
    /*std::cout << sim.GetInternalEnergy() << " " << sim.GetKineticEnergy()
              << "
                 "
              << sim.GetPotentialEnergy() << '\n';*/
  }
  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';

  // fclose(fptr);

  return 0;
}