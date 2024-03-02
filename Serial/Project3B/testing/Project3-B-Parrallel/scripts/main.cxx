#include "../include/LennardJones.h"
#include "../include/Sim.h"
#include "../include/vec.h"
#include <cmath>
#include <cstdio>

void SimpleSim();
void GenerateAnimationData();

int main() {
  SimpleSim();
  // GenerateAnimationData();
  return 0;
}

void SimpleSim() {
  // Define the simulation cell
  float L = cbrt(108.0 / 0.8);
  Vector3F SimBox = {L, L, L};

  // Define arguments to simulation class
  float h =
      0.008311612f; // Time step size = 1% of the period of oscillation for the
                    // harmonic oscillator in the lennard jones potential
  SimArgs args = {.SimulationCellBounds = SimBox,
                  .TimeStepSize = h,
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

  for (int i = 0; i < 6'000; ++i) {
    sim.VelocityVerletStep();
    /*std::cout << sim.GetVelocityDistribution() << " " << sim.GetTemperature()
              << " " << sim.GetInternalEnergy() << " " << sim.GetKineticEnergy()
              << " " << sim.GetPotentialEnergy() << '\n';*/
  }
  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';
}

void GenerateAnimationData() {
  // Create output text file to write out energies
  FILE *fptr;
  fptr = fopen("pos.txt", "w");

  // Define the simulation cell
  float L = cbrt(108.0 / 0.8);
  std::cout << "FCC cell l = " << L / 3.0 << '\n';
  Vector3F SimBox = {L, L, L};

  // Define arguments to simulation class
  float h =
      0.008311612f; // Time step size = 1% of the period of oscillation for the
                    // harmonic oscillator in the lennard jones potential
  SimArgs args = {.SimulationCellBounds = SimBox,
                  .TimeStepSize = h,
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

  for (int j = 0; j < 108; ++j) {
    Vector3F tmp = sim.GetPos(j);
    fprintf(fptr, "%f %f %f ", tmp[0], tmp[1], tmp[2]);
  }
  fprintf(fptr, "\n");
  for (int i = 0; i < 2'000; ++i) {
    sim.VelocityVerletStep();
    for (int j = 0; j < 108; ++j) {
      Vector3F tmp = sim.GetPos(j);
      fprintf(fptr, "%f %f %f ", tmp[0], tmp[1], tmp[2]);
    }
    fprintf(fptr, "\n");
    std::cout << sim.GetVelocityDistribution() << " " << sim.GetTemperature()
              << " " << sim.GetInternalEnergy() << " " << sim.GetKineticEnergy()
              << " " << sim.GetPotentialEnergy() << '\n';
  }
  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';

  fclose(fptr);
}