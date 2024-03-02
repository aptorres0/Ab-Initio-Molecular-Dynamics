#include "../include/Ar/constants.h"
#include "../include/LennardJones.h"
#include "../include/Particle.h"
#include "../include/Sim.h"
#include "../include/Vector3.h"
#include <cmath>
#include <cstdio>

void MicrocanonicalEnsembleSim();
void MacrocanonicalEnsembleSim();
void GenerateAnimationData();

int main() {
  MicrocanonicalEnsembleSim();
  MacrocanonicalEnsembleSim();
  // GenerateAnimationData();
  return 0;
}

void MicrocanonicalEnsembleSim() {

  FILE *fptr;
  fptr = fopen("MicrocanonicalEnsembleSim.txt", "w");
  fprintf(fptr, "#t E KE PE\n");

  // Define the simulation cell
  double L = cbrt(108.0 / 0.8);

  // Define arguments to simulation class
  double h =
      0.008311612f; // Time step size = 1% of the period of oscillation for the
                    // harmonic oscillator in the lennard jones potential
  SimArgs args = {.BoxLength = L,
                  .TimeStepSize = h,
                  .PotentialEnergyFunction = LennardJonesPotential,
                  .ForceFunction = LennardJonesForce};

  // Declare an instance of the simulation class
  Sim sim;
  // initialize arguments
  sim.InitSim(args);
  sim.InitPositions();
  sim.InitVelocities();
  std::cout << "Initial Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nInitial Temperature [reduced units] = "
            << sim.GetTemperature() << '\n';
  std::cout << "Velocity Distribution = " << sim.GetVelocityDistribution()
            << '\n';

  fprintf(fptr, "%d %f %f %f\n", 0, sim.GetInternalEnergy() * Ar::gEnergy,
          sim.GetKineticEnergy() * Ar::gEnergy,
          sim.GetPotentialEnergy() * Ar::gEnergy);

  for (int i = 0; i < 6'000; ++i) {
    sim.VelocityVerletStep();
    fprintf(fptr, "%f %f %f %f\n", (i + 1) * h * Ar::gTime,
            sim.GetInternalEnergy() * Ar::gEnergy,
            sim.GetKineticEnergy() * Ar::gEnergy,
            sim.GetPotentialEnergy() * Ar::gEnergy);
  }
  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';
  std::cout << "Velocity Distribution = " << sim.GetVelocityDistribution()
            << '\n';
}

void MacrocanonicalEnsembleSim() {

  // Define the simulation cell
  double L = cbrt(108.0 / 0.8);

  // Define arguments to simulation class
  double h =
      0.008311612f; // Time step size = 1% of the period of oscillation for the
                    // harmonic oscillator in the lennard jones potential
  SimArgs args = {.BoxLength = L,
                  .TimeStepSize = h,
                  .PotentialEnergyFunction = LennardJonesPotential,
                  .ForceFunction = LennardJonesForce};

  // Declare an instance of the simulation class
  Sim sim;
  // initialize arguments
  sim.InitSim(args);
  sim.InitPositions();
  sim.InitVelocities();
  std::cout << "Initial Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nInitial Temperature [reduced units] = "
            << sim.GetTemperature() << '\n';
  std::cout << "Velocity Distribution = " << sim.GetVelocityDistribution()
            << '\n';

  sim.PhaseI(6'000);

  std::cout << "Momentum after Phase I [reduced units] = "
            << sim.GetNetMomentum()
            << "\nTemperature after Phase I [reduced units] = "
            << sim.GetTemperature() << '\n';
  std::cout << "Velocity Distribution = " << sim.GetVelocityDistribution()
            << '\n';

  sim.PhaseII(2'000);

  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';
  std::cout << "Velocity Distribution = " << sim.GetVelocityDistribution()
            << '\n';
}

void GenerateAnimationData() {
  // Create output text file to write out energies
  FILE *fptr;
  fptr = fopen("pos.txt", "w");

  // Define the simulation cell
  double L = cbrt(108.0 / 0.8);
  std::cout << "FCC cell l = " << L / 3.0 << '\n';
  Vector3D SimBox = {L, L, L};

  // Define arguments to simulation class
  double h =
      0.008311612f; // Time step size = 1% of the period of oscillation for the
                    // harmonic oscillator in the lennard jones potential
  SimArgs args = {.BoxLength = L,
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
    Vector3D tmp = sim.GetPos(j);
    fprintf(fptr, "%f %f %f ", tmp[0], tmp[1], tmp[2]);
  }
  fprintf(fptr, "\n");
  for (int i = 0; i < 6'000; ++i) {
    sim.VelocityVerletStep();
    for (int j = 0; j < 108; ++j) {
      Vector3D tmp = sim.GetPos(j);
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