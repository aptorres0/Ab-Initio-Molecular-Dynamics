#include "../include/LennardJones.h"
#include "../include/Particle.h"
#include "../include/Sim.h"
#include "../include/Vector3.h"
#include <cmath>
#include <cstdio>

void SimpleSim();
void GenerateAnimationData();

Vector3F gFCC[108];
Vector3F gV0[108];
void InitFCC(float l) {
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    float xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      float ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        float zmin = k * l;
        gFCC[i * 36 + j * 12 + k * 4] = {xmin, ymin, zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 1] = {xmin, 0.5f * l + ymin,
                                             0.5f * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 2] = {0.5f * l + xmin, ymin,
                                             0.5f * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 3] = {0.5f * l + xmin, 0.5f * l + ymin,
                                             zmin};
      }
    }
  }
}
void InitV0() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless
   * variance of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // 1. set velocities from gaussian distribution
  Vector3F NetMomentum = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < 108; ++i) {
    gV0[i].set_random();
    NetMomentum += gV0[i];
  }

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  float Temp = 0.0f;
  for (int i = 0; i < 108; ++i) {
    gV0[i] -= (NetMomentum / 108.0f);
    Temp += gV0[i].sq();
  }
  Temp /= (3.0 * 107.0);

  // 3. set temperature of the system by multiplying by a correction factor
  float T_corr = 1.0 / sqrt(Temp);
  for (int i = 0; i < 108; ++i) {
    gV0[i] *= T_corr;
  }
}
Vector3F FCC(int i) { return gFCC[i]; }
Vector3F V0(int i) { return gV0[i]; }

void test() {
  ParticleListt *list = new ParticleListt();

  float L = cbrt(108.0 / 0.8) / 3.0f;

  InitFCC(L);
  InitV0();
  std::cout << "Starting to add particles to list...\n";
  for (int i = 0; i < 108; ++i) {
    list->AddParticle(new Particle(FCC(i), V0(i)));
    std::cout << FCC(i) << "\n";
  }
  std::cout << "Finished adding particles to list...\n";

  // compute temperature
  float result = 0.0;
  for (ParticleListt::Iterator it = list->begin(); it != list->end(); ++it) {
    std::cout << (*it).r << "\n";
    result += (*it).v.sq();
  }
  result /= 3.0 * (108 - 1.0);
  std::cout << "Initial Temperature [reduced units] = " << result << '\n';
}

int main() {
  // test();
  SimpleSim();
  //  GenerateAnimationData();
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
  for (int i = 0; i < 6'000; ++i) {
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