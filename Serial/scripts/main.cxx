/**
 * main.cxx
 *
 * Author: Alexander Paul Torres
 * Date: 14 APR 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Main file for Project 3. Contains the main function which implements the
 * `Sim` class and functions to generate the initial conditions for the
 * simulation.
 */
#include "../include/Ar/constants.h"
#include "../include/Sim.h"
#include "../include/Vector3.h"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

Vector3D gFCC[108];
Vector3D gV0[108];
void GenerateFCC();
void GenerateV0();

void Step7();
void MicrocanonicalEnsemble();
void CanonicalEnsemble();
void GenerateAnimationData();

int main() {
  // Step7();
  // MicrocanonicalEnsemble();
  CanonicalEnsemble();
  // GenerateAnimationData();

  return 0;
}

void Step7() {
  // Declare output file for energies
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/Step7_new.txt", "w");
  fprintf(fptr, "#t E KE PE\n");

  // Declare output file for ROOT histogram of velocity distribution
  TFile f("/Users/ajpt/School/Project3/data/Step7_new.root", "RECREATE");
  TH1D h_V("h_V", "Velocity Distribution;Velocity [r.u.]", 30, 2, 5);

  // Initialize simulation
  Sim sim;
  GenerateFCC();
  GenerateV0();
  double dt = 0.01;
  sim.InitSim(gFCC, gV0, dt);

  // Get initial energies and write to output file
  double KE = sim.GetKineticEnergy();
  double PE = sim.GetPotentialEnergy();
  std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
  fprintf(fptr, "%d %f %f %f\n", 0, (KE + PE) * Ar::gEnergy, KE * Ar::gEnergy,
          PE * Ar::gEnergy);

  // Get Velocity distribution and fill histogram
  Vector3D VelocityDistribution = sim.GetVelocityDistribution();
  h_V.Fill(VelocityDistribution.x());
  h_V.Fill(VelocityDistribution.y());
  h_V.Fill(VelocityDistribution.z());

  // Run simulation for 2'000 steps with a time step of 0.01 [r.u.]
  for (int i = 0; i < 2'000; ++i) {
    // Update atom positions, velocities, and accelerations
    sim.VelocityVerletStep();
    // Get energies and write to output file
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
    fprintf(fptr, "%f %f %f %f\n", (i + 1) * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy);
    // Get Velocity distribution and fill histogram
    VelocityDistribution = sim.GetVelocityDistribution();
    std::cout << "fv = " << VelocityDistribution * Ar::gVelocity << '\n';
    h_V.Fill(VelocityDistribution.x());
    h_V.Fill(VelocityDistribution.y());
    h_V.Fill(VelocityDistribution.z());
  }
  std::cout << "\n\tFinal Net Momentum: " << sim.GetNetMomentum() << '\n';
  fclose(fptr);

  // Write histogram to ROOT file
  h_V.Write();
  f.Close();
}

void MicrocanonicalEnsemble() {
  // Declare output file for energies
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/Microcanonical_new.txt", "w");
  fprintf(fptr, "#t[ps] E KE PE[eV]\n");
  std::cout << std::setprecision(3);

  // Initialize simulation
  Sim sim;
  GenerateFCC();
  GenerateV0();
  // double dt = 0.0008311612;
  double dt = 2e-5;
  sim.InitSim(gFCC, gV0, dt);

  // Get initial energies and write to output file
  double KE, PE;
  KE = sim.GetKineticEnergy();
  PE = sim.GetPotentialEnergy();
  fprintf(fptr, "%d %f %f %f\n", 0, (KE + PE) * Ar::gEnergy, KE * Ar::gEnergy,
          PE * Ar::gEnergy);

  // Run simulation for 60'000 steps
  for (int i = 0; i < 6'000'000; ++i) {
    // Update atom positions, velocities, and accelerations
    sim.VelocityVerletStep();
    // Get energies and write to output file
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
    fprintf(fptr, "%f %f %f %f\n", (i + 1) * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy);
  }
  fclose(fptr);
}

void CanonicalEnsemble() {
  // Declare output file for energies during phase I
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/CanonicalI_new.txt", "w");
  fprintf(fptr, "#t[ps] E KE PE[eV] T[K] P[MPa]\n");
  // Declare output file for energies during phase I and II
  FILE *fptrAll;
  fptrAll = fopen("/Users/ajpt/School/Project3/data/Canonical_All.txt", "w");
  fprintf(fptrAll, "#t[ps] E KE PE[eV] T[K] P[MPa]\n");
  // Declare output file to animate the system
  FILE *fptr_Anim;
  fptr_Anim = fopen("/Users/ajpt/School/Project3/data/pos_Canonical.txt", "w");
  std::cout << std::setprecision(3);

  // Initialize simulation
  Sim sim;
  GenerateFCC();
  GenerateV0();
  // double dt = 0.0008311612;
  double dt = 2e-5;
  sim.InitSim(gFCC, gV0, dt);

  // Declare variables for storing system properties and energy
  double KE, PE, Temp, Pressure;
  sim.ResetPressure(); // set virial to zero

  // *******************************************************
  // Start Phase 1 - Equilibration
  // *******************************************************

  // Run Equilibration for 60'000 steps
  for (int i = 0; i < 60'000; ++i) {
    // Update atom r, v, and a
    sim.VelocityVerletStep();
    // Maintain constant temperature
    sim.BerendsenThermostat(0.1);
    // sim.GaussianThermostat();
    //  Write out positions at each step for animation
    for (int j = 0; j < 108; ++j) {
      Vector3D pos = sim.GetPositionOfAtom(j);
      fprintf(fptr_Anim, "%f %f %f ", pos[0], pos[1], pos[2]);
    }
    fprintf(fptr_Anim, "\n");
    // Get energies, P, and T and write to output files
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    Temp = sim.GetTemperature();
    Pressure = sim.GetPressure();
    std::cout << std::setprecision(3);
    std::cout << i << " " << std::setw(8) << KE + PE << " " << std::setw(8)
              << KE << " " << std::setw(8) << PE << " " << Temp << " "
              << std::setw(8) << Pressure << '\n';
    fprintf(fptr, "%f %f %f %f %f %f\n", i * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy, Temp * Ar::gTemperature,
            Pressure * Ar::gPressure);
    fprintf(fptrAll, "%f %f %f %f %f %f\n", i * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy, Temp * Ar::gTemperature,
            Pressure * Ar::gPressure);
  }
  std::cout << "\n\tNet Momentum after equilibration: " << sim.GetNetMomentum()
            << "\n\n";
  fclose(fptr);
  // std::cin.get();

  // *******************************************************
  // Start Phase 2
  // *******************************************************
  // Declare output file for energies during phase II
  FILE *fptr2;
  fptr2 = fopen("/Users/ajpt/School/Project3/data/CanonicalII_new.txt", "w");
  fprintf(
      fptr2,
      "#t[ps] <E> <KE> <PE>[eV] <T>[K] dT[1] <P>[MPa] Cv[J/kg/K] MSD[nm^2]\n");

  // Declare output file for radial distribution function
  FILE *fptr3 =
      fopen("/Users/ajpt/School/Project3/data/radialdist_new.txt", "w");
  fprintf(fptr3, "#r[nm] g(r)\n");

  // Reset all of the sums used to determine the averages in phase I
  sim.ResetPressure();
  double Esum = 0, PEsum = 0, KEsum = 0, KE2sum = 0, Tempsum = 0, Temp2sum = 0;
  int CountSum = 0;
  double Eave, PEave, KEave, KE2ave, Tempave, Temp2ave, TempFluc, Cv,
      Pressureave;
  double L = cbrt(108.0 / 0.8);

  // Setup for Mean-Square Displacement (MSD) Calculations
  sim.InitMSD();
  double MSD;

  // Run simulation for 20'000 steps
  for (int i = 0; i < 20'000; ++i) {
    // Update atom positions, velocities, and accelerations
    sim.VelocityVerletStep();
    // Write out positions at each step for animation
    for (int j = 0; j < 108; ++j) {
      Vector3D pos = sim.GetPositionOfAtom(j);
      fprintf(fptr_Anim, "%f %f %f ", pos[0], pos[1], pos[2]);
    }
    fprintf(fptr_Anim, "\n");
    // Update values for MSD calculations (No PBC)
    sim.StepMSD();
    // Get values at each step
    KE = sim.GetKineticEnergy();
    PE = sim.GetPotentialEnergy();
    Temp = sim.GetTemperature();
    MSD = sim.GetMeanSquaredDisplacement();
    sim.FillRadialDistributionHisto();
    // Calculate sums for averages
    Esum += PE + KE;
    KEsum += KE;
    KE2sum += KE * KE;
    PEsum += PE;
    Tempsum += Temp;
    Temp2sum += Temp * Temp;
    CountSum += 1;
    // Calculate averages and write to output file
    Eave = Esum / CountSum;
    KEave = KEsum / CountSum;
    KE2ave = KE2sum / CountSum;
    PEave = PEsum / CountSum;
    Tempave = Tempsum / CountSum;
    Temp2ave = Temp2sum / CountSum;
    TempFluc = (Temp2ave - (Tempave * Tempave)) / (Tempave * Tempave);
    Pressureave = (108.0 * Tempave + sim.GetAveVirial() / 3.0) / (L * L * L);
    Cv = 1.0 /
         ((2.0 / (3 * 108.0)) - ((KE2ave - (KEave * KEave)) / (KEave * KEave)));
    std::cout << i << " " << Eave / 108.0 << " " << KEave / 108.0 << " "
              << PEave / 108.0 << " " << Tempave << " " << TempFluc << " "
              << Pressureave << " " << Cv << '\n';
    fprintf(fptr2, "%f %f %f %f %f %f %f %f %f\n", (i + 1) * dt,
            Eave * Ar::gEnergy / 108.0, KEave * Ar::gEnergy / 108.0,
            PEave * Ar::gEnergy / 108.0, Tempave * Ar::gTemperature, TempFluc,
            Pressureave * Ar::gPressure, Cv * Ar::gSpecificHeatCapacity,
            MSD * Ar::gLength2);
    fprintf(fptrAll, "%f %f %f %f %f %f\n", (i + 60000) * dt,
            Eave * Ar::gEnergy, KEave * Ar::gEnergy, PEave * Ar::gEnergy,
            Tempave * Ar::gTemperature, Pressureave * Ar::gPressure);
  }
  Histo RadDist = sim.GetRadialDistributionHisto();
  for (int i = 0; i < RadDist.nBins; ++i) {
    fprintf(fptr3, "%f %f\n", (RadDist.binWidth * i) * Ar::gLength,
            RadDist.bin[i]);
  }
  fclose(fptr2);
  fclose(fptr3);
  fclose(fptrAll);
  fclose(fptr_Anim);
}

void GenerateAnimationData() {
  // Create output file for animation data
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/pos.txt", "w");

  // Initialize simulation
  Sim sim;
  GenerateFCC();
  GenerateV0();
  double dt = 0.0008311612;
  sim.InitSim(gFCC, gV0, dt);

  // Write out initial positions
  for (int i = 0; i < 108; ++i) {
    Vector3D pos = sim.GetPositionOfAtom(i);
    fprintf(fptr, "%f %f %f ", pos[0], pos[1], pos[2]);
  }
  fprintf(fptr, "\n");

  // Run simulation for 60'000 time steps
  for (int i = 0; i < 60'000; ++i) {
    // Update atom positions, velocities, and accelerations
    sim.VelocityVerletStep();
    // Write out positions at each step
    for (int j = 0; j < 108; ++j) {
      Vector3D pos = sim.GetPositionOfAtom(j);
      fprintf(fptr, "%f %f %f ", pos[0], pos[1], pos[2]);
    }
    fprintf(fptr, "\n");
    std::cout << std::setprecision(3);
    std::cout << std::setw(8)
              << sim.GetKineticEnergy() + sim.GetPotentialEnergy() << " "
              << std::setw(8) << sim.GetKineticEnergy() << " " << std::setw(8)
              << sim.GetPotentialEnergy() << " " << std::setw(8)
              << sim.GetVelocityDistribution() << " " << std::setw(8)
              << sim.GetTemperature() << '\n';
  }
  std::cout << "Final Momentum [reduced units] = " << sim.GetNetMomentum()
            << "\nFinal Temperature [reduced units] = " << sim.GetTemperature()
            << '\n';

  fclose(fptr);
}

void GenerateFCC() {
  double l = cbrt(108.0 / 0.8) / 3.0;
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * l;
        gFCC[i * 36 + j * 12 + k * 4] = {xmin, ymin, zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 1] = {xmin, 0.5 * l + ymin,
                                             0.5 * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 2] = {0.5 * l + xmin, ymin,
                                             0.5 * l + zmin};
        gFCC[i * 36 + j * 12 + k * 4 + 3] = {0.5 * l + xmin, 0.5 * l + ymin,
                                             zmin};
      }
    }
  }
}

void GenerateV0() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless
   * variance of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // 1. set velocities from gaussian distribution
  Vector3D NetMomentum = {0.0, 0.0, 0.0};
  for (int i = 0; i < 108; ++i) {
    gV0[i].set_random();
    NetMomentum += gV0[i];
  }

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  double Temp = 0.0;
  for (int i = 0; i < 108; ++i) {
    gV0[i] -= (NetMomentum / 108.0);
    Temp += gV0[i].sq();
  }
  Temp /= (3.0 * 107.0);

  // 3. set temperature of the system by multiplying by a correction factor
  double T_corr = 1.0 / sqrt(Temp);
  for (int i = 0; i < 108; ++i) {
    gV0[i] *= T_corr;
  }
}