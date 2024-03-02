#include "../include/Ar/Atom.h"
#include "../include/Ar/constants.h"
#include "../include/Sim.h"
#include "../include/Vector3.h"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include <cstdio>

Vector3D gFCC[108];
Vector3D gV0[108];
Vector3D FCC(int i) { return gFCC[i]; }
Vector3D V0(int i) { return gV0[i]; }
void GenerateFCC();
void GenerateV0();

void Step7();
void MicrocanonicalEnsemble();
void CanonicalEnsemble();

int main() {
  // Step7();
  MicrocanonicalEnsemble();
  CanonicalEnsemble();
  return 0;
}

void Step7() {
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/Step7.txt", "w");
  fprintf(fptr, "#t E KE PE\n");

  GenerateFCC();
  GenerateV0();
  NewSim sim;
  sim.InitSim(gFCC, gV0);
  double PE = sim.GetPotentialEnergy();
  double KE = sim.GetKineticEnergy();
  double dt = sim.GetTimeStepSize();
  Vector3D VelocityDistribution = sim.GetVelocityDistribution();
  TFile f("/Users/ajpt/School/Project3/data/Step7.root", "RECREATE");
  TH1D h_V("h_V", "Velocity Distribution;Velocity [r.u.]", 30, 2, 5);
  h_V.Fill(VelocityDistribution.x());
  h_V.Fill(VelocityDistribution.y());
  h_V.Fill(VelocityDistribution.z());
  std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
  /*fprintf(fptr, "%d %f %f %f\n", 0, (KE + PE), KE, PE);*/
  fprintf(fptr, "%d %f %f %f\n", 0, (KE + PE) * Ar::gEnergy, KE * Ar::gEnergy,
          PE * Ar::gEnergy);
  for (int i = 0; i < 2'000; ++i) {
    sim.VelocityVerletStep();
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
    Vector3D VelocityDistribution = sim.GetVelocityDistribution();
    fprintf(fptr, "%f %f %f %f\n", (i + 1) * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy);
    /*fprintf(fptr, "%f %f %f %f\n", (i + 1.0), (KE + PE), KE, PE);*/
    VelocityDistribution = sim.GetVelocityDistribution();
    h_V.Fill(VelocityDistribution.x());
    h_V.Fill(VelocityDistribution.y());
    h_V.Fill(VelocityDistribution.z());
    std::cout << "fv = " << VelocityDistribution * Ar::gLength / Ar::gTime
              << '\n';
  }
  std::cout << "\n\tFinal Net Momentum: " << sim.GetNetMomentum() << '\n';
  h_V.Write();
  f.Close();
}

void MicrocanonicalEnsemble() {
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/Microcanonical.txt", "w");
  fprintf(fptr, "#t[ps] E KE PE[eV]\n");
  TFile f("ThermalWaveLength.root", "RECREATE");
  TH1D h_Newton("h_Newton",
                "Mean Separation Distance over Thermal de Broglie Wavelength",
                100, 0, 10);
  double a = cbrt(3.0 / (4.0 * M_PI * Ar::gNumberDensity_r));

  GenerateFCC();
  GenerateV0();
  NewSim sim;
  sim.InitSim(gFCC, gV0);
  double KE, PE;
  KE = sim.GetKineticEnergy();
  PE = sim.GetPotentialEnergy();
  fprintf(fptr, "%d %f %f %f\n", 0, (KE + PE) * Ar::gEnergy, KE * Ar::gEnergy,
          PE * Ar::gEnergy);
  double dt = sim.GetTimeStepSize();
  for (int i = 0; i < 60'000; ++i) {
    sim.VelocityVerletStep();
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    fprintf(fptr, "%f %f %f %f\n", (i + 1) * dt, (KE + PE) * Ar::gEnergy,
            KE * Ar::gEnergy, PE * Ar::gEnergy);
    h_Newton.Fill(a / sim.GetThermalWavelength());
    std::cout << KE + PE << " " << KE << " " << PE << " " << '\n';
  }
  h_Newton.Write();
  f.Close();
  fclose(fptr);
}

void CanonicalEnsemble() {
  FILE *fptr;
  fptr = fopen("/Users/ajpt/School/Project3/data/CanonicalI.txt", "w");
  fprintf(fptr, "#t[ps] E KE PE[eV] T[K] P[MPa]\n");
  double L = cbrt(108.0 / 0.8); // length of simulation cell [r.u.]

  GenerateFCC();
  GenerateV0();
  NewSim sim;
  sim.InitSim(gFCC, gV0);
  double KE, PE, Temp, Pressure;
  double dt = sim.GetTimeStepSize();
  sim.ResetPressure();
  for (int i = 0; i < 60'000; ++i) {
    sim.RescaleVelocities();
    PE = sim.GetPotentialEnergy();
    KE = sim.GetKineticEnergy();
    Temp = sim.GetTemperature();
    Pressure = sim.GetPressure();
    std::cout << KE + PE << " " << KE << " " << PE << " " << Temp << " "
              << Pressure << '\n';
    fprintf(fptr, "%f %f %f %f %f %f\n", (i + 1) * dt, (KE + PE) * Ar::gEnergy,
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
  FILE *fptr2;
  fptr2 = fopen("/Users/ajpt/School/Project3/data/CanonicalII.txt", "w");
  fprintf(fptr2,
          "#t[ps] <E> <KE> <PE>[eV] <T> dT[K] <P>[MPa] Cv[J/kg/K] MSD[nm^2]\n");
  // 0a. Reset all of the sums used to determine the averages in phase I
  double Esum = PE + KE, PEsum = PE, KEsum = KE, KE2sum = KE * KE,
         Tempsum = Temp, Temp2sum = Temp * Temp;
  int CountSum = 1;
  double Eave, PEave, KEave, KE2ave, Tempave, Temp2ave, TempFluc, Cv,
      Pressureave;
  sim.ResetPressure();
  // 0b. Setup for Mean-Square Displacement (MSD) Calculations
  sim.InitMSD();
  double MSD;
  // 0c. Setup radial distribution function histogram
  sim.BookRadialDistributionHisto();
  // 1. Run the simulation for 2000 steps
  for (int i = 0; i < 20'000; ++i) {
    sim.VelocityVerletStep();
    sim.StepMSD();
    sim.FillRadialDistributionHisto();
    CountSum += 1;
    PE = sim.GetPotentialEnergy();
    PEsum += PE;
    PEave = PEsum / CountSum;
    KE = sim.GetKineticEnergy();
    KEsum += KE;
    KE2sum += KE * KE;
    KEave = KEsum / CountSum;
    KE2ave = KE2sum / CountSum;
    Esum += PE + KE;
    Eave = Esum / CountSum;
    Temp = sim.GetTemperature();
    Tempsum += Temp;
    Temp2sum += Temp * Temp;
    Tempave = Tempsum / CountSum;
    Temp2ave = Temp2sum / CountSum;
    TempFluc = (Temp2ave - (Tempave * Tempave)) / (Tempave * Tempave);
    Pressureave =
        (108.0 * Tempave + sim.GetVirialAverage() / 3.0) / (L * L * L);
    Cv = 1.0 /
         ((2.0 / (3 * 108.0)) - ((KE2ave - (KEave * KEave)) / (KEave * KEave)));
    MSD = sim.GetMeanSquaredDisplacement();
    std::cout << Eave / 108.0 << " " << KEave / 108.0 << " " << PEave / 108.0
              << " " << Tempave << " " << TempFluc << " " << Pressureave << " "
              << Cv << '\n';
    fprintf(fptr2, "%f %f %f %f %f %f %f %f %f\n", (i + 1) * dt,
            Eave * Ar::gEnergy / 108.0, KEave * Ar::gEnergy / 108.0,
            PEave * Ar::gEnergy / 108.0, Tempave * Ar::gTemperature,
            TempFluc * Ar::gTemperature, Pressureave * Ar::gPressure,
            Cv * Ar::gSpecificHeatCapacity, MSD * Ar::gLength2);
  }
  Histo RadDist = sim.GetRadialDistributionHisto();
  fclose(fptr2);
  FILE *fptr3 = fopen("/Users/ajpt/School/Project3/data/radialdist.txt", "w");
  fprintf(fptr3, "#r[nm] g(r)\n");
  for (int i = 0; i < RadDist.nBins; ++i) {
    fprintf(fptr3, "%f %f\n", (RadDist.binWidth * i) * Ar::gLength,
            RadDist.bin[i]);
  }
  fclose(fptr3);
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