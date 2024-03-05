#include "../include/Sim.h"
#include "../include/Ar/constants.h"
#include <cstdio>

Vector3D gFCC[108];
Vector3D gV0[108];
Vector3D FCC(int i) { return gFCC[i]; }
Vector3D V0(int i) { return gV0[i]; }
void InitFCC(double l) {
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

void Sim::InitSim(SimArgs args) {
  fBoxLength = args.BoxLength;
  h = args.TimeStepSize;
  U_ij = args.PotentialEnergyFunction;
  F_ij = args.ForceFunction;
  fSimBox = Vector3D(fBoxLength, fBoxLength, fBoxLength);
  fSimBoxHalf = 0.5 * fSimBox;
  InitPositions();
  InitVelocities();
}

// Initialize position of atoms to FCC crystal lattice
void Sim::InitPositions() { InitFCC(fBoxLength / 3.0); }

void Sim::InitVelocities() {
  InitV0();
  fParticleHead = new Particle(FCC(0), V0(0), 0);
  Particle *NewParticle = new Particle(FCC(1), V0(1), 1);
  fParticleHead->next = NewParticle;
  for (int i = 2; i < fNumAtoms; ++i) {
    NewParticle = (NewParticle->next = new Particle(FCC(i), V0(i), i));
    NewParticle->next = nullptr;
  }
}

Vector3D Sim::GetNetMomentum() const {
  Vector3D result = {0.0, 0.0, 0.0};
  // Loop over all particles in list
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v;
  }
  return result;
}

double Sim::GetTemperature() const {
  double result = 0.0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v.sq();
  }
  result /= 3.0 * (fNumAtoms - 1.0);
  return result;
}

void Sim::ApplyPeriodicBC() {
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    if (IsInSimBox(entry->r))
      continue;
    entry->r -= (copysign_vector3(fSimBoxHalf, entry->r - fSimBox) +
                 copysign_vector3(fSimBoxHalf, entry->r));
  }
}

inline Vector3D Sim::copysign_vector3(const Vector3D &_lhs,
                                      const Vector3D &_rhs) {
  return Vector3D(copysign(_lhs[0], _rhs[0]), copysign(_lhs[1], _rhs[1]),
                  copysign(_lhs[2], _rhs[2]));
}

bool Sim::IsInSimBox(const Vector3D &check) {
  for (int i = 0; i < 3; ++i) {
    if (check[i] > fSimBox[i] || check[i] < 0)
      return false;
  }
  return true;
}

inline bool Sim::ImageIsCloser(const Vector3D &rij) {
  return (fabs(rij[0]) > 0.5 * fBoxLength || fabs(rij[1]) > 0.5 * fBoxLength ||
          fabs(rij[2]) > 0.5 * fBoxLength);
}

// this function gets called a lot, so it's worth optimizing
const Vector3D &Sim::MinimumImage(const Vector3D &ri, const Vector3D &rj) {
  rij[0] = ri[0] - rj[0];
  rij[1] = ri[1] - rj[1];
  rij[2] = ri[2] - rj[2];
  if (ImageIsCloser(rij)) {
    rij[0] -= copysign(fSimBoxHalf[0], rij[0] - fSimBoxHalf[0]) +
              copysign(fSimBoxHalf[0], rij[0] + fSimBoxHalf[0]);
    rij[1] -= copysign(fSimBoxHalf[1], rij[1] - fSimBoxHalf[1]) +
              copysign(fSimBoxHalf[1], rij[1] + fSimBoxHalf[1]);
    rij[2] -= copysign(fSimBoxHalf[2], rij[2] - fSimBoxHalf[2]) +
              copysign(fSimBoxHalf[2], rij[2] + fSimBoxHalf[2]);
  }
  return rij;
}

void Sim::ComputeForces() {
  // zero them out
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    entry->a.assign(0.0);
  }
  for (entry = fParticleHead; entry; entry = entry->next) {
    Particle *entry2 = entry->next;
    while (entry2) {
      _force.assign_to(F_ij(MinimumImage(entry->r, entry2->r)));
      entry->a += _force;
      entry2->a -= _force;
      entry2 = entry2->next;
    }
  }
}

void Sim::ComputeForcesNoPBC() {
  // std::cout << "inside computeforcesnoPBC\n";
  Particle *entry;
  for (entry = fParticleHeadNoPBC; entry; entry = entry->next) {
    entry->a.assign(0.0);
  }
  for (entry = fParticleHeadNoPBC; entry; entry = entry->next) {
    Particle *entry2 = entry->next;
    while (entry2) {
      _force.assign_to(F_ij(entry->r - entry2->r));
      /*std::cout << "r1 = " << entry->r << " r2 = " << entry2->r << "\n";
      std::cout << "\t force = " << _force << "\n";
      std::cin.get();*/
      entry->a += _force;
      entry2->a -= _force;
      entry2 = entry2->next;
    }
  }
}

void Sim::VelocityVerletStep() {
  ComputeForces();
  Particle *entry = fParticleHead;
  while (entry) {
    entry->r += h * entry->v + 0.5 * h * h * entry->a;
    entry->v += 0.5 * h * entry->a;
    entry = entry->next;
  }
  ApplyPeriodicBC();
  ComputeForces();
  entry = fParticleHead;
  while (entry) {
    entry->v += 0.5 * h * entry->a;
    entry = entry->next;
  }
}

void Sim::VelocityVerletStepNoPBC() {
  // std::cout << "inside VelocityVerletStepNoPBC" << std::endl;
  ComputeForcesNoPBC();
  // std::cout << "computed forces\n";
  Particle *entry = fParticleHeadNoPBC;
  while (entry) {
    /*std::cout << "entry->r = " << entry->r << " entry->v = " << entry->v
              << " entry->a = " << entry->a << std::endl;
    std::cin.get();*/
    entry->r += h * entry->v + 0.5 * h * h * entry->a;
    entry->v += 0.5 * h * entry->a;
    entry = entry->next;
  }
  ComputeForcesNoPBC();
  entry = fParticleHeadNoPBC;
  while (entry) {
    entry->v += 0.5 * h * entry->a;
    entry = entry->next;
  }
}

double Sim::GetPotentialEnergy() {
  double result = 0.0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    Particle *entry2;
    for (entry2 = entry->next; entry2; entry2 = entry2->next) {
      result += U_ij(MinimumImage(entry->r, entry2->r));
    }
  }
  return result;
}

double Sim::GetKineticEnergy() {
  double result = 0.0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v.sq();
  }
  return 0.5 * result;
}

double Sim::GetInternalEnergy() {
  return GetKineticEnergy() + GetPotentialEnergy();
}

void Sim::CorrectTemperature() {
  double T_corr = 1.0 / sqrt(GetTemperature());
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    entry->v *= T_corr;
  }
}

Vector3D Sim::GetMeanVelocity() const {
  Vector3D result = {0.0, 0.0, 0.0};
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    result += entry->v;
  }
  result /= 108.0;
  return result;
}

double Sim::Getfvelo(int j) const {
  double numerator = 0.0;
  double denominator = 0.0;
  Vector3D mean = GetMeanVelocity();
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    double tmpn = entry->v[j] - mean[j];
    numerator += tmpn * tmpn * tmpn * tmpn;
    denominator += tmpn * tmpn;
  }
  numerator /= 108.0;
  denominator /= 108.0;
  return numerator / (denominator * denominator);
}

void Sim::Equilibrate(int nsteps) {
  CorrectTemperature();
  for (int i = 0; i < nsteps; ++i) {
    VelocityVerletStep();
    CorrectTemperature();
  }
}

void Sim::ComputeVirial() {
  double result = 0.0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next) {
    Particle *entry2;
    for (entry2 = entry->next; entry2; entry2 = entry2->next) {
      rij = MinimumImage(entry->r, entry2->r);
      result += 24.0 * (2.0 / pow6(rij) - 1.0) / pow6(rij);
    }
  }
  fVirial = result;
  fViralSum += result;
}

void Sim::UpdatePressure(int nsteps) {
  ComputeVirial();
  double VirialAve = fViralSum / nsteps;

  fPressure =
      (fNumAtoms * GetTemperature() + VirialAve / 3.0) / pow(fBoxLength, 3.0);
}

void Sim::PhaseI(int nsteps) {
  fViralSum = 0;
  FILE *fptr;
  fptr = fopen("CanoncialEnsembleI.txt", "w");
  fprintf(fptr, "#i E KE PE T P\n");

  CorrectTemperature();
  ComputeVirial();
  double Pressure = (108.0 + fViralSum / 3.0) / pow(fBoxLength, 3.0);
  fprintf(fptr, "%f %f %f %f %f %f\n", 0.0, GetInternalEnergy() * Ar::gEnergy,
          GetKineticEnergy() * Ar::gEnergy, GetPotentialEnergy() * Ar::gEnergy,
          GetTemperature() * Ar::gTemperatureArgon, Pressure * Ar::gPressure);

  // Equilibrate
  for (int i = 0; i < nsteps; ++i) {
    VelocityVerletStep();
    CorrectTemperature();
    ComputeVirial();
    Pressure = (108.0 * GetTemperature() + fViralSum / (3.0 * (i + 2.0))) /
               pow(fBoxLength, 3.0);
    fprintf(fptr, "%f %f %f %f %f %f\n", (i + 1) * h * Ar::gTime,
            GetInternalEnergy() * Ar::gEnergy, GetKineticEnergy() * Ar::gEnergy,
            GetPotentialEnergy() * Ar::gEnergy,
            GetTemperature() * Ar::gTemperatureArgon, Pressure * Ar::gPressure);
  }
  fclose(fptr);
}

void Sim::InitR0() {
  int i = 0;
  Particle *entry;
  for (entry = fParticleHead; entry; entry = entry->next, i++) {
    fR0[i] = entry->r;
  }

  // init fParticleHeadNoPBC
  Particle *NewParticle = new Particle(
      fParticleHead->next->r, fParticleHead->next->v, fParticleHead->next->ID);
  fParticleHeadNoPBC =
      new Particle(fParticleHead->r, fParticleHead->v, fParticleHead->ID);
  fParticleHeadNoPBC->next = NewParticle;
  for (entry = fParticleHead->next->next; entry; entry = entry->next) {
    NewParticle =
        (NewParticle->next = new Particle(entry->r, entry->v, entry->ID));
    NewParticle->next = nullptr;
  }
}

void Sim::PhaseII(int nsteps) {
  // SETUP:
  fViralSum = 0.0;
  InitR0();

  // output file
  FILE *fptr;
  fptr = fopen("CanoncialEnsembleII.txt", "w");
  fprintf(fptr, "#i E KE PE T dT P Cv R2\n");

  double Cv = 0.5 * 3 * fNumAtoms;
  double KEsum, KEsum2, KEave, KEave2;
  double Tsum, Tsum2, Tave, Tave2;
  KEsum = GetKineticEnergy();
  KEsum2 = KEsum * KEsum;
  Tsum = GetTemperature();
  Tsum2 = Tsum * Tsum;
  ComputeVirial();
  double Pave =
      (fNumAtoms * GetTemperature() + fVirial / 3.0) / pow(fBoxLength, 3.0);
  double R2;

  // output initial values
  fprintf(fptr, "%d %f %f %f %f %f %f %f %f\n", 0,
          Ar::gEnergy * GetInternalEnergy() / 108.0,
          Ar::gEnergy * GetKineticEnergy() / 108.0,
          Ar::gEnergy * GetPotentialEnergy() / 108.0,
          Ar::gTemperatureArgon * Tsum, 0., Pave * Ar::gPressure, Cv, 0.0);

  // Compute Time Averages
  for (int i = 0; i < nsteps; ++i) {
    VelocityVerletStep();
    KEsum += GetKineticEnergy();
    KEsum2 += GetKineticEnergy() * GetKineticEnergy();
    KEave = KEsum / (i + 2.0);
    KEave2 = KEsum2 / (i + 2.0);
    Tsum += GetTemperature();
    Tsum2 += GetTemperature() * GetTemperature();
    Tave = Tsum / (i + 2.0);
    Tave2 = Tsum2 / (i + 2.0);
    Cv = 1.0 /
         (2.0 / (3.0 * fNumAtoms) - (KEave2 - KEave * KEave) / (KEave * KEave));
    ComputeVirial();
    VelocityVerletStepNoPBC();
    R2 = GetMeanSquareDisplacement();
    std::cout << R2 / (6.0 * (i + 1.0) * h * Ar::gTime) << std::endl;
    Pave = ((fNumAtoms * Tave + fViralSum / (3.0 * (i + 2.0))) /
            pow(fBoxLength, 3.0));
    fprintf(fptr, "%f %f %f %f %f %f %f %f %f\n", (i + 1) * h * Ar::gTime,
            Ar::gEnergy * GetInternalEnergy() / 108.0,
            Ar::gEnergy * GetKineticEnergy() / 108.0,
            Ar::gEnergy * GetPotentialEnergy() / 108.0,
            Ar::gTemperatureArgon * Tave, (Tave2 - Tave * Tave) / (Tave * Tave),
            Pave * Ar::gPressure, Cv, R2 * Ar::gLength * Ar::gLength);
  }
  fclose(fptr);
}

double Sim::GetMeanSquareDisplacement() {
  double result = 0.0;
  Particle *entry;
  for (entry = fParticleHeadNoPBC; entry; entry = entry->next) {
    result += (entry->r - fR0[entry->ID]).sq();
  }
  return result / fNumAtoms;
}