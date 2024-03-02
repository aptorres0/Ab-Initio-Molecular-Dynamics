#include "../include/Sim.h"

void Sim::InitSim(SimArgs args) {
  fSimBox = args.SimulationCellBounds;
  fSimBoxHalf = 0.5f * fSimBox;
  LHalf = fSimBoxHalf[0];
  h = args.TimeStepSize;
  U_ij = args.PotentialEnergyFunction;
  F_ij = args.ForceFunction;
}

// Initialize position of atoms to FCC crystal lattice
void Sim::InitPositions() {
  float l = fSimBox[0] / 3.0;
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    float xmin = i * fSimBox[0] / 3.0;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      float ymin = j * fSimBox[1] / 3.0;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        float zmin = k * fSimBox[2] / 3.0;
        fEnsemble.r(i * 36 + j * 12 + k * 4) = {xmin, ymin, zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 1) = {xmin, 0.5f * l + ymin,
                                                    0.5f * l + zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 2) = {0.5f * l + xmin, ymin,
                                                    0.5f * l + zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 3) = {0.5f * l + xmin,
                                                    0.5f * l + ymin, zmin};
      }
    }
  }
}

void Sim::InitVelocities() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless
   * variance of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // 1. set velocities from gaussian distribution
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i).set_random();
  }

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  Vector3F NetMomentum = GetNetMomentum();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) -= NetMomentum / ((float)fNumAtoms);
  }

  // 3. set temperature of the system by multiplying by a correction factor
  /*float T_corr = 1.0 / sqrt(GetTemperature());
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) *= T_corr;
  }*/
  CorrectTemperature();
}

Vector3F Sim::GetNetMomentum() const {
  Vector3F result = {0.0, 0.0, 0.0};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i);
  }
  return result;
}

float Sim::GetTemperature() const {
  float result = 0.0;
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i).sq();
  }
  result /= 3.0 * (fNumAtoms - 1.0);
  return result;
}

void Sim::ApplyPeriodicBC() {
  for (int i = 0; i < fNumAtoms; ++i) {
    if (IsInSimBox(fEnsemble.r(i)))
      continue;
    fEnsemble.r(i) -= (copysign_vector3(fSimBoxHalf, fEnsemble.r(i) - fSimBox) +
                       copysign_vector3(fSimBoxHalf, fEnsemble.r(i)));
  }
}

inline Vector3F Sim::copysign_vector3(const Vector3F &_lhs,
                                      const Vector3F &_rhs) {
  /*Vector3F result{0.0, 0.0, 0.0};
  result[0] = copysign(_lhs[0], _rhs[0]);
  result[1] = copysign(_lhs[1], _rhs[1]);
  result[2] = copysign(_lhs[2], _rhs[2]);
  return result;*/
  return Vector3F(copysign(_lhs[0], _rhs[0]), copysign(_lhs[1], _rhs[1]),
                  copysign(_lhs[2], _rhs[2]));
}

bool Sim::IsInSimBox(const Vector3F &check) {
  for (int i = 0; i < 3; ++i) {
    if (check[i] > fSimBox[i] || check[i] < 0)
      return false;
  }
  return true;
}

inline bool Sim::ImageIsCloser(const Vector3F &rij) {
  return (fabs(rij[0]) > LHalf || fabs(rij[1]) > LHalf || fabs(rij[2]) > LHalf);
}

// this function gets called a lot, so it's worth optimizing
const Vector3F &Sim::MinimumImage(const Vector3F &ri, const Vector3F &rj) {
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
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.a(i).assign(0.0); // avoid a copy
  }
  for (int i = 0; i < fNumAtoms - 1; ++i) {
    for (int j = i + 1; j < fNumAtoms; ++j) {
      _force.assign_to(
          F_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)))); // avoid a copy
      fEnsemble.a(i) += _force;
      fEnsemble.a(j) -= _force;
    }
  }
}

void Sim::VelocityVerletStep() {
  ComputeForces();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.r(i) += h * fEnsemble.v(i) + 0.5f * h * h * fEnsemble.a(i);
    fEnsemble.v(i) += 0.5f * h * fEnsemble.a(i);
  }
  ApplyPeriodicBC();
  ComputeForces();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) += 0.5f * h * fEnsemble.a(i);
  }
}

void Sim::PrintSample() {
  std::cout << "r[0] = " << fEnsemble.r(0) << "v[0] = " << fEnsemble.v(0)
            << "a[0] = " << fEnsemble.a(0) << '\n';
  std::cout << "r[1] = " << fEnsemble.r(1) << "v[1] = " << fEnsemble.v(1)
            << "a[1] = " << fEnsemble.a(1) << '\n';
}

float Sim::GetPotentialEnergy() {
  float result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    for (int j = i + 1; j < fNumAtoms; ++j) {
      result += U_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
    }
  }
  return result;
}

float Sim::GetKineticEnergy() {
  float result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i).sq();
  }
  return result / 2.0;
}

float Sim::GetInternalEnergy() {
  return GetKineticEnergy() + GetPotentialEnergy();
}

void Sim::CorrectTemperature() {
  float T_corr = 1.0 / sqrt(GetTemperature());
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) *= T_corr;
  }
}

Vector3F Sim::GetMeanVelocity() const {
  Vector3F result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i);
  }
  result /= 108.0f;
  return result;
}

float Sim::Getfvelo(int j) const {
  float numerator = 0.0;
  float denominator = 0.0;
  Vector3F mean = GetMeanVelocity();
  for (int i = 0; i < fNumAtoms; ++i) {
    float tmpn = fEnsemble.v(i)[j] - mean[j];
    numerator += tmpn * tmpn * tmpn * tmpn;
    denominator += tmpn * tmpn;
  }
  numerator /= 108.0f;
  denominator /= 108.0f;
  return numerator / (denominator * denominator);
}