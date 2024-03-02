#include "../include/Sim.h"

void Sim::InitSim(SimArgs args) {
  fDIM = args.NumberOfDimensions;
  fNumAtoms = args.NumberOfAtoms;
  fSimBox = args.SimulationCellBounds;
  h = args.TimeStepSize;
  U_ij = args.PotentialEnergyFunction;
  F_ij = args.ForceFunction;
  // initialize pos, vec, and acc vectors to 0
  /*for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.r(i) = Vector3D(0.0);
    fEnsemble.v(i) = Vector3D(0.0);
    fEnsemble.a(i) = Vector3D(0.0);
  }*/
}

void Sim::InitSim(int _DIM, int _NumAtoms, double (*_U)(const Vector3D &),
                  Vector3D (*_F)(const Vector3D &), Vector3D _SimBox,
                  double _h) {
  fDIM = _DIM;
  fNumAtoms = _NumAtoms;
  U_ij = _U;
  F_ij = _F;
  fSimBox = _SimBox;
  h = _h;

  // initialize pos, vec, and acc vectors to 0
  /*for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.r(i) = Vector3D(0.0);
    fEnsemble.v(i) = Vector3D(0.0);
    fEnsemble.a(i) = Vector3D(0.0);
  }*/
}

// Initialize position of atoms to FCC crystal lattice
void Sim::InitPositions() {
  double l = fSimBox[0] / 3.0;
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * fSimBox[0] / 3.0;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * fSimBox[1] / 3.0;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * fSimBox[2] / 3.0;
        fEnsemble.r(i * 36 + j * 12 + k * 4) = {xmin, ymin, zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 1) = {xmin, 0.5 * l + ymin,
                                                    0.5 * l + zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 2) = {0.5 * l + xmin, ymin,
                                                    0.5 * l + zmin};
        fEnsemble.r(i * 36 + j * 12 + k * 4 + 3) = {0.5 * l + xmin,
                                                    0.5 * l + ymin, zmin};
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
    fEnsemble.v(i).set_rnd();
  }

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  Vector3D NetMomentum = GetNetMomentum();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) -= NetMomentum / ((double)fNumAtoms);
  }

  // 3. set temperature of the system by multiplying by a correction factor
  double T_corr = 1.0 / sqrt(GetTemperature());
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) *= T_corr;
  }
}

Vector3D Sim::GetNetMomentum() const {
  Vector3D result = {0.0, 0.0, 0.0};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i);
  }
  return result;
}

double Sim::GetTemperature() const {
  double result = 0.0;
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
    fEnsemble.r(i) -= (copysign_v(0.5 * fSimBox, fEnsemble.r(i) - fSimBox) +
                       copysign_v(0.5 * fSimBox, fEnsemble.r(i)));
  }
}

Vector3D Sim::copysign_v(const Vector3D &_lhs, const Vector3D &_rhs) {
  Vector3D result{0.0, 0.0, 0.0};
  for (int i = 0; i < fDIM; ++i) {
    result[i] = copysign(_lhs[i], _rhs[i]);
  }
  return result;
}

bool Sim::IsInSimBox(const Vector3D &check) {
  for (int i = 0; i < fDIM; ++i) {
    if (check[i] > fSimBox[i] || check[i] < 0)
      return false;
  }
  return true;
}

const Vector3D &Sim::MinimumImage(const Vector3D &ri, const Vector3D &rj) {
  for (int i = 0; i < fDIM; ++i) {
    rij[i] = ri[i] - rj[i] -
             copysign(0.5 * fSimBox[i], ri[i] - rj[i] - 0.5 * fSimBox[i]) -
             copysign(0.5 * fSimBox[i], ri[i] - rj[i] + 0.5 * fSimBox[i]);
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
    fEnsemble.r(i) += h * fEnsemble.v(i) + 0.5 * h * h * fEnsemble.a(i);
    fEnsemble.v(i) += 0.5 * h * fEnsemble.a(i);
  }
  ApplyPeriodicBC();
  ComputeForces();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) += 0.5 * h * fEnsemble.a(i);
  }
}

void Sim::PrintSample() {
  std::cout << "r[0] = " << fEnsemble.r(0) << "v[0] = " << fEnsemble.v(0)
            << "a[0] = " << fEnsemble.a(0) << '\n';
  std::cout << "r[1] = " << fEnsemble.r(1) << "v[1] = " << fEnsemble.v(1)
            << "a[1] = " << fEnsemble.a(1) << '\n';
}

double Sim::GetPotentialEnergy() {
  double result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    for (int j = i + 1; j < fNumAtoms; ++j) {
      result += U_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
    }
  }
  return result;
}

double Sim::GetKineticEnergy() {
  double result{};
  for (int i = 0; i < fNumAtoms; ++i) {
    result += fEnsemble.v(i).sq();
  }
  return result / 2.0;
}

double Sim::GetInternalEnergy() {
  return GetKineticEnergy() + GetPotentialEnergy();
}