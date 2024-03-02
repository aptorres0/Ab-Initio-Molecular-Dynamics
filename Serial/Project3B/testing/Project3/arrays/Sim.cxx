#include "Sim.h"

void Sim::InitSim(int _DIM, int _NumAtoms, double (*_U)(const VectorD &),
                  const VectorD (*_F)(const VectorD &), VectorD _SimBox,
                  double _h) {
  fDIM = _DIM;
  fNumAtoms = _NumAtoms;
  U_ij = _U;
  F_ij = _F;
  fSimBox = _SimBox;
  h = _h;

  // initialize pos, vec, and acc vectors to 0
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.r(i) = VectorD(fDIM);
    fEnsemble.v(i) = VectorD(fDIM);
    fEnsemble.a(i) = VectorD(fDIM);
  }
  /*for (int j = 0; j < 5779; ++j)
    Forces[j] = VectorD(3);*/
}

void Sim::InitPositions() {
  double l = fSimBox[0] / 3.0;
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * fSimBox[0] / 3.0;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * fSimBox[1] / 3.0;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * fSimBox[2] / 3.0;
        // Define one unit cell at a time
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

  // 1. set velocities from gaussian and compute the resulting momentum
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i).assign_rnd(fDIM);
  }

  // 2. correct for momentum = 0 by subtracting the mean from each atom's
  // value
  VectorD NetMomentum = GetNetMomentum();
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) -= NetMomentum / ((double)fNumAtoms);
  }

  // 3. set temperature of the system by multiplying by a correction factor
  // double T_corr = 1.0 / sqrt(GetTemperature());
  double T_corr = sqrt(1.0 / GetTemperature());
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.v(i) *= T_corr;
  }
}

VectorD Sim::GetNetMomentum() const {
  VectorD result = {0.0, 0.0, 0.0};
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

const VectorD &Sim::copysign_v(const VectorD &_lhs, const VectorD &_rhs) {
  for (int i = 0; i < fDIM; ++i) {
    copysign_vect[i] = copysign(_lhs[i], _rhs[i]);
  }
  return copysign_vect;
}

bool Sim::IsInSimBox(const VectorD &check) {
  for (int i = 0; i < fDIM; ++i) {
    if (check[i] > fSimBox[i] || check[i] < 0)
      return false;
  }
  return true;
}

bool Sim::ImageIsCloser(const VectorD &check) {
  for (int i = 0; i < fDIM; ++i) {
    if (abs(check[i]) > 0.5 * fSimBox[i])
      return true;
  }
  return false;
}

const VectorD &Sim::MinimumImage(const VectorD &ri, const VectorD &rj) {
  for (int i = 0; i < fDIM; ++i) {
    rij[i] = ri[i] - rj[i] -
             copysign(0.5 * fSimBox[i], ri[i] - rj[i] - 0.5 * fSimBox[i]) -
             copysign(0.5 * fSimBox[i], ri[i] - rj[i] + 0.5 * fSimBox[i]);
  }
  return rij;
}

void Sim::ComputeForces() {
  // int cnt = 0;
  for (int i = 0; i < fNumAtoms - 1; ++i) {
    fEnsemble.a(i).assign(fDIM, 0.0); // avoid a copy
    for (int j = i + 1; j < fNumAtoms; ++j) {
      if (i == j)
        continue;
      fEnsemble.a(i) += F_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
      // cnt++;
    }
  }
  // std::cout << cnt << '\n';
}

void Sim::NewComputeForces() {
  // zero them out
  for (int i = 0; i < fNumAtoms; ++i) {
    fEnsemble.a(i).assign(fDIM, 0.0); // avoid a copy
  }
  for (int i = 0; i < fNumAtoms - 1; ++i) {
    for (int j = i + 1; j < fNumAtoms; ++j) {
      _Force = F_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
      fEnsemble.a(i) += _Force;
      fEnsemble.a(j) -= _Force;
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
  for (int i = 0; i < fNumAtoms - 1; ++i) {
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