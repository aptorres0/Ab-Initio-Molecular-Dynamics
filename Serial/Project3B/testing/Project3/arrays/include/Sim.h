#ifndef _SIM_H
#define _SIM_H

#include "Ar/Ensemble.h"
#include "Ar/constants.h"
#include "MyRNG.h"
#include "vec.h"
#include <stdarg.h>

class Sim {
private:
  double (*U_ij)(const VectorD &); // Pointer to potential energy function
  // void (*F_ij)(VectorD &, const VectorD &); // Pointer to force function
  VectorD (*F_ij)(const VectorD &); // Pointer to force function
  int fDIM;                         // dimensions of system
  int fNumAtoms;                    // number of atoms
  VectorD fSimBox;                  // lengths of each side of simulation box
  Ar::Ensemble fEnsemble;
  VectorD rij{fDIM};
  VectorD _force{fDIM};
  double h;

public:
  Sim() : fDIM{3}, fNumAtoms{108}, h{0.01} {}

  void InitSim(int _DIM, int _NumAtoms, double (*_U)(const VectorD &),
               VectorD (*_F)(const VectorD &), VectorD _SimBox, double _h) {
    fDIM = _DIM;
    fNumAtoms = _NumAtoms;
    U_ij = _U;
    F_ij = _F;
    fSimBox = _SimBox;
    h = _h;

    // initialize pos, vec, and acc vectors to 0
    for (int i = 0; i < fNumAtoms; ++i) {
      fEnsemble.r(i) = VectorD(3);
      fEnsemble.v(i) = VectorD(3);
      fEnsemble.a(i) = VectorD(3);
    }
  }

  void InitPositions() {
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

  void InitVelocities() {
    /**
     * Steps to initializing the velocities:
     *  1. Set to random values from gaussian centered at 0 with unitless
     * variance of 1
     *  2. Apply correction to set system net momentum to 0
     *  3. Apply correction to set system temperature to 120 K (or 1 in
     * characteristic units)
     */

    // initialize the rng
    // init_random();

    // 1. set velocities from gaussian and compute the resulting momentum
    for (int i = 0; i < fNumAtoms; ++i) {
      fEnsemble.v(i).assign_rnd(fDIM);
    }
    // save that random state
    // save_random();

    VectorD NetMomentum = GetNetMomentum();

    // 2. correct for momentum = 0 by subtracting the mean from each value and
    // compute the resulting temperature
    for (int i = 0; i < fNumAtoms; ++i) {
      fEnsemble.v(i) -= NetMomentum / ((double)fNumAtoms);
    }

    // 3. set temperature of the system by multiplying by a correction factor
    double T_corr = sqrt(1.0 / GetTemperature());
    for (int i = 0; i < fNumAtoms; ++i) {
      fEnsemble.v(i) *= T_corr;
    }
  }

  VectorD GetNetMomentum() const {
    VectorD result = {0.0, 0.0, 0.0};
    for (int i = 0; i < fNumAtoms; ++i) {
      result += fEnsemble.v(i);
    }
    return result;
  }

  double GetTemperature() const {
    double result = 0.0;
    for (int i = 0; i < fNumAtoms; ++i) {
      result += fEnsemble.v(i).sq();
    }
    result /= 3.0 * (fNumAtoms - 1.0);
    return result;
  }

  void ApplyPeriodicBC() {
    for (int i = 0; i < fNumAtoms; ++i) {
      if (IsInSimBox(fEnsemble.r(i)))
        continue;
      fEnsemble.r(i) -= (copysign_v(0.5 * fSimBox, fEnsemble.r(i) - fSimBox) +
                         copysign_v(0.5 * fSimBox, fEnsemble.r(i)));
    }
  }

  VectorD copysign_v(const VectorD &_lhs, const VectorD &_rhs) {
    if (fDIM == 2)
      return VectorD(copysign(_lhs[0], _rhs[0]), copysign(_lhs[1], _rhs[1]));
    return VectorD(copysign(_lhs[0], _rhs[0]), copysign(_lhs[1], _rhs[1]),
                   copysign(_lhs[2], _rhs[2]));
  }

  bool IsInSimBox(const VectorD &check) {
    for (int i = 0; i < fDIM; ++i) {
      if (check[i] > fSimBox[i] || check[i] < 0)
        return false;
    }
    return true;
  }

  bool ImageIsCloser(const VectorD &check) {
    for (int i = 0; i < fDIM; ++i) {
      if (abs(check[i]) > 0.5 * fSimBox[i])
        return true;
    }
    return false;
  }

  const VectorD &MinimumImage(const VectorD &ri, const VectorD &rj) {
    for (int i = 0; i < fDIM; ++i) {
      rij[i] = ri[i] - rj[i] -
               copysign(0.5 * fSimBox[i], ri[i] - rj[i] - 0.5 * fSimBox[i]) -
               copysign(0.5 * fSimBox[i], ri[i] - rj[i] + 0.5 * fSimBox[i]);
    }
    return rij;
  }

  void ComputeForces() {
    // zero them out
    for (int i = 0; i < fNumAtoms; ++i) {
      fEnsemble.a(i).assign(fDIM, 0.0); // avoid a copy
    }
    for (int i = 0; i < fNumAtoms - 1; ++i) {
      for (int j = i + 1; j < fNumAtoms; ++j) {
        _force.assign_to(F_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j))),
                         fDIM); // avoid a copy
        /*F_ij(_force,
             MinimumImage(fEnsemble.r(i), fEnsemble.r(j))); // avoid a copy*/
        fEnsemble.a(i) += _force;
        fEnsemble.a(j) -= _force;
      }
    }
  }

  void VelocityVerletStep() {
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

  void PrintSample() {
    std::cout << "r[0] = " << fEnsemble.r(0) << "v[0] = " << fEnsemble.v(0)
              << "a[0] = " << fEnsemble.a(0) << '\n';
    std::cout << "r[1] = " << fEnsemble.r(1) << "v[1] = " << fEnsemble.v(1)
              << "a[1] = " << fEnsemble.a(1) << '\n';
  }

  double GetPotentialEnergy() {
    double result{};
    // VectorD rij{fDIM};
    for (int i = 0; i < fNumAtoms; ++i) {
      for (int j = i + 1; j < fNumAtoms; ++j) {
        // rij = MinimumImage(fEnsemble.r(i), fEnsemble.r(j));
        result += U_ij(MinimumImage(fEnsemble.r(i), fEnsemble.r(j)));
      }
    }
    return result;
  }

  double GetKineticEnergy() {
    double result{};
    for (int i = 0; i < fNumAtoms; ++i) {
      result += fEnsemble.v(i).sq();
    }
    return result / 2.0;
  }

  double GetInternalEnergy() {
    return GetKineticEnergy() + GetPotentialEnergy();
  }
};

#endif