#include "Ar/system.h"
#include "MyRNG.h"
#include "unit_cell.h"
#include <cmath>

namespace Ar {

Vector3D &System::q(int i) { return atom[i].position; }
Vector3D &System::p(int i) { return atom[i].velocity; }
Vector3D System::q(int i) const { return atom[i].position; }
Vector3D System::p(int i) const { return atom[i].velocity; }
Vector3D &System::F(int i) { return atom[i].Force; }
Vector3D System::F(int i) const { return atom[i].Force; }

void System::InitPositions() {
  // Generate position vectors for the 3D Volume of 9 unit cells
  for (int i = 0; i < 3; i++) { // loop over cells in x
    double xmin = i * l;
    for (int j = 0; j < 3; j++) { // loop over cells in y
      double ymin = j * l;
      for (int k = 0; k < 3; k++) { // loop over cells in z
        double zmin = k * l;
        // NOTE: this should be the same as q(i+j+k) = ...
        /*atom[i * 36 + j * 12 + k * 4].position = Vector3D(xmin, ymin, zmin);
        atom[i * 36 + j * 12 + k * 4 + 1].position =
            Vector3D(xmin, 0.5 * l + ymin, 0.5 * l + zmin);
        atom[i * 36 + j * 12 + k * 4 + 2].position =
            Vector3D(0.5 * l + xmin, ymin, 0.5 * l + zmin);
        atom[i * 36 + j * 12 + k * 4 + 3].position =
            Vector3D(0.5 * l + xmin, 0.5 * l + ymin, zmin);*/
        q(i * 36 + j * 12 + k * 4) = Vector3D(xmin, ymin, zmin);
        q(i * 36 + j * 12 + k * 4 + 1) =
            Vector3D(xmin, 0.5 * l + ymin, 0.5 * l + zmin);
        q(i * 36 + j * 12 + k * 4 + 2) =
            Vector3D(0.5 * l + xmin, ymin, 0.5 * l + zmin);
        q(i * 36 + j * 12 + k * 4 + 3) =
            Vector3D(0.5 * l + xmin, 0.5 * l + ymin, zmin);
      }
    }
  }
}

Vector3D System::GetNetMomentum() const {
  Vector3D result(0, 0, 0);
  for (int i = 0; i < nAtoms; ++i) {
    // result += atom[i].velocity;
    result += p(i);
  }

  return result;
}

double System::GetTemperature() const {
  double result = 0.0;
  for (int i = 0; i < nAtoms; ++i) {
    // result += atom[i].velocity.sq();
    result += p(i).sq();
  }
  result /= 3.0 * (nAtoms - 1.0);
  return result;
}

void System::InitVelocities() {
  /**
   * Steps to initializing the velocities:
   *  1. Set to random values from gaussian centered at 0 with unitless variance
   *     of 1
   *  2. Apply correction to set system net momentum to 0
   *  3. Apply correction to set system temperature to 120 K (or 1 in
   * characteristic units)
   */

  // initialize the rng
  init_random();

  // 1. set velocities from gaussian and compute the resulting momentum
  Vector3D lNetMomentum(0.0, 0.0, 0.0);
  for (int i = 0; i < nAtoms; ++i) {
    /*atom[i].velocity = Vector3D(rnd(), rnd(), rnd());
    lNetMomentum += atom[i].velocity;*/
    p(i) = Vector3D(rnd(), rnd(), rnd());
    lNetMomentum += p(i);
  }
  // save that random state
  save_random();

  // 2. correct for momentum = 0 by subtracting the mean from each value and
  // compute the resulting temperature
  double Temp = 0;
  for (int i = 0; i < nAtoms; ++i) {
    /*atom[i].velocity -= lNetMomentum / ((double)nAtoms);
    Temp += atom[i].velocity.sq();*/
    p(i) -= lNetMomentum / ((double)nAtoms);
    Temp += p(i).sq();
  }
  Temp /= 3.0 * (nAtoms - 1.0);

  // 3. set temperature of the system by multiplying by a correction factor
  double T_corr = sqrt(1.0 / Temp);
  for (int i = 0; i < nAtoms; ++i) {
    // atom[i].velocity *= T_corr;
    p(i) *= T_corr;
  }
}

void System::ApplyPeriodicBoundaryCondition() {
  for (int i = 0; i < nAtoms; ++i) {
    double xyz[3] = {q(i).GetX(), q(i).GetY(), q(i).GetZ()};
    xyz[0] -= (copysign(0.5 * L, xyz[0] - L) + copysign(0.5 * L, xyz[0]));
    xyz[1] -= (copysign(0.5 * L, xyz[1] - L) + copysign(0.5 * L, xyz[1]));
    xyz[2] -= (copysign(0.5 * L, xyz[2] - L) + copysign(0.5 * L, xyz[2]));
    q(i) = Vector3D(xyz[0], xyz[1], xyz[2]);
  }
  // auto BC = [L](double &x) {
  //   x -= copysign(0.5 * L, x - L) + copysign(0.5 * L, x);
  // };
}

Vector3D System::MIC(Vector3D r12) {
  double xyz[3] = {r12.GetX(), r12.GetY(), r12.GetZ()};
  // std::cout << xyz[0] << "," << xyz[1] << "," << xyz[2] << '\n';
  xyz[0] -= (copysign(0.5 * L, xyz[0] - 0.5 * L) +
             copysign(0.5 * L, xyz[0] + 0.5 * L));
  xyz[1] -= (copysign(0.5 * L, xyz[1] - 0.5 * L) +
             copysign(0.5 * L, xyz[1] + 0.5 * L));
  xyz[2] -= (copysign(0.5 * L, xyz[2] - 0.5 * L) +
             copysign(0.5 * L, xyz[2] + 0.5 * L));
  // std::cout << xyz[0] << "," << xyz[1] << "," << xyz[2] << '\n';
  return Vector3D(xyz[0], xyz[1], xyz[2]);
}

void System::ComputeForces() {
  for (int i = 0; i < nAtoms; ++i) {
    F(i) = Vector3D(0, 0, 0);
    for (int j = 0; j < nAtoms; ++j) {
      if (j == i)
        continue;
      F(i) += (1.0 / MIC(q(i) - q(j)).pow7() -
               0.5 * 1.0 / MIC(q(i) - q(j)).pow4()) *
              MIC(q(i) - q(j));
    }
    F(i) *= 48.0;
  }
}

void System::Step() {
  double h = 0.01;

  ComputeForces();

  for (int i = 0; i < nAtoms; ++i) {
    q(i) += (h * p(i) + 0.5 * h * h * F(i));
    p(i) += 0.5 * h * F(i);
  }

  ApplyPeriodicBoundaryCondition();

  ComputeForces();

  for (int i = 0; i < nAtoms; ++i)
    p(i) += 0.5 * h * F(i);
}

double System::GetU() {
  double res = 0;
  for (int i = 0; i < nAtoms; ++i) {
    for (int j = i + 1; j < nAtoms; ++j) {
      res += 1.0 / MIC(q(i) - q(j)).pow6() - 1.0 / MIC(q(i) - q(j)).pow3();
    }
  }
  res *= 4;
  return res;
}

double System::GetT() {
  double res = 0;
  for (int i = 0; i < nAtoms; ++i) {
    res += p(i).sq();
  }
  res *= 0.5;
  return res;
}

double System::GetE() { return (GetT() + GetU()); }

} // namespace Ar