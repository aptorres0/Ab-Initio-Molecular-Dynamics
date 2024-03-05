#include "velocity-verlet.h"
#include "Vector3D.h"
#include "cell.h"
#include <cmath>

void VelocityVerlet::InitSystem() {
  c.InitPositions();
  c.InitVelocities();
}

void VelocityVerlet::UpdateSystem() {
  // 1. Calculate F(t):
  ComputeForce();
  // 2. Compute r(t+h):
  for (int i = 0; i < 4; ++i) {
    c.p(i) += h * c.GetVel(i) + 0.5 * h * h * F[i];
    c.PeriodicBoundaryCondition();
    c.v(i) += 0.5 * h * F[i];
  }
  // 3. Compute F(t+h):
  ComputeForce();
  // 4. Compute v(t+h):
  for (int i = 0; i < 4; ++i) {
    c.v(i) += 0.5 * h * F[i];
  }
  // 5. Compute the energies
  ComputeU();
  ComputeT();
}

void VelocityVerlet::ComputeForce() {
  for (int i = 0; i < 4; ++i) {
    Vector3D tmp(0.0, 0.0, 0.0);
    for (int j = 0; j < 4; ++j) {
      if (j == i)
        continue;
      Vector3D rij = c.MIC(c.GetPos(i) - c.GetPos(j));
      tmp +=
          rij * (1.0 / pow(rij.sq(), 3.0) - 0.5) * (1.0 / pow(rij.sq(), 4.0));
    }
    F[i] = tmp * 48.0;
  }
}

void VelocityVerlet::ComputeU() {
  U = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < j; ++i) {
      Vector3D rij = c.MIC(c.GetPos(i) - c.GetPos(j));
      U += 1.0 / pow(rij.sq(), 6.0);
      U -= 1.0 / pow(rij.sq(), 3.0);
    }
  }
}
void VelocityVerlet::ComputeT() {
  T = 0;
  for (int i = 0; i < 4; ++i) {
    T += c.GetVel(i).sq();
  }
}

void VelocityVerlet::PrintState() {
  c.PrintState();
  std::cout << "\n\tF1 = ";
  F[0].Print();
  std::cout << "\n\tF2 = ";
  F[1].Print();
  std::cout << "\n\tF3 = ";
  F[2].Print();
  std::cout << "\n\tF4 = ";
  F[3].Print();
}