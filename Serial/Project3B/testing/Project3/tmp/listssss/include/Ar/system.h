#ifndef _SYSTEM_H
#define _SYSTEM_H

#include "Ar/constants.h"
#include "unit_cell.h"
#include <cmath> // cbrt

namespace Ar {

// TODO: Check profiler to see if array of atoms is faster than std::vector
class System {
private:
  int nAtoms = 108; // Number of atoms in cell
  double L;         // Dimensionless length of the Volume
  double l;         // Dimensionless length of the unit cell

public:
  Atom *atom; // This is an array!
  System() : atom{new Atom[nAtoms]} {
    // dimensionless number density of argon atoms
    double n = Ar::gNumberDensity_r;
    L = cbrt(nAtoms / n);
    l = L / 3.0;
  }
  ~System() { delete[] atom; }

  // Function to return position vector for atom `i` by reference
  Vector3D &q(int i);

  // Function to return velocity vector for atom `i` by reference
  Vector3D &p(int i);

  Vector3D &F(int i);
  Vector3D F(int i) const;

  // Return position vector for atom `i` by value
  Vector3D q(int i) const;
  // Return momentum vector for atom `i` by value
  Vector3D p(int i) const;

  Vector3D GetPosition(int i) const { return atom[i].position; }

  void InitPositions();

  void InitVelocities();

  Vector3D GetNetMomentum() const;

  double GetTemperature() const;

  void ApplyPeriodicBoundaryCondition();

  void ComputeForces();

  Vector3D MIC(Vector3D r12);

  void VelocityVerletIntegration();

  void Step();

  double GetE();

  double GetU();

  double GetT();
};

} // namespace Ar

#endif