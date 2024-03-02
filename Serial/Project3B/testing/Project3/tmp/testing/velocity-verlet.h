#ifndef _VELOCITY_VERLET_H
#define _VELOCITY_VERLET_H

#include "Vector3D.h"
#include "cell.h"
#include <cmath>

class VelocityVerlet {
  Vector3D *F; // this is an array!
  double h = 0.01;
  Cell c;
  double U, T;

public:
  VelocityVerlet() : F{new Vector3D[4]}, U(0.0), T(0.0) {}
  ~VelocityVerlet() { delete[] F; }

  void InitSystem();

  void UpdateSystem();

  void ComputeForce();

  void ComputeU();
  void ComputeT();

  inline double GetE() { return 4 * U + 0.5 * T; }
  inline double GetU() { return 4 * U; }
  inline double GetT() { return 0.5 * T; }

  void PrintState();

  Vector3D GetP(int i) { return c.GetPos(i); }
};

#endif