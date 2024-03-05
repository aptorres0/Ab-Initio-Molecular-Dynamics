#ifndef _UNITCELL_H
#define _UNITCELL_H

#include "Vector3D.h"

typedef struct {
  Vector3D position;
  Vector3D velocity;
  Vector3D Force;
} Atom;

class UnitCell {
  size_t nAtoms = 4; // Number of atoms in cell
  Atom *atom;        // This is an array!

public:
  UnitCell() : atom{new Atom[nAtoms]} {}
  ~UnitCell() { delete[] atom; }

  // Function to set the position vectors of the
  // 4 atoms in the unit cell from an array of double with 12 elements organized
  // by x0,y0,z0,x1,y2,...,y3,z3
  void SetPositions(double p[12]);

  // Function to to set the velocity vectors of
  // the 4 atoms in the unit cell
  void SetVelocities(double v[12]);

  // Function to return position vector for atom `i` by reference
  Vector3D &q(int i);

  // Function to return velocity vector for atom `i` by reference
  Vector3D &p(int i);
};

#endif