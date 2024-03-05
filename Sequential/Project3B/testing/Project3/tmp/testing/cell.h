#ifndef _CELL_H
#define _CELL_H

#include "Vector3D.h"

typedef struct {
  Vector3D position;
  Vector3D velocity;
} Atom;

class Cell {
  size_t nAtoms = 4; // Number of atoms in cell
  Atom *atom;        // This is an array!

public:
  Cell() : atom{new Atom[nAtoms]} { /*PrintState();*/
  }
  ~Cell() { delete[] atom; }

  inline Vector3D GetPos(int i) const { return atom[i].position; }
  inline Vector3D GetVel(int i) const { return atom[i].velocity; }
  void SetPos(int i, Vector3D v_to) { atom[i].position = v_to; }
  void SetVel(int i, Vector3D v_to) { atom[i].velocity = v_to; }
  inline Vector3D &p1() { return atom[0].position; }
  inline Vector3D &p2() { return atom[1].position; }
  inline Vector3D &p3() { return atom[2].position; }
  inline Vector3D &p4() { return atom[3].position; }
  inline Vector3D &v1() { return atom[0].velocity; }
  inline Vector3D &v2() { return atom[1].velocity; }
  inline Vector3D &v3() { return atom[2].velocity; }
  inline Vector3D &v4() { return atom[3].velocity; }
  inline Vector3D &p(int i) { return atom[i].position; }
  inline Vector3D &v(int i) { return atom[i].velocity; }

  void PrintState();

  // Compute the unitless momentum vector for the 4 atoms in the finite
  // simulation cell
  inline Vector3D Momentum3D() {
    return atom[0].velocity + atom[1].velocity + atom[2].velocity +
           atom[3].velocity;
  }

  inline double Temperature() {
    return (atom[0].velocity.sq() + atom[1].velocity.sq() +
            atom[2].velocity.sq() + atom[3].velocity.sq()) /
           9.0;
  }

  void InitPositions();

  void InitVelocities();

  void PeriodicBoundaryCondition();

  Vector3D MIC(Vector3D r12);
};

#endif