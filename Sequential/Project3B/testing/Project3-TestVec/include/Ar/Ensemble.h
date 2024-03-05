#ifndef _ENSEMBLE_H
#define _ENSEMBLE_H

#include "../vec.h"
#include "Atom.h"

namespace Ar {
class Ensemble {
private:
  const int fNumAtoms; // number of atoms
  Atom *fAtom;         // This is an array of atoms

public:
  Ensemble(int _NumAtoms = 108)
      : fNumAtoms{_NumAtoms}, fAtom{new Atom[fNumAtoms]} {}
  ~Ensemble() { delete[] fAtom; }

  // Get address of position vector for atom `i` by reference for editing
  Vector3D &r(int i) { return fAtom[i].position; }
  const Vector3D &r(int i) const { return fAtom[i].position; } // const version

  // get address of velocity vector for atom `i` by reference for editing
  Vector3D &v(int i) { return fAtom[i].velocity; }
  const Vector3D &v(int i) const { return fAtom[i].velocity; } // const version

  // get address of acceleration vector for atom `i` by reference for editing
  Vector3D &a(int i) { return fAtom[i].acceleration; }
  const Vector3D &a(int i) const {
    return fAtom[i].acceleration;
  } // const version
};

} // namespace Ar

#endif