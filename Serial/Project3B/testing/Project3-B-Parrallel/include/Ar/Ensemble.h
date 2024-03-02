#ifndef _ENSEMBLE_H
#define _ENSEMBLE_H

#include "../Vector3.h"
#include "Atom.h"

namespace Ar {
class Ensemble {
private:
  // const int fNumAtoms = 108; // number of atoms
  Atom fAtom[108]; // This is an array of atoms

public:
  Ensemble() {}

  // Get address of position vector for atom `i` by reference for editing
  Vector3F &r(int i) { return fAtom[i].position; }
  const Vector3F &r(int i) const { return fAtom[i].position; } // const version

  // get address of velocity vector for atom `i` by reference for editing
  Vector3F &v(int i) { return fAtom[i].velocity; }
  const Vector3F &v(int i) const { return fAtom[i].velocity; } // const version

  // get address of acceleration vector for atom `i` by reference for editing
  Vector3F &a(int i) { return fAtom[i].acceleration; }
  const Vector3F &a(int i) const {
    return fAtom[i].acceleration;
  } // const version
};

} // namespace Ar

#endif