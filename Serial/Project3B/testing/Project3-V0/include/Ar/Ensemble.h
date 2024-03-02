#ifndef _ENSEMBLE_H
#define _ENSEMBLE_H

#include "Atom.h"
#include "vec.h"

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
  VectorD &r(int i) { return fAtom[i].position; }
  const VectorD &r(int i) const { return fAtom[i].position; } // const version

  // get address of velocity vector for atom `i` by reference for editing
  VectorD &v(int i) { return fAtom[i].velocity; }
  const VectorD &v(int i) const { return fAtom[i].velocity; } // const version

  // get address of acceleration vector for atom `i` by reference for editing
  VectorD &a(int i) { return fAtom[i].acceleration; }
  const VectorD &a(int i) const {
    return fAtom[i].acceleration;
  } // const version
};

} // namespace Ar

#endif