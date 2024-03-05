#ifndef _FCC_H
#define _FCC_H

#include "particle.h"

class Cell {
  size_t _nAtoms;
  Atom *_atom; // This is an array!!

public:
  Cell(size_t N = 4) : _nAtoms{N}, _atom{new Atom[N]} { InitRandomState(); }
}

// Need a destructor to free memory when object
// goes out of scope
~Cell() {
  delete[] _atom;
}

void InitRandomState();

#endif