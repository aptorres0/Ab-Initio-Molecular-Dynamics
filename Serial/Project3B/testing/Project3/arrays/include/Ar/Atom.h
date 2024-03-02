#ifndef _ATOM_H
#define _ATOM_H

#include "vec.h"

namespace Ar {

typedef struct {
  VectorD position;
  VectorD velocity;
  VectorD acceleration;
} Atom;

} // namespace Ar

#endif