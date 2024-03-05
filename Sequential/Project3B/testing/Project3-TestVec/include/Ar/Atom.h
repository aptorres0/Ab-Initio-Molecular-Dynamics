#ifndef _ATOM_H
#define _ATOM_H

#include "../vec.h"

namespace Ar {

typedef struct {
  Vector3D position;
  Vector3D velocity;
  Vector3D acceleration;
} Atom;

} // namespace Ar

#endif