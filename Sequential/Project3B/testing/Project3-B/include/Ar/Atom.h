#ifndef _ATOM_H
#define _ATOM_H

#include "../Vector3.h"

namespace Ar {

typedef struct {
  Vector3F position;
  Vector3F velocity;
  Vector3F acceleration;
} Atom;

} // namespace Ar

#endif