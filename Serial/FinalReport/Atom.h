#ifndef _ATOM_H
#define _ATOM_H

#include "../Vector3.h"

namespace Ar {

class Atom {
public:
  Vector3D r, v, a;
  int ID;
  Atom *next;

  Atom() : r{}, v{}, a{}, next{nullptr} {}
  Atom(Vector3D _r, Vector3D _v, int _ID, Atom *next = nullptr)
      : r{_r}, v{_v}, a{}, ID{_ID}, next{next} {}
  Atom(Vector3D _r, Vector3D _v, Vector3D _a, int _ID, Atom *next = nullptr)
      : r{_r}, v{_v}, a{_a}, ID{_ID}, next{next} {}
};

} // namespace Ar

#endif