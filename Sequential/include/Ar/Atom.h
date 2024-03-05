#ifndef _ATOM_H
#define _ATOM_H
/**
 * Atom.h
 *
 * Author: Alexander Paul Torres
 * Date: 14 APR 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Header file with declaration of the Atom class used to represent the
 * velocity, position, and acceleration vectors of an atom in the simulation.
 * Intended to be used as a linked list.
 */

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