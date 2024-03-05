#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "Vector3D.h"

typedef struct {
  Vector3D position;
  Vector3D velocity;

  void Print() {
    std::cout << "q = (" << position.GetX() << ", " << position.GetY() << ", "
              << position.GetZ() << ")\n";
    std::cout << "v = (" << velocity.GetX() << ", " << velocity.GetY() << ", "
              << velocity.GetZ() << ")\n";
  }
} Atom;

class Particle {
private:
  Vector3D _position;
  Vector3D _velocity;
  Vector3D _force;

public:
  Particle() : _position{}, _velocity{}, _force{} {}

  inline Vector3D position() const { return _position; }
  inline Vector3D &position() { return _position; }
  inline Vector3D velocity() const { return _velocity; }
  inline Vector3D &velocity() { return _velocity; }
  inline Vector3D force() const { return _force; }
  inline Vector3D &force() { return _force; }

  void SetPosition(Vector3D r) { _position = r; }
  void SetVelocity(Vector3D v) { _velocity = v; }
  void SetForce(Vector3D f) { _force = f; }

  void Print() {
    std::cout << "\tP = ";
    _position.Print();
    std::cout << "\tV = ";
    _velocity.Print();
  }
};

#endif