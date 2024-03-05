#include "unit_cell.h"

void UnitCell::SetPositions(double p[12]) {
  for (int i = 0; i < 4; ++i) {
    atom[i].position = Vector3D(p[i * 3], p[i * 3 + 1], p[i * 3 + 2]);
  }
}

void UnitCell::SetVelocities(double v[12]) {
  for (int i = 0; i < 4; ++i) {
    atom[i].velocity = Vector3D(v[i * 3], v[i * 3 + 1], v[i * 3 + 2]);
  }
}

Vector3D &UnitCell::q(int i) { return atom[i].position; }
Vector3D &UnitCell::p(int i) { return atom[i].velocity; }