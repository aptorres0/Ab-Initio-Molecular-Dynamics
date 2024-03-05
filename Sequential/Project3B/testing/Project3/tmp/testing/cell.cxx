#include "cell.h"
#include "MyRNG.h"
#include "Vector3D.h"
#include <cmath>
#include <vector>

void Cell::InitPositions() {
  double l = cbrt(4.0 / 0.8) / 3.0;
  atom[0].position = Vector3D(0.0, 0.0, 0.0);
  atom[1].position = Vector3D(0.0, 0.5 * l, 0.5 * l);
  atom[2].position = Vector3D(0.5 * l, 0.0, 0.5 * l);
  atom[3].position = Vector3D(0.5 * l, 0.5 * l, 0.0);
}

void Cell::InitVelocities() {
  init_random();
  // initialize from gaussian
  atom[0].velocity = Vector3D(rnd(), rnd(), rnd());
  atom[1].velocity = Vector3D(rnd(), rnd(), rnd());
  atom[2].velocity = Vector3D(rnd(), rnd(), rnd());
  atom[3].velocity = Vector3D(rnd(), rnd(), rnd());

  // set net momentum to 0
  Vector3D p0 = Momentum3D();
  atom[0].velocity -= (0.25 * p0);
  atom[1].velocity -= (0.25 * p0);
  atom[2].velocity -= (0.25 * p0);
  atom[3].velocity -= (0.25 * p0);

  // set temperature to 120 K
  double T = Temperature();
  atom[0].velocity *= sqrt(1. / T);
  atom[1].velocity *= sqrt(1. / T);
  atom[2].velocity *= sqrt(1. / T);
  atom[3].velocity *= sqrt(1. / T);

  save_random();
}

void Cell::PeriodicBoundaryCondition() {
  double l = cbrt(4.0 / 0.8) / 3.0;
  double x, y, z, xp, yp, zp;
  for (int i = 0; i < 4; i++) {
    xp = atom[i].position.GetX();
    yp = atom[i].position.GetY();
    zp = atom[i].position.GetZ();
    x = xp;
    y = yp;
    z = zp;

    x = x - copysign(0.5 * l, x - l) - copysign(0.5 * l, x);
    y = y - copysign(0.5 * l, y - l) - copysign(0.5 * l, y);
    z = z - copysign(0.5 * l, z - l) - copysign(0.5 * l, z);

    if ((x > l) || (x < 0) || (y > l) || (y < 0) || (z > l) || (z < 0)) {
      std::cout << "PERIODIC BC IS BROKEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      std::cout << "x', y', z' = " << xp << ", " << yp << ", " << zp << '\n';
      std::cout << "x, y, z = " << x << ", " << y << ", " << z << '\n';
      // exit(0);
    }

    atom[i].position = Vector3D(x, y, z);
  }
}

Vector3D Cell::MIC(Vector3D r12) {
  double l = cbrt(4.0 / 0.8) / 3.0;
  double x = r12.GetX();
  double y = r12.GetY();
  double z = r12.GetZ();

  x = x - copysign(0.5 * l, x - 0.5 * l) - copysign(0.5 * l, x + 0.5 * l);
  y = y - copysign(0.5 * l, y - 0.5 * l) - copysign(0.5 * l, y + 0.5 * l);
  z = z - copysign(0.5 * l, z - 0.5 * l) - copysign(0.5 * l, z + 0.5 * l);
  return Vector3D(x, y, z);
}

void Cell::PrintState() {
  std::cout << "System:\n\tP = ";
  Momentum3D().Print();
  std::cout << "\tT = " << Temperature() << '\n';
  std::cout << "Atom 1:\n\tq = ";
  atom[0].position.Print();
  std::cout << "\tv = ";
  atom[0].velocity.Print();
  std::cout << "Atom 2:\n\tq = ";
  atom[1].position.Print();
  std::cout << "\tv = ";
  atom[1].velocity.Print();
  std::cout << "Atom 3:\n\tq = ";
  atom[2].position.Print();
  std::cout << "\tv = ";
  atom[2].velocity.Print();
  std::cout << "Atom 4:\n\tq = ";
  atom[3].position.Print();
  std::cout << "\tv = ";
  atom[3].velocity.Print();
}
