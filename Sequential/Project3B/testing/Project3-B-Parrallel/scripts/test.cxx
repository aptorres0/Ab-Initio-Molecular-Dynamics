#include "../include/vec.h"

Vector3D _tmpf{3};

const Vector3D &testforce(const Vector3D &rij) {
  std::cout << "Computing force\n";
  _tmpf = (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
  std::cout << "\tReturning result\n";
  return _tmpf;
}

Vector3D testforce1(const Vector3D &rij) {
  std::cout << "Computing force\n";
  _tmpf = rij * (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5);
  std::cout << "\tReturning result\n";
  return _tmpf;
}

Vector3D testforce2(const Vector3D &rij) {
  std::cout << "\tReturning result\n";
  return rij * (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5);
}

Vector3D testforce3(const Vector3D &rij) {
  std::cout << "\tReturning result\n";
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}

int main() {
  /*Vector3D a = {1.0, 2.0, 3.0};
  std::cout << "a = " << a << "\n";
  Vector3D b = a;
  std::cout << "b = " << b << "\n";*/

  std::cout << "Declaring b\n";
  Vector3D b = {1.0, 2.0, 3.0};
  std::cout << "\n\n";

  Vector3D f0{3};
  Vector3D f1{3};
  Vector3D f2{3};

  std::cout << "f0 = testforce(b)\n";
  f0 = testforce(b);
  std::cout << "\n\n";

  std::cout << "f1 = testforce1(b)\n";
  f1 = testforce1(b);
  std::cout << "\n\n";

  std::cout << "f2 = testforce2(b)\n";
  f2 = testforce2(b);
  std::cout << "\n\n";

  std::cout << "f3 = testforce3(b)\n";
  f2 = testforce3(b);
  std::cout << "\n\n";

  /*std::cout << "Declaring _force and tmp\n";
  Vector3D _force{3};
  Vector3D tmp{3};
  std::cout << "\n\n";

  std::cout << "Calling: _force = testforce(b)\n";
  _force = testforce(b);
  std::cout << "\n\n";

  std::cout << "Calling tmp.assing_to(_force,3)\n";
  tmp.assign_to(_force, 3);
  std::cout << "\n\n";

  std::cout << "Calling tmp.assing_to(testforce(b),3)\n";
  tmp.assign_to(testforce(b), 3);
  std::cout << "\n\n";

  std::cout << "Calling Vector3D c = testforce(b)\n";
  Vector3D c = testforce(b);*/

  /*std::cout << "Declaring arr\n";
  Vector3D *arr = new Vector3D[3];
  std::cout << "Calling testforce\n";
  arr[0] = testforce(b);
  arr[1] = Vector3D(3);
  std::cout << "Using assign_to\n";
  arr[1].assign_to(_force, 3);*/

  return 0;
}