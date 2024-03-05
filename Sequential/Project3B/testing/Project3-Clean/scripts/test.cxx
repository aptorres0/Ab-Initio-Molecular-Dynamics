#include "../include/vec.h"

VectorD _tmpf{3};

const VectorD &testforce(const VectorD &rij) {
  std::cout << "Computing force\n";
  _tmpf = (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
  std::cout << "\tReturning result\n";
  return _tmpf;
}

VectorD testforce1(const VectorD &rij) {
  std::cout << "Computing force\n";
  _tmpf = rij * (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5);
  std::cout << "\tReturning result\n";
  return _tmpf;
}

VectorD testforce2(const VectorD &rij) {
  std::cout << "\tReturning result\n";
  return rij * (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5);
}

VectorD testforce3(const VectorD &rij) {
  std::cout << "\tReturning result\n";
  return (48.0 / pow8(rij)) * (1.0 / pow6(rij) - 0.5) * rij;
}

int main() {
  /*VectorD a = {1.0, 2.0, 3.0};
  std::cout << "a = " << a << "\n";
  VectorD b = a;
  std::cout << "b = " << b << "\n";*/

  std::cout << "Declaring b\n";
  VectorD b = {1.0, 2.0, 3.0};
  std::cout << "\n\n";

  VectorD f0{3};
  VectorD f1{3};
  VectorD f2{3};

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
  VectorD _force{3};
  VectorD tmp{3};
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

  std::cout << "Calling VectorD c = testforce(b)\n";
  VectorD c = testforce(b);*/

  /*std::cout << "Declaring arr\n";
  VectorD *arr = new VectorD[3];
  std::cout << "Calling testforce\n";
  arr[0] = testforce(b);
  arr[1] = VectorD(3);
  std::cout << "Using assign_to\n";
  arr[1].assign_to(_force, 3);*/

  return 0;
}