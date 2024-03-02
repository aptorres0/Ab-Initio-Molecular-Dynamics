#include "cell.h"
#include "particle.h"
#include <random>

void SimulationCell::Test() {
  _particle[0].position() = Vector3D(1, 2, 3);
  _particle[1].position() = Vector3D(4, 5, 6);
  _particle[0].position().Print();
  _particle[1].position().Print();
}

void SimulationCell::InitRandomState() {
  double l = 1.0;

  std::random_device rd;
  std::mt19937 rng(rd());
  std::normal_distribution<double> v(0.0, 1.0); // to use: d(rng)

  // Initialize positions
  _particle[0].position() = Vector3D(0, 0, 0);
  _particle[1].position() = Vector3D(0, 0.5 * l, 0.5 * l);
  _particle[2].position() = Vector3D(0.5 * l, 0, 0.5 * l);
  _particle[3].position() = Vector3D(0.5 * l, 0.5 * l, 0);

  // Initialize velocities
  _particle[0].velocity() = Vector3D(v(rng), v(rng), v(rng));
  _particle[1].velocity() = Vector3D(v(rng), v(rng), v(rng));
  _particle[2].velocity() = Vector3D(v(rng), v(rng), v(rng));
  _particle[3].velocity() = Vector3D(v(rng), v(rng), v(rng));

  Vector3D p = _particle[0].velocity() + _particle[1].velocity() +
               _particle[2].velocity() + _particle[3].velocity();
  std::cout << "Before Correcting for Momentum: P = ";
  p.Print();
  for (int i = 0; i < 4; i++) {
    std::cout << "\nParticle " << i << ":\n";
    _particle[i].Print();
  }
  _particle[0].velocity() -= 0.25 * p;
  _particle[1].velocity() -= 0.25 * p;
  _particle[2].velocity() -= 0.25 * p;
  _particle[3].velocity() -= 0.25 * p;

  p = _particle[0].velocity() + _particle[1].velocity() +
      _particle[2].velocity() + _particle[3].velocity();
  std::cout << "\nAfter Correcting for Momentum: P = ";
  p.Print();
  for (int i = 0; i < 4; i++) {
    std::cout << "Particle " << i << ":\n";
    _particle[i].Print();
  }
}