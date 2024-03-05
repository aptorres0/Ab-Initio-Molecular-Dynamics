#include "FCC.h"
#include "particle.h"
#include <random>

void Cell::InitRandomState() {
  std::random_device rd;
  std::mt19937 rng(rd());
  std::normal_distribution<double> v(0.0, 1.0); // to use: d(rng)
}