#ifndef _MYRNG_H
#define _MYRNG_H

#include <algorithm> // std::generate_n
#include <array>
#include <fstream>
#include <functional> // std::bind
#include <random>

using std::mt19937;
using std::normal_distribution;
using std::uniform_real_distribution;

static mt19937 rng; // mersenne twistor random number generator
static auto rnd = std::bind(
    normal_distribution<double>(0.0, 1.0),
    std::ref(
        rng)); // random number drawn from gaussian with mean = 0 and variance 1

static auto rndu = std::bind(
    uniform_real_distribution<double>(0.0, 1.0),
    std::ref(
        rng)); // random number drawn from uniform distribution between 0 and 1

void init_random() {
  // This function tries to load the previous saved state of the the random
  // number generator. Otherwise it seeds the generator from a non-deterministic
  // source. Call this function at the start of the program.
  std::ifstream frand("/Users/ajpt/.local/state/rng_internal_state.txt");
  if (frand.is_open())
    frand >> rng;
  else {
    std::array<int, 624> seed_data;
    std::random_device r;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq);
  }
  frand.close();
}

void save_random() {
  // This function saves internal state of the random number generator.
  // Call it just before the program terminates.
  std::ofstream frand("/Users/ajpt/.local/state/rng_internal_state.txt");
  if (frand.is_open())
    frand << rng;
  frand.close();
}

#endif