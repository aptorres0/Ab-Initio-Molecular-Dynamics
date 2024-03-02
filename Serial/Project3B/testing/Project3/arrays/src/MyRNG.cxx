#include "MyRNG.h"

void init_random() {
  // This function tries to load the previous saved state of the the random
  // number generator. Otherwise it seeds the generator from a non-deterministic
  // source. Call this function at the start of the program.
  std::ifstream frand("mtrand_internal_state.dat");
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
  std::ofstream frand("mtrand_internal_state.dat");
  if (frand.is_open())
    frand << rng;
  frand.close();
}