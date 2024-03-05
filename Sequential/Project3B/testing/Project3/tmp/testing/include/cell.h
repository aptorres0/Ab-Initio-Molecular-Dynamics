#ifndef _CELL_H
#define _CELL_H

#include "particle.h"

class SimulationCell {
private:
  Particle *_particle; // This is an array!

public:
  SimulationCell() : _particle{new Particle[4]} {}

  // Destructor needed to deallocate memory
  ~SimulationCell() { delete[] _particle; }

  void InitRandomState();
  void Test();
};

#endif