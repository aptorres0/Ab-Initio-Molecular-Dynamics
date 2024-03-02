#ifndef _SIM_H
#define _SIM_H

#ifndef _VEC_H
#include "vec.h"
#endif

#ifndef _PARTICLELIST_H
#include "ParticleList.h"
#endif

class Sim {
private:
  double refreshVerletLists(bool, bool);
  double velocityVerletForce();

  double (*pot)(vec), dt, t, r_iter, r_verlet, eKin, ePot, ePotMin,
      histogramLength;
  vec (*f)(vec), simBox, radial;
  int verletSteps, verletUpdate, dim, histogramResolution;
  bool periodic;

public:
  Sim(int, double, vec(vec), double(vec), double, double, int);
  ~Sim() {}

  void initSim(bool, vec, vec(int), vec(int), int, double);
  void velocityVerletStep(bool);
  void updateGraphs();
  void resetGraphs();
  void resetRadialDistribution();
  vec getRadialDistribution();
  int getHistogramResolution();
  void resetDirectionalDistribution();
  mat getDirectionalDistribution();
  double getDt();
  double getT();
  int getDim();
  double getEPotMin();

  ParticleList *particles;
  bool pause;
};

#endif // _SIM_H