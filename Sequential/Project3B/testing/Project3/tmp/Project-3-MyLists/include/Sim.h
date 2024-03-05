#ifndef _SIM_H
#define _SIM_H

#include "ParticleList.h"
#include "Vector3D.h"
#include "vec.h"

class Sim {
public:
  Sim(double _dt, Vector3D _SimBox, Vector3D (*Force_)(Vector3D),
      double (*Potential_)(Vector3D));
  ~Sim() {}

  void InitPositions();
  void InitVelocities();
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
  double getEPotMin();

  ParticleList *particles;
  bool pause;

private:
  double refreshVerletLists(bool, bool);
  double velocityVerletForce();

  double (*pot)(vec), dt, t, r_iter, r_verlet, eKin, ePot, ePotMin,
      histogramLength;
  vec (*f)(vec), simBox, radial;
  int verletSteps, verletUpdate, dim, histogramResolution;
};

double Sim::refreshVerletLists(bool calc, bool countRadial) {
  int i;
  double r_abs, pot_ = 0;
  Vector3D r_, Force;
  ParticleListEntry *entry;
  if (calc) // Set acceleration for all particles to zero
    for (entry = this->particles->getFirst(); entry; entry = entry->GetNext())
      entry->GetThis()->a *= 0;
  for (entry = this->particles->getFirst(); entry; entry = entry->GetNext()) {
    entry->GetThis()->verletList->clear();
    ParticleListEntry *verlet;
  }
}

#endif // _SIM_H