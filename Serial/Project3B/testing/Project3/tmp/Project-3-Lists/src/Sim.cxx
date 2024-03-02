#include "Sim.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticleListEntry.h"
#include "vec.h"

void Sim::initSim(bool periodic_, vec simBox_, vec (*r0)(int) = 0,
                  vec (*v0)(int) = 0, int histogramResolution_ = 201,
                  double histogramLength_ = 0.5) {
  // restore starting conditions
  // **************************************************
  int i, j, k;
  this->simBox = simBox_;
  MDParticleListEntry *entry;
  if (r0) {
    i = 0;
    for (entry = this->particles->getFirst(); entry; entry = entry->GetNext())
      entry->GetThis()->r = r0(i++);
  } else {
    double n = ceil(pow((double)this->particles->getLength(), 1.0 / this->dim));
    entry = this->particles->getFirst();
    if (this->dim == 2)
      for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
          if (entry) {
            MDParticle *particle = entry->getThis();
            particle->r[1] = (double)i * this->simBox[1] / n;
            if (j & 1)
              particle->r[1] += 1.0 / (2.0 * n);
            particle->r[2] = (double)j * this->simBox[2] / n;
            entry = entry->getNext();
          }
    if (this->dim == 3)
      for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
          for (k = 0; k < n; k++)
            if (entry) {
              MDParticle *particle = entry->getThis();
              particle->r[1] = (double)i * this->simBox[1] / n;
              particle->r[2] = (double)j * this->simBox[2] / n;
              particle->r[3] = (double)k * this->simBox[3] / n;
              entry = entry->getNext();
            }
  }
  if (v0) {
    i = 0;
    for (entry = this->particles->getFirst(); entry; entry = entry->getNext())
      entry->getThis()->v = v0(i++);
  } else {
    for (entry = this->particles->getFirst(); entry; entry = entry->getNext())
      entry->getThis()->v = vec(0, this->dim);
  }
  this->periodic = periodic_;
  this->t = 0;
  this->radial = vec(0, this->histogramResolution);
  this->directional = mat(this->histogramResolution, this->histogramResolution);
  this->ePot = (this->ePotMin = this->refreshVerletLists(true, true));
  // calculate kinetic energy and total momentum
  // vec p_ges(0, this->dim);
  this->eKin = 0;
  for (entry = this->particles->getFirst(); entry; entry = entry->getNext()) {
    // p_ges += entry->getThis()->v;
    this->eKin += entry->getThis()->v * entry->getThis()->v;
  }
  this->eKin *= 0.5;
  this->resetGraphs();
  this->updateGraphs();
  this->pause = true;
}