#ifndef _PARTICLELIST_H
#define _PARTICLELIST_H

#ifndef _PARTICLE_H
#include "Particle.h"
#endif // _PARTICLE_H

#ifndef _PARTICLELISTENTRY_H
#include "ParticleListEntry.h"
#endif // _PARTICLELISTENTRY_H

class ParticleList {
  friend class ParticleListEntry;

private:
  int length;
  ParticleListEntry *last;
  ParticleListEntry *first;

public:
  ParticleList() : first{0}, last{0}, length{0} {}
  ~ParticleList();

  int GetLength() { return length; }
  ParticleListEntry *GetFirst() { return first; }
  ParticleListEntry *GetLast() { return last; }
  void AddParticle(Particle *particle) {
    if (this->last)
      this->last = (this->last->next =
                        new ParticleListEntry(particle, this, this->last));
    else
      this->last =
          (this->first = new ParticleListEntry(particle, this, this->last));
    this->length++;
  }
  void clear() { // ParticleListEntry destructor handles everything
    while (this->first)
      delete this->first;
  }
};

#endif // _PARTICLELIST_H