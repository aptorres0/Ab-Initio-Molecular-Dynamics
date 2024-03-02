#ifndef _PARTICLENODE_H
#define _PARTICLENODE_H

#include "Particle.h"
#include "ParticleList.h"

/**
 * Doubly Linked list node class
 */
class ParticleNode {
  friend class ParticleList;

private:
  // node data element
  Particle *particle;
  // Doubly linked nodes
  ParticleNode *next;
  ParticleNode *prev;

public:
  ParticleListEntry(Particle *_thisParticle, ParticleListEntry *_prior = 0,
                    ParticleListEntry *_next = 0)
      : particle{_thisParticle}, prior{_prior}, next{_next} {}
  ~ParticleListEntry() {
    if (this->prior)
      this->prior->next = this->next;
    if (this->next)
      this->next->prior = this->prior;
    if (!this->prior)
      list->first = this->next;
    if (!this->next)
      list->last = this->prior;
    this->list->length--;
  }

  ParticleNode *GetNext() { return next; }
  ParticleNode *GetPrev() { return prev; }
  Particle *GetParticle() { return particle; }
  bool IsLast() { return !next; }
  bool IsFirst() { return !prior; }
};

#endif