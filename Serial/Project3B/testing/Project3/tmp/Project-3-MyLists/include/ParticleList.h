#ifndef _PARTICLELIST_H
#define _PARTICLELIST_H

#include "Particle.h"
#include "ParticleNode.h"

/**
 * Doubly Linked List class
 */
class ParticleList {
  friend class ParticleNode;

private:
  int length;
  ParticleNode *last;
  ParticleNode *first;

public:
  ParticleList() : first{0}, last{0}, length{0} {}
  ~ParticleList() { this->clear(); }

  int GetLength() { return length; }
  ParticleNode *GetFirst() { return first; }
  ParticleNode *GetLast() { return last; }
  void AddParticle(Particle *particle) {
    if (this->last)
      this->last =
          (this->last->next = new ParticleNode(particle, this, this->last));
    else
      this->last = (this->first = new ParticleNode(particle, this, this->last));
    this->length++;
  }
  void clear() { // ParticleListEntry destructor handles everything
    while (this->first)
      delete this->first;
  }
};

#endif // _PARTICLELIST_H