#ifndef _PARTICLENODE_H
#define _PARTICLENODE_H

#include "Particle.h"
#include "ParticleList.h"

class ParticleNode {
  friend class ParticleList;

private:
  Particle *particle;
  ParticleNode *next;
  ParticleNode *prev;

public:
  ParticleNode() : particle{new Particle()}, next{nullptr}, prev{nullptr} {}
  ParticleNode(Particle *_particle)
      : particle{_particle}, next{NULL}, prev{NULL} {}
  ParticleNode(Particle *_particle, ParticleNode *_prev = NULL,
               ParticleNode *_next = NULL)
      : particle{_particle}, next{_next}, prev{_prev} {}

  ~ParticleNode();

  ParticleNode *GetNext() { return next; }
  ParticleNode *GetPrev() { return prev; }
  Particle *GetParticle() { return particle; }
  Vector3F &r() { return particle->r; }
  Vector3F &v() { return particle->v; }
  Vector3F &a() { return particle->a; }
  const Vector3F &r() const { return particle->r; }
  const Vector3F &v() const { return particle->v; }
  const Vector3F &a() const { return particle->a; }
  bool IsLast() { return !next; }
  bool IsFirst() { return !prev; }
};

ParticleNode::~ParticleNode() {
  if (this->prev)
    this->prev->next = this->next;
  if (this->next)
    this->next->prev = this->prev;
}

#endif