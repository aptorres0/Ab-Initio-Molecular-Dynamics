#include "Particle.h"
#include "ParticleList.h"

class ParticleNode {
  friend class ParticleList;

private:
  Particle *particle;
  ParticleNode *next;
  ParticleNode *prev;
  ParticleList *list;

public:
  ParticleNode() {}
  ParticleNode(Particle *_particle)
      : particle{_particle}, list{NULL}, next{NULL}, prev{NULL} {}
  ParticleNode(Particle *_particle, ParticleList *_list,
               ParticleNode *_prev = NULL, ParticleNode *_next = NULL)
      : particle{_particle}, list{_list}, prev{_prev}, next{_next} {}

  ~ParticleNode();

  ParticleNode *GetNext() { return next; }
  ParticleNode *GetPrev() { return prev; }
  Particle *GetParticle() { return particle; }
  bool IsLast() { return !next; }
  bool IsFirst() { return !prev; }
};

ParticleNode::~ParticleNode() {
  if (this->prev)
    this->prev->next = this->next;
  if (this->next)
    this->next->prev = this->prev;
  if (!this->prior)
    list->first = this->next;
  if (!this->next)
    list->last = this->prev;
  this->list->length--;
}