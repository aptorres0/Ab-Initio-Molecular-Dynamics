#ifndef _PARTICLELISTENTRY_H
#define _PARTICLELISTENTRY_H

class Particle;
class ParticleList;

class ParticleListEntry {
  friend class ParticleList;

private:
  ParticleList *list;
  Particle *thisParticle;
  ParticleListEntry *next;
  ParticleListEntry *prior;

public:
  ParticleListEntry(Particle *thisParticle_, ParticleList *list_,
                    ParticleListEntry *prior_ = 0, ParticleListEntry *next_ = 0)
      : thisParticle(thisParticle_), list(list_), prior(prior_), next(next_) {}
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

  ParticleListEntry *GetNext() { return next; }
  ParticleListEntry *GetPrior() { return prior; }
  Particle *GetThis() { return thisParticle; }
  bool IsLast() { return !next; }
  bool IsFirst() { return !prior; }
};

#endif // _PARTICLELISTENTRY_H