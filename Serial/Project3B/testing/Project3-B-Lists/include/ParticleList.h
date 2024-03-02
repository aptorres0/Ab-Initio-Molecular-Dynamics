#ifndef _PARTICLELIST_H
#define _PARTICLELIST_H

#include "Particle.h"
#include "ParticleNode.h"

class ParticleList {
  friend class ParticleNode;

private:
  ParticleNode *first;
  ParticleNode *last;
  int length = 0;

public:
  ParticleList() : first{nullptr}, last{nullptr}, length{0} {}
  ~ParticleList() { this->clear(); }

  int GetLength() { return length; }
  ParticleNode *GetFirst() { return first; }
  ParticleNode *GetLast() { return last; }
  void AddParticle(Particle *in);
  void clear();
};
/**
 * Add particle node to the linked list
 */
void ParticleList::AddParticle(Particle *in) {
  /*std::cout << "Inside ParticleList::AddParticle\n";
  std::cout << "this->length = " << this->GetLength() << std::endl;
  if (GetLast() == nullptr)
    std::cout << "\t this->last = nullptr" << std::endl;
  else
    std::cout << "\t this->last != nullptr" << std::endl;*/
  if (this->last) { // if `last` is not NULL, append to list
    this->last = (this->last->next = new ParticleNode(in, this->last));
    /*std::cout << "added particle to list with r = " << in->r << "v = " <<
       in->v
              << "a = " << in->a << std::endl;*/
  } else { // if `tail` is NULL
    this->last = (this->first = new ParticleNode(in, this->last));
    /*std::cout << "-----> Successfully added particle to list with r = " <<
       in->r
              << "v = " << in->v << "a = " << in->a << std::endl;*/
  }
  this->length++;
}

void ParticleList::clear() {
  while (this->first)
    delete this->first;
}

#endif