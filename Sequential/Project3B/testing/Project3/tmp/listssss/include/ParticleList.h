#include "ParticleNode.h"

class ParticleList {
private:
  int length;
  ParticleNode *first;
  ParticleNode *last;

public:
  ParticleList() : first{NULL}, last{NULL}, length{0} {}
  ~ParticleList() { this->clear(); }

  int GetLength() { return length; }
  ParticleNode *GetFirst() { return first; }
  ParticleNode *GetLast() { return last; }
  void AddParticle(Particle *in);
  void clear();
};
/**
 * Add particle node to the doubly linked list
 */
void ParticleList::AddParticle(Particle *in) {
  ParticleNode *NewNode = new ParticleNode(in); // create new node

  /*// check if list is empty:
  if (this->last) {
    // adjust the links:
    // 0. store the address of the last node to `prev` of the new node
    NewNode->prev = this->last;
    // 1. store the address of the head node to `next` of the new node
    NewNode->next = this->last->next;
    // 2. point the current last node to the new node
    this->last->next = NewNode;
    // 3. make the new node the last node
    this->last = NewNode;
  } else {
    // if the list is empty, make this node first and last and store the address
    // of the head node to `next` of the last node
    // 0. store the address of the last node to `prev` of the new node
    NewNode->prev = this->last;
    // 1. assign new node to head of the list
    this->first = NewNode;
    // 2. assign new node to tail of the list
    this->last = NewNode;
  }*/

  if (this->last) // if `last` is not NULL, append to list
    this->last = (this->last->next = new ParticleNode(in, this, this->last));
  else // if `tail` is NULL
    this->last = (this->first = new ParticleNode(in, this, this->last));
  this->length++;
}

void ParticleList::clear() {
  while (this->first)
    delete this->first;
}