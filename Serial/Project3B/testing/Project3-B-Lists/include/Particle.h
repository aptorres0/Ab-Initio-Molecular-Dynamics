#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "Vector3.h"

class Particle {
public:
  Vector3F r, v, a;
  Particle *next; // pointer to the next particle in the list

  Particle() : r{}, v{}, a{}, next{nullptr} {}
  Particle(Vector3F _r, Vector3F _v) : r{_r}, v{_v}, a{}, next{nullptr} {}
};

// see The Cherno "Writing an iterator in c++" video

class ParticleIterator {
public:
  using ValueType = Particle;
  using PointerType = Particle *;
  using ReferenceType = Particle &;

public:
  ParticleIterator(PointerType ptr) : m_Ptr{ptr} {}

  ParticleIterator &operator++() {
    m_Ptr++;
    return *this;
  }

  // copy
  ParticleIterator operator++(int) {
    ParticleIterator iterator = *this;
    ++(*this);
    return iterator;
  }

  ParticleIterator &operator--() {
    m_Ptr--;
    return *this;
  }

  ParticleIterator operator--(int) {
    ParticleIterator iterator = *this;
    --(*this);
    return iterator;
  }

  ReferenceType operator[](int index) { return *(m_Ptr + index); }

  PointerType operator->() { return m_Ptr; }

  ReferenceType operator*() { return *m_Ptr; }

  bool operator==(const ParticleIterator &other) const {
    return m_Ptr == other.m_Ptr;
  }

  bool operator!=(const ParticleIterator &other) const {
    return !(*this == other);
  }

private:
  PointerType m_Ptr;
};

class ParticleListt {
public:
  using Iterator = ParticleIterator;

public:
  ParticleListt() : head{nullptr}, fNumParticles{0} {}
  void clear() {
    while (head != nullptr) {
      Particle *current = head;
      head = head->next;
      delete current;
    }
  }
  ~ParticleListt() { clear(); }

  void AddParticle(Particle *p) {
    if (head == nullptr) {
      head = p;
    } else {
      Particle *current = head;
      while (current->next != nullptr) {
        current = current->next;
      }
      current->next = p;
    }
    fNumParticles++;
  }
  void RemoveParticle(Particle *p) {
    if (head == nullptr) {
      return;
    } else if (head == p) {
      head = head->next;
    } else {
      Particle *current = head;
      while (current->next != nullptr) {
        if (current->next == p) {
          current->next = current->next->next;
          break;
        }
        current = current->next;
      }
    }
    fNumParticles--;
  }
  int GetNumParticles() const { return fNumParticles; }
  Iterator begin() { return ParticleIterator{head}; }
  Iterator end() { return ParticleIterator{head + fNumParticles}; }

private:
  Particle *head;
  int fNumParticles;
};

#endif