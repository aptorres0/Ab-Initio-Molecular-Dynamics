#ifndef _VEC_H
#define _VEC_H

#include "MyRNG.h"
#include <algorithm> // std::copy
#include <cmath>
#include <cstdarg>
#include <iostream>

template <typename T, int DIM> class vec {
private:
  T fV[DIM];

public:
  // default constructor
  vec() {
    for (int i = 0; i < DIM; ++i) {
      fV[i] = 0;
    }
  }
  vec(T x, T y) {
    fV[0] = x;
    fV[1] = y;
  }
  vec(T x, T y, T z) {
    fV[0] = x;
    fV[1] = y;
    fV[2] = z;
  }
  // initialize to constant value
  vec(const T &a) {
    for (int i = 0; i < DIM; ++i) {
      fV[i] = a;
    }
  }
  // initialize from array
  vec(const T *a) {
    for (int i = 0; i < DIM; ++i) {
      fV[i] = *a++;
    }
  }
  // copy-constructor
  vec(const vec &other) { std::copy(other.fV, other.fV + DIM, fV); }
  // move constructor
  // NOTE: the move constructor handles temporary objects, such as the result of
  // a function call.
  //  "move" a temporary object into an existing object, rather than copying it.
  vec(vec &&other) noexcept
      : vec() // initialize via default constructor
  {
    swap(*this, other);
  }
  friend void swap(vec &first, vec &second) { // nothrow
    using std::swap;
    swap(first.fV, second.fV);
  }

  // assignment operator
  vec &operator=(vec other) {
    swap(*this, other);
    return *this;
  }

  T mag() const;
  T sq() const;

  // compound assignment operations
  vec &operator+=(const vec &rhs) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] += rhs.fV[i];
    }
    return *this;
  }
  vec &operator-=(const vec &rhs) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] -= rhs.fV[i];
    }
    return *this;
  }
  vec &operator*=(const T &factor) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] *= factor;
    }
    return *this;
  }
  vec &operator/=(const T &factor) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] /= factor;
    }
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fV[i]; }
  inline const T &operator[](const int i) const { return fV[i]; }

  // assign values
  void assign(const T &a) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] = a;
    }
  }
  void assign(const T *a) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] = *a++;
    }
  }
  void set_rnd() {
    init_random();
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] = rnd();
    }
    save_random();
  }
  void assign_to(const vec &source) {
    for (int i = 0; i < DIM; ++i) {
      this->fV[i] = source.fV[i];
    }
  }

  constexpr int Size() const { return DIM; }
  inline int length() const { return DIM; }
};

template <typename T, int DIM> inline T vec<T, DIM>::sq() const {
  T result = 0;
  for (int i = 0; i < DIM; ++i) {
    result += fV[i] * fV[i];
  }
  return result;
}

template <typename T, int DIM> inline T vec<T, DIM>::mag() const {
  return sqrt(sq());
}

template <typename T, int DIM> inline T pow6(const vec<T, DIM> &vec) {
  T val = vec.sq();
  return val * val * val;
}

template <typename T, int DIM> inline T pow8(const vec<T, DIM> &vec) {
  return (vec.sq() * pow6(vec));
}

// binary operations
// NOTE: the `+`, `-`, `/`, and `*` operators are returning copies of lhs, so
// lhs is not actually altered
template <typename T, int DIM>
inline vec<T, DIM> operator+(vec<T, DIM> lhs, const vec<T, DIM> &rhs) {
  lhs += rhs;
  return lhs;
}
template <typename T, int DIM>
inline vec<T, DIM> operator-(vec<T, DIM> lhs, const vec<T, DIM> &rhs) {
  lhs -= rhs;
  return lhs;
}
// vector * scalar
template <typename T, int DIM>
inline const vec<T, DIM> operator*(vec<T, DIM> lhs, const double factor) {
  lhs *= factor;
  return lhs;
}
// scalar * vector
template <typename T, int DIM>
inline const vec<T, DIM> operator*(const double factor, vec<T, DIM> rhs) {
  rhs *= factor;
  return rhs;
}
// vector / scalar
template <typename T, int DIM>
inline const vec<T, DIM> operator/(vec<T, DIM> lhs, const double factor) {
  lhs /= factor;
  return lhs;
}

// Other useful math
template <typename T, int DIM>
inline T dot(const vec<T, DIM> &a, const vec<T, DIM> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T, int DIM>
std::ostream &operator<<(std::ostream &os, const vec<T, DIM> &v) {
  for (int i = 0; i < v.Size(); ++i)
    os << v[i] << " ";
  return os;
}

// typedef vec<double, 3> Vector3D;

#endif