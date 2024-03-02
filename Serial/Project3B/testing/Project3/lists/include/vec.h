/**
 * vec.h
 *
 * Author: Alexander Torres
 * Date: 21 MAR 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Generic static vector class for holding numerical arrays and performing
 *  arithmetic. Code is adapted for Press et al. NRvector class and from the
 *  following stackoverflow.com posts:
 *  https://stackoverflow.com/questions/4421706/what-are-the-basic-rules-and-idioms-for-operator-overloading
 *  https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
 *
 */

#ifndef _VEC_H
#define _VEC_H

#include "MyRNG.h"
#include <algorithm> // std::copy
#include <cmath>
#include <cstdarg>
#include <iostream>

template <typename T> class vec {
private:
  int fDIM;
  T *fV;

public:
  // default constructor
  vec() : fDIM{0}, fV{NULL} {}
  // construct 0 array of size n
  explicit vec(int n) : fDIM{n}, fV{fDIM ? new T[n] : nullptr} {}
  // initialize to constant value
  vec(int n, const T &a);
  // initialize from array
  vec(int n, const T *a);
  // variable length vector
  vec(int _DIM, T...);
  // 2D vector
  vec(T x, T y);
  // 3D vector
  vec(T x, T y, T z);
  // copy-constructor
  vec(const vec &other) : fDIM{other.fDIM}, fV{fDIM ? new T[fDIM] : nullptr} {
    std::copy(other.fV, other.fV + fDIM, fV);
  }
  // move constructor
  vec(vec &&other) noexcept
      : vec() // initialize via default constructor
  {
    swap(*this, other);
  }
  friend void swap(vec &first, vec &second) { // nothrow
    using std::swap;
    swap(first.fDIM, second.fDIM);
    swap(first.fV, second.fV);
  }

  // destructor
  ~vec() { delete[] fV; }

  // assignment operator
  vec &operator=(vec other) {
    swap(*this, other);
    return *this;
  }

  T mag() const;
  T sq() const;

  // compound assignment operations
  vec &operator+=(const vec &rhs) {
    for (int i = 0; i < fDIM; ++i) {
      this->fV[i] += rhs.fV[i];
    }
    return *this;
  }
  vec &operator-=(const vec &rhs) {
    for (int i = 0; i < fDIM; ++i) {
      this->fV[i] -= rhs.fV[i];
    }
    return *this;
  }
  vec &operator*=(const T &factor) {
    for (int i = 0; i < fDIM; ++i) {
      this->fV[i] *= factor;
    }
    return *this;
  }
  vec &operator/=(const T &factor) {
    for (int i = 0; i < fDIM; ++i) {
      this->fV[i] /= factor;
    }
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fV[i]; }
  inline const T &operator[](const int i) const { return fV[i]; }

  // assign values
  void assign(int _DIM, const T &a) {
    for (int i = 0; i < _DIM; ++i) {
      this->fV[i] = a;
    }
  }
  void assign(int _DIM, const T *a) {
    for (int i = 0; i < _DIM; ++i) {
      this->fV[i] = *a++;
    }
  }
  void assign_rnd(int _DIM) {
    init_random();
    for (int i = 0; i < _DIM; ++i) {
      this->fV[i] = rnd();
    }
    save_random();
  }

  inline int length() const { return this->fDIM; }
};

template <typename T>
vec<T>::vec(int _DIM, T...) : fDIM{_DIM}, fV{fDIM ? new T[fDIM] : nullptr} {
  va_list valist;
  va_start(valist, fDIM);
  for (int i = 0; i < fDIM; ++i) {
    fV[i] = va_arg(valist, T);
  }
  va_end(valist);
}

template <typename T> vec<T>::vec(T x, T y) : fDIM{2}, fV{new T[fDIM]} {
  fV[0] = x;
  fV[1] = y;
}
// 3D vector
template <typename T> vec<T>::vec(T x, T y, T z) : fDIM{3}, fV{new T[3]} {
  fV[0] = x;
  fV[1] = y;
  fV[2] = z;
}
template <typename T>
vec<T>::vec(int n, const T &a) : fDIM{n}, fV{fDIM ? new T[fDIM] : nullptr} {
  for (int i = 0; i < fDIM; ++i)
    fV[i] = a;
}
template <typename T>
vec<T>::vec(int n, const T *a) : fDIM{n}, fV{fDIM ? new T[fDIM] : nullptr} {
  for (int i = 0; i < fDIM; ++i)
    fV[i] = *a++;
}

template <typename T> inline T vec<T>::sq() const {
  T result = 0;
  for (int i = 0; i < this->fDIM; ++i) {
    result += fV[i] * fV[i];
  }
  return result;
}

template <typename T> inline T vec<T>::mag() const { return sqrt(sq()); }

template <typename T> inline T pow6(const vec<T> &vec) {
  T val = vec.sq();
  return val * val * val;
}

template <typename T> inline T pow8(const vec<T> &vec) {
  return (vec.sq() * pow6(vec));
}

// binary operations
// NOTE: the `+`, `-`, `/`, and `*` operators are returning copies of lhs, so
// lhs is not actually altered
template <typename T> inline vec<T> operator+(vec<T> lhs, const vec<T> &rhs) {
  lhs += rhs;
  return lhs;
}
template <typename T> inline vec<T> operator-(vec<T> lhs, const vec<T> &rhs) {
  lhs -= rhs;
  return lhs;
}
// vector * scalar
template <typename T>
inline const vec<T> operator*(vec<T> lhs, const double factor) {
  lhs *= factor;
  return lhs;
}
// scalar * vector
template <typename T>
inline const vec<T> operator*(const double factor, vec<T> rhs) {
  rhs *= factor;
  return rhs;
}
// vector / scalar
template <typename T>
inline const vec<T> operator/(vec<T> lhs, const double factor) {
  lhs /= factor;
  return lhs;
}

// Other useful math
template <typename T> inline T dot(const vec<T> &a, const vec<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T> inline T cube(T x) { return x * x * x; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const vec<T> &v) {
  for (int i = 0; i < v.length(); ++i)
    os << v[i] << " ";
  return os;
}

// For checking boundary conditions:
template <typename T>
inline bool operator<(const vec<T> &lhs, const vec<T> &rhs) {
  for (int i = 0; i < lhs.length(); ++i) {
    if (lhs[i] < rhs[i])
      return true;
  }
  return false;
}
template <typename T>
inline bool operator>(const vec<T> &lhs, const vec<T> &rhs) {
  return operator<(rhs, lhs);
}
template <typename T>
inline bool operator<(const vec<T> &lhs, const double rhs) {
  for (int i = 0; i < lhs.length(); ++i) {
    if (lhs[i] < rhs)
      return true;
  }
  return false;
}
template <typename T>
inline bool operator>(const double lhs, const vec<T> &rhs) {
  return operator<(rhs, lhs);
}

typedef vec<double> VectorD;
typedef vec<float> VectorF;
typedef vec<int> VectorI;

#endif