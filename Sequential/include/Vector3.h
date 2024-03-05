#ifndef _VECTRO3_H
#define _VECTRO3_H

#include "MyRNG.h"
#include <cmath>
#include <cstdarg>
#include <cstring> // memcpy
#include <iostream>

template <typename T> class Vector3 {
private:
  T fV[3];

public:
  // default constructor
  Vector3() {
    /*fV[0] = 0.0;
    fV[1] = 0.0;
    fV[2] = 0.0;*/
  }
  Vector3(T x, T y, T z) {
    fV[0] = x;
    fV[1] = y;
    fV[2] = z;
  }
  // initialize to constant value
  Vector3(const T &a) {
    fV[0] = a;
    fV[1] = a;
    fV[2] = a;
  }
  // copy-constructor
  Vector3(const Vector3 &other) {
    // std::copy(other.fV, other.fV + 3, fV);
    memcpy(fV, other.fV, 3 * sizeof(T));
  }
  // move constructor
  // NOTE: the move constructor handles temporary objects, such as the result of
  // a function call.
  //  "move" a temporary object into an existing object, rather than copying it.
  Vector3(Vector3 &&other) noexcept
      : Vector3() // initialize via default constructor
  {
    swap(*this, other);
  }
  friend void swap(Vector3 &first, Vector3 &second) { // nothrow
    using std::swap;
    swap(first.fV, second.fV);
  }
  // assignment operator
  Vector3 &operator=(Vector3 other) {
    swap(*this, other);
    return *this;
  }

  inline T sq() const {
    T result = 0;
    result += fV[0] * fV[0];
    result += fV[1] * fV[1];
    result += fV[2] * fV[2];
    return result;
  }

  inline void sq_inplace() {
    fV[0] *= fV[0];
    fV[1] *= fV[1];
    fV[2] *= fV[2];
  }

  inline T mag() const { return sqrt(sq()); }

  Vector3 unit_vec() const {
    T m = mag();
    return Vector3(fV[0] / m, fV[1] / m, fV[2] / m);
  }

  // compound assignment operations
  Vector3 &operator+=(const Vector3 &rhs) {
    this->fV[0] += rhs.fV[0];
    this->fV[1] += rhs.fV[1];
    this->fV[2] += rhs.fV[2];
    return *this;
  }
  Vector3 &operator-=(const Vector3 &rhs) {
    this->fV[0] -= rhs.fV[0];
    this->fV[1] -= rhs.fV[1];
    this->fV[2] -= rhs.fV[2];
    return *this;
  }
  Vector3 &operator*=(const T &factor) {
    this->fV[0] *= factor;
    this->fV[1] *= factor;
    this->fV[2] *= factor;
    return *this;
  }
  Vector3 &operator/=(const T &factor) {
    this->fV[0] /= factor;
    this->fV[1] /= factor;
    this->fV[2] /= factor;
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fV[i]; }
  inline const T &operator[](const int i) const { return fV[i]; }
  inline T &x() { return fV[0]; }
  inline T &y() { return fV[1]; }
  inline T &z() { return fV[2]; }
  inline const T &x() const { return fV[0]; }
  inline const T &y() const { return fV[1]; }
  inline const T &z() const { return fV[2]; }

  // assign values
  void assign(const T &a) {
    this->fV[0] = a;
    this->fV[1] = a;
    this->fV[2] = a;
  }
  void assign(const T *a) {
    this->fV[0] = *a++;
    this->fV[1] = *a++;
    this->fV[2] = *a++;
  }
  void set_random() {
    init_random();
    this->fV[0] = rnd();
    this->fV[1] = rnd();
    this->fV[2] = rnd();
    save_random();
  }
  void assign_to(const Vector3 &source) {
    this->fV[0] = source.fV[0];
    this->fV[1] = source.fV[1];
    this->fV[2] = source.fV[2];
  }
};

template <typename T> inline T pow6(const Vector3<T> &vec) {
  T val = vec.sq();
  return val * val * val;
}

template <typename T> inline T pow8(const Vector3<T> &vec) {
  T val = vec.sq();
  return val * val * val * val;
}

// binary operations
// NOTE: the `+`, `-`, `/`, and `*` operators are returning copies of lhs, so
// lhs is not actually altered
template <typename T>
inline Vector3<T> operator+(Vector3<T> lhs, const Vector3<T> &rhs) {
  lhs += rhs;
  return lhs;
}
template <typename T>
inline Vector3<T> operator-(Vector3<T> lhs, const Vector3<T> &rhs) {
  lhs -= rhs;
  return lhs;
}
// vector * scalar
template <typename T>
inline const Vector3<T> operator*(Vector3<T> lhs, const T factor) {
  lhs *= factor;
  return lhs;
}
// scalar * vector
template <typename T>
inline const Vector3<T> operator*(const T factor, Vector3<T> rhs) {
  rhs *= factor;
  return rhs;
}
// vector / scalar
template <typename T>
inline const Vector3<T> operator/(Vector3<T> lhs, const T factor) {
  lhs /= factor;
  return lhs;
}
// for checking boundary conditions
template <typename T>
inline bool operator<(const Vector3<T> &lhs, const Vector3<T> &rhs) {
  return (lhs[0] < rhs[0] && lhs[1] < rhs[1] && lhs[2] < rhs[2]);
}
template <typename T>
inline bool operator>(const Vector3<T> &lhs, const Vector3<T> &rhs) {
  return (rhs < lhs);
}

// Other useful math
template <typename T> inline T dot(const Vector3<T> &a, const Vector3<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
  os << v[0] << " " << v[1] << " " << v[2] << " ";
  return os;
}

typedef Vector3<double> Vector3D;
typedef Vector3<float> Vector3F;

#endif