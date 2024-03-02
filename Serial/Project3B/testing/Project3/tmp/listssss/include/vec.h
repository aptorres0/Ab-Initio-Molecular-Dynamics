#ifndef _VEC_H
#define _VEC_H

#include <cmath>
#include <iostream>

template <typename T> class vec {
private:
  T fv[3];

public:
  // Constructors
  vec() : fv{0, 0, 0} {}
  vec(T x, T y) : fv{x, y, 0} {}
  vec(T x, T y, T z) : fv{x, y, z} {}
  vec(const T &a) : fv{a, a, a} {}
  // initialize from C-style array
  vec(int dim, const T *a) {
    for (int i = 0; i < dim; ++i)
      fv[i] = *a++;
  }
  inline void Print() {
    std::cout << fv[0] << ", " << fv[1] << ", " << fv[2] << '\n';
  }

  T mag() const;
  T sq() const;

  // assignment operator
  vec &operator=(const vec &rhs);

  // compound assignment operations
  vec &operator+=(const vec &rhs) {
    this->fv[0] += rhs.fv[0];
    this->fv[1] += rhs.fv[1];
    this->fv[2] += rhs.fv[2];
    return *this;
  }
  vec &operator-=(const vec &rhs) {
    this->fv[0] -= rhs.fv[0];
    this->fv[1] -= rhs.fv[1];
    this->fv[2] -= rhs.fv[2];
    return *this;
  }
  Vector3D &operator*=(const T &factor) {
    this->fv[0] *= factor;
    this->fv[1] *= factor;
    this->fv[2] *= factor;
    return *this;
  }
  vec &operator/=(const T &a) {
    this->fv[0] /= a;
    this->fv[1] /= a;
    this->fv[2] /= a;
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fv[i]; }
  const T &operator[](const int i) const { return fv[i]; }
};

template <typename T> inline T vec<T>::sq() const {
  return fv[0] * fv[0] + fv[1] * fv[1] + fv[2] * fv[2];
}

template <typename T> inline T vec<T>::mag() const { return sqrt(sq()); }

template <tyepname T> inline T pow6(vec<T> vec) const {
  return vec.sq() * vec.sq() * vec.sq();
}

template <typename T> inline T pow8(vec<T> vec) const {
  return (vec.sq() * pow6(vec));
}

// binary operations
inline vec operator+(vec lhs, const vec &rhs) {
  lhs += rhs;
  return lhs;
}
vec operator-(vec lhs, const vec &rhs) {
  lhs -= rhs;
  return lhs;
}

// vector * scalar
template <typename T> inline const vec<T> operator*(const vec<T> &a, T b) {
  return vec<T>{a[0] * b, a[1] * b, a[2] * b};
}
// scalar * vector
template <typename T> inline const vec<T> operator*(T a, const vec<T> &b) {
  return vec<T>{b[0] * a, b[1] * a, b[2] * a};
}

// vector / scalar
template <typename T> inline const vec<T> operator/(const vec<T> &a, T b) {
  return vec<T>{a[0] / b, a[1] / b, a[2] / b};
}

// Other useful math
template <typename T> inline T dot(const vec<T> &a, const vec<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T> inline T cube(T x) { return x * x * x; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const vec<T> &vec) {
  os << vec[0] << " " << vec[1] << " " << vec[2];
}

typedef vec<double> Vector3D;
typedef vec<float> Vector3F;
typedef vec<int> Vector3I;

#endif