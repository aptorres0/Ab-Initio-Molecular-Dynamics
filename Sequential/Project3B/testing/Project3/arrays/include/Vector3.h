#ifndef _VECTOR3_H
#define _VECTOR3_H

#include <cmath>
#include <iostream>

template <typename T> class Vector3 {
private:
  T fv[3];

public:
  // Constructors
  Vector3() : fv{} {}
  Vector3(T x, T y) : fv{x, y, 0} {}
  Vector3(T x, T y, T z) : fv{x, y, z} {}
  Vector3(const T &a) : fv{a, a, a} {}
  // initialize from C-style array
  Vector3(const T *a) {
    for (int i = 0; i < 3; ++i)
      fv[i] = *a++;
  }
  inline void Print() {
    std::cout << fv[0] << ", " << fv[1] << ", " << fv[2] << '\n';
  }

  T mag() const;
  T sq() const;

  // compound assignment operations
  Vector3 &operator+=(const Vector3 &rhs) {
    this->fv[0] += rhs.fv[0];
    this->fv[1] += rhs.fv[1];
    this->fv[2] += rhs.fv[2];
    return *this;
  }
  Vector3 &operator-=(const Vector3 &rhs) {
    this->fv[0] -= rhs.fv[0];
    this->fv[1] -= rhs.fv[1];
    this->fv[2] -= rhs.fv[2];
    return *this;
  }
  Vector3tor3D &operator*=(const T &factor) {
    this->fv[0] *= factor;
    this->fv[1] *= factor;
    this->fv[2] *= factor;
    return *this;
  }
  Vector3 &operator/=(const T &a) {
    this->fv[0] /= a;
    this->fv[1] /= a;
    this->fv[2] /= a;
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fv[i]; }
  const T &operator[](const int i) const { return fv[i]; }
};

template <typename T> inline T Vector3<T>::sq() const {
  return fv[0] * fv[0] + fv[1] * fv[1] + fv[2] * fv[2];
}

template <typename T> inline T Vector3<T>::mag() const { return sqrt(sq()); }

template <tyepname T> inline T pow6(Vector3<T> vec) const {
  return vec.sq() * vec.sq() * vec.sq();
}

template <typename T> inline T pow8(Vector3<T> vec) const {
  return (vec.sq() * pow6(vec));
}

// binary operations
inline Vector3 operator+(Vector3 lhs, const Vector3 &rhs) {
  lhs += rhs;
  return lhs;
}
Vector3 operator-(Vector3 lhs, const Vector3 &rhs) {
  lhs -= rhs;
  return lhs;
}

// vector * scalar
template <typename T>
inline const Vector3<T> operator*(const Vector3<T> &a, T b) {
  return Vector3<T>{a[0] * b, a[1] * b, a[2] * b};
}
// scalar * vector
template <typename T>
inline const Vector3<T> operator*(T a, const Vector3<T> &b) {
  return Vector3<T>{b[0] * a, b[1] * a, b[2] * a};
}

// vector / scalar
template <typename T>
inline const Vector3<T> operator/(const Vector3<T> &a, T b) {
  return Vector3<T>{a[0] / b, a[1] / b, a[2] / b};
}

// Other useful math
template <typename T> inline T dot(const Vector3<T> &a, const Vector3<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T> inline T cube(T x) { return x * x * x; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector3<T> &vec) {
  os << vec[0] << " " << vec[1] << " " << vec[2];
}

typedef Vector3<double> Vector3D;
typedef Vector3<float> Vector3F;
typedef Vector3<int> Vector3I;

#endif