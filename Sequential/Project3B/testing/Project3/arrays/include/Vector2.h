#ifndef _VECTOR2_H
#define _VECTOR2_H

#include <cmath>
#include <iostream>

template <typename T> class Vector2 {
private:
  T fv[2];

public:
  // Constructors
  Vector2() : fv{} {}
  Vector2(T x, T y) : fv{x, y} {}
  Vector2(const T &a) : fv{a, a} {}
  // initialize from C-style array
  Vector2(const T *a) {
    for (int i = 0; i < 2; ++i)
      fv[i] = *a++;
  }
  inline void Print() { std::cout << fv[0] << ", " << fv[1] << ", " << '\n'; }

  T mag() const;
  T sq() const;

  // compound assignment operations
  Vector2 &operator+=(const Vector3 &rhs) {
    this->fv[0] += rhs.fv[0];
    this->fv[1] += rhs.fv[1];
    return *this;
  }
  Vector2 &operator-=(const Vector3 &rhs) {
    this->fv[0] -= rhs.fv[0];
    this->fv[1] -= rhs.fv[1];
    return *this;
  }
  Vector2tor3D &operator*=(const T &factor) {
    this->fv[0] *= factor;
    this->fv[1] *= factor;
    return *this;
  }
  Vector2 &operator/=(const T &a) {
    this->fv[0] /= a;
    this->fv[1] /= a;
    return *this;
  }

  // subscript operations
  inline T &operator[](const int i) { return fv[i]; }
  const T &operator[](const int i) const { return fv[i]; }
};

template <typename T> inline T Vector2<T>::sq() const {
  return fv[0] * fv[0] + fv[1] * fv[1];
}

template <typename T> inline T Vector2<T>::mag() const { return sqrt(sq()); }

template <tyepname T> inline T pow6(Vector2<T> vec) const {
  return vec.sq() * vec.sq() * vec.sq();
}

template <typename T> inline T pow8(Vector2<T> vec) const {
  return (vec.sq() * pow6(vec));
}

// binary operations
inline Vector2 operator+(Vector3 lhs, const Vector3 &rhs) {
  lhs += rhs;
  return lhs;
}
Vector2 operator-(Vector3 lhs, const Vector3 &rhs) {
  lhs -= rhs;
  return lhs;
}

// vector * scalar
template <typename T>
inline const Vector2<T> operator*(const Vector3<T> &a, T b) {
  return Vector2<T>{a[0] * b, a[1] * b};
}
// scalar * vector
template <typename T>
inline const Vector2<T> operator*(T a, const Vector3<T> &b) {
  return Vector2<T>{b[0] * a, b[1] * a};
}

// vector / scalar
template <typename T>
inline const Vector2<T> operator/(const Vector3<T> &a, T b) {
  return Vector2<T>{a[0] / b, a[1] / b};
}

// Other useful math
template <typename T> inline T dot(const Vector2<T> &a, const Vector3<T> &b) {
  return a[0] * b[0] + a[1] * b[1];
}

template <typename T> inline T cube(T x) { return x * x * x; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector2<T> &vec) {
  os << vec[0] << " " << vec[1];
}

typedef Vector2<double> Vector2D;
typedef Vector2<float> Vector2F;
typedef Vector2<int> Vector2I;

#endif