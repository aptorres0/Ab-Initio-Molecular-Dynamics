#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <cmath>
#include <iostream>

class Vector3D {
private:
  double fX, fY, fZ;

public:
  // Constructors
  Vector3D() : fX{}, fY{}, fZ{} {}
  Vector3D(double x, double y, double z) : fX{x}, fY{y}, fZ{z} {}
  // copy-constructor
  inline void Print() { std::cout << fX << ", " << fY << ", " << fZ << '\n'; }
  double mag() const;
  double magsq() const;

  Vector3D &operator+=(const Vector3D &rhs) {
    this->fX += rhs.fX;
    this->fY += rhs.fY;
    this->fZ += rhs.fZ;
    return *this;
  }
  Vector3D &operator-=(const Vector3D &rhs) {
    this->fX -= rhs.fX;
    this->fY -= rhs.fY;
    this->fZ -= rhs.fZ;
    return *this;
  }

  // Getters
  inline double GetX() const { return fX; }
  inline double GetY() const { return fY; }
  inline double GetZ() const { return fZ; }
};

inline double Vector3D::magsq() const { return fX * fX + fY * fY + fZ * fZ; }
inline double Vector3D::mag() const { return sqrt(magsq()); }

/*Vector3D &Vector3D::operator=(const Vector3D &rhs) {
  // check for self-assignment
  if (this == &rhs)
    return *this;
  fX = rhs.fX;
  fY = rhs.fY;
  fZ = rhs.fZ;
  return *this;
}*/

inline const Vector3D operator+(const Vector3D &a, const Vector3D &b) {
  return Vector3D{a.GetX() + b.GetX(), a.GetY() + b.GetY(),
                  a.GetZ() + b.GetZ()};
}

inline const Vector3D operator-(const Vector3D &a, const Vector3D &b) {
  return Vector3D{a.GetX() - b.GetX(), a.GetY() - b.GetY(),
                  a.GetZ() - b.GetZ()};
}

// Vector * scalar
inline const Vector3D operator*(const Vector3D &a, double b) {
  return Vector3D{a.GetX() * b, a.GetY() * b, a.GetZ() * b};
}

// scalar * Vector
inline const Vector3D operator*(double a, const Vector3D &b) {
  return Vector3D{a * b.GetX(), a * b.GetY(), a * b.GetZ()};
}

// Other useful math
inline double dot(const Vector3D &a, const Vector3D &b) {
  return a.GetX() * b.GetX() + a.GetY() * b.GetY() + a.GetZ() * b.GetZ();
}

inline Vector3D cross(const Vector3D &a, const Vector3D &b) {
  return Vector3D{a.GetY() * b.GetZ() - a.GetZ() * b.GetY(),
                  a.GetZ() * b.GetX() - a.GetX() * b.GetZ(),
                  a.GetX() * b.GetY() - a.GetY() * b.GetX()};
}

#endif