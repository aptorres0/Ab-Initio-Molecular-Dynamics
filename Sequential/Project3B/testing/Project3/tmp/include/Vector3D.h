#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <cmath>
#include <iostream>

class Vector3D {
protected:
  double fX, fY, fZ;

public:
  // Constructors
  Vector3D() : fX{}, fY{}, fZ{} {}
  Vector3D(double x, double y, double z) : fX{x}, fY{y}, fZ{z} {}
  inline void Print() { std::cout << fX << ", " << fY << ", " << fZ << '\n'; }
  double mag() const;
  double magsq() const;
  double sq() const;
  double pow3() const;
  double pow4() const;
  double pow6() const;
  double pow7() const;

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
  Vector3D &operator*=(const double factor) {
    this->fX *= factor;
    this->fY *= factor;
    this->fZ *= factor;
    return *this;
  }
  double &operator[](const int i) {
    if (i == 0)
      return this->fX;
    else if (i == 1)
      return this->fY;
    else
      return this->fZ;
  }

  // Getters
  inline double GetX() const { return fX; }
  inline double GetY() const { return fY; }
  inline double GetZ() const { return fZ; }
  inline double *GetXYZ() const {
    static double res[3] = {fX, fY, fZ};
    return res;
  }

  // Setters
  inline void SetX(double _x) { fX = _x; }
  inline void SetY(double _y) { fY = _y; }
  inline void SetZ(double _z) { fZ = _z; }
  inline void SetXYZ(double _x, double _y, double _z) {
    fX = _x;
    fY = _y;
    fZ = _z;
  }
  inline void SetXYZ(double xyz[3]) {
    fX = xyz[0];
    fY = xyz[1];
    fZ = xyz[2];
  }
};

inline double Vector3D::magsq() const { return fX * fX + fY * fY + fZ * fZ; }
inline double Vector3D::sq() const { return fX * fX + fY * fY + fZ * fZ; }
inline double Vector3D::mag() const { return sqrt(magsq()); }
inline double Vector3D::pow3() const { return sq() * mag(); }
inline double Vector3D::pow4() const { return sq() * sq(); }
inline double Vector3D::pow6() const { return sq() * sq() * sq(); }
inline double Vector3D::pow7() const { return sq() * sq() * sq() * mag(); }

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

// Vector / scalar
inline const Vector3D operator/(const Vector3D &a, double b) {
  return Vector3D{a.GetX() / b, a.GetY() / b, a.GetZ() / b};
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