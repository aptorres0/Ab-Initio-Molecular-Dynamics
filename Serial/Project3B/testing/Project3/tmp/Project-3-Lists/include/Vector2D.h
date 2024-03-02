#ifndef _VECTOR2D_H
#define _VECTOR2D_H

#include <cmath>
#include <iostream>

class Vector2D {
protected:
  double fX, fY;

public:
  // Constructors
  Vector2D() : fX{}, fY{} {}
  Vector2D(double x, double y) : fX{x}, fY{y} {}
  inline void Print() { std::cout << fX << ", " << fY << '\n'; }
  double mag() const;
  double magsq() const;
  double sq() const;

  Vector2D &operator+=(const Vector2D &rhs) {
    this->fX += rhs.fX;
    this->fY += rhs.fY;
    return *this;
  }
  Vector2D &operator-=(const Vector2D &rhs) {
    this->fX -= rhs.fX;
    this->fY -= rhs.fY;
    return *this;
  }
  Vector2D &operator*=(const double factor) {
    this->fX *= factor;
    this->fY *= factor;
    return *this;
  }
  double &operator[](const int i) {
    if (i == 0)
      return this->fX;
    else if (i == 1)
      return this->fY;
  }

  // Getters
  inline double GetX() const { return fX; }
  inline double GetY() const { return fY; }

  // Setters
  inline void SetX(double _x) { fX = _x; }
  inline void SetY(double _y) { fY = _y; }
  inline void SetXY(double _x, double _y) {
    fX = _x;
    fY = _y;
  }
};

inline double Vector2D::magsq() const { return fX * fX + fY * fY; }
inline double Vector2D::sq() const { return fX * fX + fY * fY; }
inline double Vector2D::mag() const { return sqrt(magsq()); }

inline const Vector2D operator+(const Vector2D &a, const Vector2D &b) {
  return Vector2D{a.GetX() + b.GetX(), a.GetY() + b.GetY()};
}

inline const Vector2D operator-(const Vector2D &a, const Vector2D &b) {
  return Vector3D{a.GetX() - b.GetX(), a.GetY() - b.GetY()};
}

// Vector * scalar
inline const Vector2D operator*(const Vector2D &a, double b) {
  return Vector2D{a.GetX() * b, a.GetY() * b};
}

// scalar * Vector
inline const Vector2D operator*(double a, const Vector2D &b) {
  return Vector3D{a * b.GetX(), a * b.GetY()};
}

// Vector / scalar
inline const Vector2D operator/(const Vector2D &a, double b) {
  return Vector2D{a.GetX() / b, a.GetY() / b};
}

// Other useful math
inline double dot(const Vector2D &a, const Vector2D &b) {
  return a.GetX() * b.GetX() + a.GetY() * b.GetY();
}

#endif