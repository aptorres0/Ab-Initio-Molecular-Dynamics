// #include "vec.h"
// #include <algorithm> // std::copy
// #include <stdarg.h>
/*
// ---> TESTING
template <typename T> vec<T>::vec() : fDIM(0), fV(NULL) {}
// construct 0 array of size n
template <typename T>
vec<T>::vec(int n) : fDIM{n}, fV{fDIM ? new T[n] : nullptr} {}
template <typename T>
vec<T>::vec(int _DIM, T...) : fDIM{_DIM}, fV{fDIM ? new T[fDIM] : nullptr} {
  // initialize argument list variable
  va_list valist;
  va_start(valist, fDIM); // initialize valist for num number of arguments
  for (int i = 0; i < fDIM; ++i) { // access all the arguments assigned to
                                   // valist
    fV[i] = va_arg(valist, T);
  }
}
// <--- TESTING
template <typename T> vec<T>::vec(T x, T y) {
  fV[0] = x;
  fV[1] = y;
}
// 3D vector
template <typename T> vec<T>::vec(T x, T y, T z) : fDIM{3}, fV{new T[3]} {
  fV[0] = x;
  fV[1] = y;
  fV[2] = z;
}
// copy-constructor
template <typename T>
vec<T>::vec(const vec<T> &other)
    : fDIM{other.fDIM}, fV{fDIM ? new T[fDIM] : nullptr} {
  std::copy(other.fV, other.fV + fDIM, fV);
}
// move constructor
template <typename T>
vec<T>::vec(vec<T> &&other) noexcept
    : vec() // initialize via default constructor
{
  swap(*this, other);
}
template <typename T> void swap(vec<T> &first, vec<T> &second) { // nothrow
  using std::swap;
  swap(first.fDIM, second.fDIM);
  swap(first.fV, second.fV);
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
* /