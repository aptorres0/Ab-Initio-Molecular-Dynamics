#ifndef _VEC_H
#define _VEC_H

template <typename T> class vec {
private:
  size_t fDim; // Size of array, indices 0,...,nDim-1
  T *fv;       // Pointer to data array
public:
  // Constructors
  vec();                  // Default constructor
  explicit vec(int n);    // n-Dim array of 0s
  vec(int n, const T &a); // initialize to value a
  vec(int n, const T *a); // initialize from C-style array
  vec(T _X, T _Y);        // initilize to 2d array
  vec(T _X, T _Y, T _Z);  // initilize to 3d array
  //  Destructor
  ~vec();

  // assignment operator
  vec &operator=(const vec &rhs);

  // compound assignment operations
  vec &operator+=(const vec &rhs) {}
  vec &operator-=(const vec &rhs);
  vec &operator*=(const T &a);
  vec &operator/=(const T &a);

  // binary operations
  vec operator+(vec lhs, const vec &rhs);
  vec operator-(vec lhs, const vec &rhs);
  vec operator*(vec lhs, const vec &rhs);

  // subscript operations
  inline T &operator[](const int i) { return fv[i]; }
  const T &operator[](const int i) const { return fv[i]; }

  // Make T available
  typedef T value_type;

  // Return size of vector
  inline int size() const { return fDim; }
};

template <typename T>
vec<T>::vec(int n, const T &a) : fDim{n}, fv{n > 0 ? new T[n] : NULL} {
  for (int i = 0; i < n; ++i)
    fv[i] = a;
}

template <typename T>
vec<T>::vec(int n) : fDim{n}, fv{n > 0 ? new T[n] : NULL} {}

template <typename T>
vec<T>::vec(int n, const T &a) : fDim{n}, fv{n > 0 ? new T[n] : NULL} {
  for (int i = 0; i < n; ++i)
    fv[i] = a;
}

template <typename T>
vec<T>::vec(int n, const T *a) : fDim{n}, fv{n > 0 ? new T[n] : NULL} {
  for (int i = 0; i < n; ++i)
    fv[i] = *(a++);
}

template <typename T> vec<T>::vec(T _X, T _Y) : fDim{2}, fv{new T[2] : NULL} {
  fv[0] = _X;
  fv[1] = _Y;
}

template <typename T>
vec<T>::vec(T _X, T _Y, T _Z) : fDim{3}, fv{new T[3] : NULL} {
  fv[0] = _X;
  fv[1] = _Y;
  fv[2] = _Z;
}

//  Destructor
template <typename T> vec<T>::~vec() {
  if (fv != NULL)
    delete[] (fv);
}

// Assignment operator
template <typename T> vec<T> &vec<T>::operator=(const vec<T> &rhs) {
  // postcondition: normal assignment via copying has been performed;
  // if vector and rhs were different sizes, vector has been resized to match
  // the size of rhs
  if (this != &rhs) { // check for self assignment
    if (fDim != rhs.fDim) {
      if (fv != NULL)
        delete[] (fv);
      fDim = rhs.fDim;
      v = fDim > 0 ? new T[fDim] : NULL:
    }
    for (int i = 0; i < fDim; ++i)
      v[i] = rhs[i];
  }
  return *this;
}

// compound assignment operations
template <typename T>
vec<T>::

    template <typename T>
    inline const T vec<T> operator+(const vec<T> &a, const vec<T> &b) {}

// return element i
inline T &operator[](const int i) { return v[i]; } // allows for modifying value
inline const T &operator[](const int i) const { return v[i]; } // const version

#endif