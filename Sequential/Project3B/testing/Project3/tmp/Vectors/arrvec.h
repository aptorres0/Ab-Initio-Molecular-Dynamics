namespace mine {
template <size_t nDim, typename T> struct vec {
private:
  T _elem[nDim];

public:
  vec();
  vec(const T &a); // initilize to value a
  vec(const T *a); // initilize from c-style array
};
} // namespace mine