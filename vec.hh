#ifndef VEC_HH_
#define VEC_HH_ 1

#include <istream>
#include <ostream>

template<int N, typename T>
class vec {
 private:
  T array[N];

 public:
  vec() {}
  vec(const T &u)   { for (int n = 0; n < N; ++n) array[n] = u; }
  vec(const vec &v) { for (int n = 0; n < N; ++n) array[n] = v[n]; }
  vec(const T p[N]) { for (int n = 0; n < N; ++n) array[n] = p[n]; }
  ~vec() {}

  static const int size = N;
  typedef T value_type;

  template<typename U>
  operator vec<N, U>() const {
    vec<N, U> v;
    for (int n = 0; n < N; ++n) v[n] = static_cast<U>(array[n]);
    return v;
  }
  operator T*()             { return array; }
  operator const T*() const { return array; }

  T &operator[](int n)             { return array[n]; }
  const T &operator[](int n) const { return array[n]; }

#define VEC_HH__DEFINE_ASSIGN_OP(op) \
  vec &operator op(const T &u) { \
    for (int n = 0; n < N; ++n) array[n] op u; \
    return *this; \
  } \
  vec &operator op(const vec &v) { \
    for (int n = 0; n < N; ++n) array[n] op v[n]; \
    return *this; \
  }
  VEC_HH__DEFINE_ASSIGN_OP(=)
  VEC_HH__DEFINE_ASSIGN_OP(+=)
  VEC_HH__DEFINE_ASSIGN_OP(-=)
  VEC_HH__DEFINE_ASSIGN_OP(*=)
  VEC_HH__DEFINE_ASSIGN_OP(/=)
  VEC_HH__DEFINE_ASSIGN_OP(%=)

  VEC_HH__DEFINE_ASSIGN_OP(^=)
  VEC_HH__DEFINE_ASSIGN_OP(&=)
  VEC_HH__DEFINE_ASSIGN_OP(|=)
  VEC_HH__DEFINE_ASSIGN_OP(>>=)
  VEC_HH__DEFINE_ASSIGN_OP(<<=)
#undef VEC_HH__DEFINE_ASSIGN_OP

#define VEC_HH__DEFINE_BINARY_OP(op) \
  friend vec operator op(const vec &x, const vec &y) { \
    vec v; \
    for (int n = 0; n < N; ++n) v[n] = (x[n] op y[n]); \
    return v; \
  } \
  friend vec operator op(const T &u, const vec &y) { \
    vec v; \
    for (int n = 0; n < N; ++n) v[n] = (u op y[n]); \
    return v; \
  } \
  friend vec operator op(const vec &x, const T &u) { \
    vec v; \
    for (int n = 0; n < N; ++n) v[n] = (x[n] op u); \
    return v; \
  }
  VEC_HH__DEFINE_BINARY_OP(<)
  VEC_HH__DEFINE_BINARY_OP(>)
  VEC_HH__DEFINE_BINARY_OP(<=)
  VEC_HH__DEFINE_BINARY_OP(>=)
  VEC_HH__DEFINE_BINARY_OP(==)
  VEC_HH__DEFINE_BINARY_OP(!=)

  VEC_HH__DEFINE_BINARY_OP(+)
  VEC_HH__DEFINE_BINARY_OP(-)
  VEC_HH__DEFINE_BINARY_OP(*)
  VEC_HH__DEFINE_BINARY_OP(/)
  VEC_HH__DEFINE_BINARY_OP(%)

  VEC_HH__DEFINE_BINARY_OP(^)
  VEC_HH__DEFINE_BINARY_OP(&)
  VEC_HH__DEFINE_BINARY_OP(|)
  VEC_HH__DEFINE_BINARY_OP(>>)
  VEC_HH__DEFINE_BINARY_OP(<<)

  VEC_HH__DEFINE_BINARY_OP(&&)
  VEC_HH__DEFINE_BINARY_OP(||)
#undef VEC_HH__DEFINE_BINARY_OP

#define VEC_HH__DEFINE_UNARY_OP(op) \
  vec operator op() const { \
    vec v; \
    for (int n = 0; n < N; ++n) v[n] = op array[n]; \
    return v; \
  }
  VEC_HH__DEFINE_UNARY_OP(~)
  VEC_HH__DEFINE_UNARY_OP(!)
  VEC_HH__DEFINE_UNARY_OP(-)
  VEC_HH__DEFINE_UNARY_OP(+)
#undef VEC_HH__DEFINE_UNARY_OP

#define VEC_HH__DEFINE_PREFIX_OP(op) \
  vec &operator op() { \
    for (int n = 0; n < N; ++n) op array[n]; \
    return *this; \
  }
  VEC_HH__DEFINE_PREFIX_OP(++)
  VEC_HH__DEFINE_PREFIX_OP(--)
#undef VEC_HH__DEFINE_PREFIX_OP

#define VEC_HH__DEFINE_POSTFIX_OP(op) \
  vec operator op(int) { \
    vec v; \
    for (int n = 0; n < N; ++n) v[n] = array[n] op; \
    return v; \
  }
  VEC_HH__DEFINE_POSTFIX_OP(++)
  VEC_HH__DEFINE_POSTFIX_OP(--)
#undef VEC_HH__DEFINE_POSTFIX_OP

  T min() const {
    T min = array[0];
    for (int n = 1; n < N; ++n) { if (array[n] < min) min = array[n]; }
    return min;
  }
  T max() const {
    T max = array[0];
    for (int n = 1; n < N; ++n) { if (array[n] > max) max = array[n]; }
    return max;
  }
  T sum() const {
    T sum = array[0];
    for (int n = 1; n < N; ++n) { sum += array[n]; }
    return sum;
  }
  T mul() const {
    T mul = array[0];
    for (int n = 1; n < N; ++n) { mul *= array[n]; }
    return mul;
  }
  T norm2() const {
    T norm2 = array[0] * array[0];
    for (int n = 1; n < N; ++n) { norm2 += array[n] * array[n]; }
    return norm2;
  }

  friend T sprod(const vec &x, const vec &y) {
    T prod = x[0] * y[0];
    for (int n = 1; n < N; ++n) { prod += x[n] * y[n]; }
    return prod;
  }

  friend std::ostream &operator<<(std::ostream &os, const vec &v) {
    os << v[0];
    for (int n = 1; n < N; ++n) { os << " " << v[n]; }
    return os;
  }
  friend std::istream &operator>>(std::istream &is, vec &v) {
    is >> v[0];
    for (int n = 1; n < N; ++n) { is.ignore(); is >> v[n]; }
    return is;
  }
};

#endif  // VEC_HH_
