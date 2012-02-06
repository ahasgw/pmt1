#ifndef VEC_HH_
#define VEC_HH_ 1

#include <iostream>
#include <istream>
#include <ostream>

template<int N, typename T>
class vec {
 private:
  T array[N];

 public:
  vec() {}
  vec(const T &u)   { for (int i = 0; i < N; ++i) array[i] = u; }
  vec(const vec &v) { for (int i = 0; i < N; ++i) array[i] = v[i]; }
  vec(const T p[N]) { for (int i = 0; i < N; ++i) array[i] = p[i]; }
  ~vec() {}

  static const int size = N;
  typedef T value_type;

  template<typename U>
  operator vec<N,U>() {
    vec<N,U> v;
    for (int i = 0; i < N; ++i) v[i] = static_cast<U>(array[i]);
    return v;
  }
  operator T*()             { return array; }
  operator const T*() const { return array; }

  T &operator[](int i)             { return array[i]; }
  const T &operator[](int i) const { return array[i]; }

#define VEC_HH__DEFINE_COMPARISON_OP(op1, op2) \
  friend bool operator op1(const vec &x, const vec &y) { \
    for (int i = 0; i < N; ++i) if (!(x[i] op1 y[i])) return false; \
    return true; \
  } \
  friend bool operator op2(const vec &x, const vec &y) { \
    for (int i = 0; i < N; ++i) if (!(x[i] op2 y[i])) return false; \
    return true; \
  }
  VEC_HH__DEFINE_COMPARISON_OP(<, >)
  VEC_HH__DEFINE_COMPARISON_OP(==, !=)
  VEC_HH__DEFINE_COMPARISON_OP(<=, >=)
#undef VEC_HH__DEFINE_COMPARISON_OP

#define VEC_HH__DEFINE_ASSIGN_OP(op) \
  vec &operator op(const T &u) { \
    for (int i = 0; i < N; ++i) array[i] op u; \
    return *this; \
  } \
  vec &operator op(const vec &v) { \
    for (int i = 0; i < N; ++i) array[i] op v[i]; \
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
    for (int i = 0; i < N; ++i) v[i] = (x[i] op y[i]); \
    return v; \
  } \
  friend vec operator op(const T &u, const vec &y) { \
    vec v; \
    for (int i = 0; i < N; ++i) v[i] = (u op y[i]); \
    return v; \
  } \
  friend vec operator op(const vec &x, const T &u) { \
    vec v; \
    for (int i = 0; i < N; ++i) v[i] = (x[i] op u); \
    return v; \
  }
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

  vec operator~() const {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = ~array[i];
    return v;
  }
  vec operator!() const {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = !array[i];
    return v;
  }

  vec operator-() const {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = -array[i];
    return v;
  }
  vec operator+() const {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = array[i];
    return v;
  }

  vec &operator++() {
    for (int i = 0; i < N; ++i) ++array[i];
    return *this;
  }
  vec &operator--() {
    for (int i = 0; i < N; ++i) --array[i];
    return *this;
  }

  vec operator++(int) {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = array[i]++;
    return v;
  }
  vec operator--(int) {
    vec v;
    for (int i = 0; i < N; ++i) v[i] = array[i]--;
    return v;
  }

  friend std::ostream &operator<<(std::ostream &os, const vec &v) {
    if (N > 0) {
      os << v[0];
      for (int i = 1; i < N; ++i) { os << " " << v[i]; }
    }
    return os;
  }
  friend std::istream &operator>>(std::istream &is, vec &v) {
    if (N > 0) {
      is >> v[0];
      for (int i = 1; i < N; ++i) { is.ignore(); is >> v[i]; }
    }
    return is;
  }
};

typedef vec<3,int> v3i;
typedef vec<3,double> v3d;

typedef vec<4,int> v4i;
typedef vec<4,double> v4d;

#endif  // VEC_HH_
