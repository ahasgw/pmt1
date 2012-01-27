#ifndef VEC_HH_
#define VEC_HH_ 1

#include <iostream>

template<int N, typename T>
class vec {
 private:
  T array[N];

 public:
  vec() {}

  vec(const T &u) {
    for (int i = 0; i < N; ++i) array[i] = u;
  }

  vec(const vec &v) {
    for (int i = 0; i < N; ++i) array[i] = v[i];
  }

  ~vec() {}


  operator T*() {
    return array;
  }

  operator const T*() const {
    return array;
  }

  T &operator[](int i) {
    return array[i];
  }

  const T &operator[](int i) const {
    return array[i];
  }

  
  friend bool operator==(const vec &lhs, const vec &rhs) {
    for (int i = 0; i < N; ++i)
      if (lhs.array[i] != rhs.array[i]) return false;
    return true;
  }

  friend bool operator!=(const vec &lhs, const vec &rhs) {
    for (int i = 0; i < N; ++i)
      if (lhs.array[i] != rhs.array[i]) return true;
    return false;
  }
  

  template<typename U>
  vec &operator=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] = u;
    return *this;
  }

  template<typename U>
  vec &operator+=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] += u;
    return *this;
  }

  template<typename U>
  vec &operator-=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] -= u;
    return *this;
  }

  template<typename U>
  vec &operator*=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] *= u;
    return *this;
  }

  template<typename U>
  vec &operator/=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] /= u;
    return *this;
  }

  template<typename U>
  vec &operator%=(const U &u) {
    for (int i = 0; i < N; ++i) array[i] %= u;
    return *this;
  }

  
  template<typename U>
  vec &operator=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] = v[i];
    return *this;
  }

  template<typename U>
  vec &operator+=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] += v[i];
    return *this;
  }

  template<typename U>
  vec &operator-=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] -= v[i];
    return *this;
  }

  template<typename U>
  vec &operator*=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] *= v[i];
    return *this;
  }

  template<typename U>
  vec &operator/=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] /= v[i];
    return *this;
  }

  template<typename U>
  vec &operator%=(const vec<N,U> &v) {
    for (int i = 0; i < N; ++i) array[i] %= v[i];
    return *this;
  }


  friend vec operator+(const T &val, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = val + rhs.array[i];
    return v;
  }

  friend vec operator-(const T &val, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = val - rhs.array[i];
    return v;
  }

  friend vec operator*(const T &val, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = val * rhs.array[i];
    return v;
  }

  friend vec operator/(const T &val, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = val / rhs.array[i];
    return v;
  }

  friend vec operator%(const T &val, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = val % rhs.array[i];
    return v;
  }


  friend vec operator+(const vec &rhs, const T &val) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = rhs.array[i] + val;
    return v;
  }

  friend vec operator-(const vec &rhs, const T &val) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = rhs.array[i] - val;
    return v;
  }

  friend vec operator*(const vec &rhs, const T &val) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = rhs.array[i] * val;
    return v;
  }

  friend vec operator/(const vec &rhs, const T &val) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = rhs.array[i] / val;
    return v;
  }

  friend vec operator%(const vec &rhs, const T &val) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = rhs.array[i] % val;
    return v;
  }


  friend vec operator+(const vec &lhs, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = lhs.array[i] + rhs.array[i];
    return v;
  }

  friend vec operator-(const vec &lhs, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = lhs.array[i] - rhs.array[i];
    return v;
  }

  friend vec operator*(const vec &lhs, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = lhs.array[i] * rhs.array[i];
    return v;
  }

  friend vec operator/(const vec &lhs, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = lhs.array[i] / rhs.array[i];
    return v;
  }

  friend vec operator%(const vec &lhs, const vec &rhs) {
    vec v;
    for (int i = 0; i < N; ++i) v.array[i] = lhs.array[i] % rhs.array[i];
    return v;
  }


  friend std::ostream &operator<<(std::ostream &os, const vec &v) {
    if (N > 0) {
      os << v.array[0];
      for (int i = 1; i < N; ++i) { os << " " << v.array[i]; }
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, vec &v) {
    if (N > 0) {
      is >> v.array[0];
      for (int i = 1; i < N; ++i) { is.ignore(); is >> v.array[i]; }
    }
    return is;
  }
};

typedef vec<3,int> v3i;
typedef vec<3,double> v3d;

typedef vec<4,int> v4i;
typedef vec<4,double> v4d;

#endif  // VEC_HH_
