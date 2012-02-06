#ifndef RANDOM_HH_
#define RANDOM_HH_ 1

#include <cmath>
#include <cstdlib>

inline double Rand() {
  return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
}

inline double Gaussian(double mu = 0.0, double sigma = 1.0) {
  double x, y, r2;
  do {
    x = 2.0 * Rand() - 1.0;
    y = 2.0 * Rand() - 1.0;
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0.0);
  return sigma * y * sqrt(-2.0 * log(r2) / r2);
}

#endif  // RANDOM_HH_
