#include "node.hh"
#include <cmath>
#include <cstdlib>

namespace {

inline double Rand() {
  return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
}

double Gaussian(double sgm, double myu) {
  static int binary = 0;
  static double g = myu;
  if ((binary++ & 1) == 0) {
    double a1, a2, b;
    do {
      a1 = 2.0 * Rand() - 1.0;
      a2 = 2.0 * Rand() - 1.0;
      b = a1 * a1 + a2 * a2;
    } while (b >= 1.0);
    b = std::sqrt((-2.0 * std::log(b)) / b);
    g = a1 * b * sgm + myu;
    return a2 * b * sgm + myu;
  }
  return g;
}

}  // namespace

void Node::GenerateParticles(int total_ptcl) {
#if 0
  for (int i = 0; i < total_ptcl; ++i) {
    PtclVelIdAttr p;
    for (int j = 0; j < 3; ++j) {
      p.vel[j] = Gaussian() * mvf * min_sys_size;
      p.crd[j] = Rand() * sys_size + sys_ofst;
    }
  }
#endif
}

void Node::StepForward(int t) {
}
