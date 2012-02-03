#include "cart.hh"
#include "conf.hh"

#if 0
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
#endif


CartNode::CartNode(Conf &conf): conf_(conf) {
  // setup timer
  t_init.Label("cart init").Comm(conf_.cart_comm).Start();
  t_step.Label("cart step").Comm(conf_.cart_comm);

  MPI_Comm_rank(conf_.cart_comm, &cart_rank);
  MPI_Comm_size(conf_.cart_comm, &cart_size);
  MPI_Cart_coords(conf_.cart_comm, cart_rank, 3, cart_pos);

  // topology

  // setup node
  div_min = div_max = conf_.sys_min;
  div_min += v3d(cart_pos    ) * conf_.sys_size / v3d(conf_.node_num);
  div_max += v3d(cart_pos + 1) * conf_.sys_size / v3d(conf_.node_num);

#if 1
  std::cout << "<" << cart_rank << ">\t" << div_min << "\t"
      << div_max << std::endl;
#endif

  GenerateParticles();

  t_init.Stop();
}

CartNode::~CartNode() {
  // print timer
  if (conf_.verbose > 1) t_step.PrintAll("# ", conf_.max_step);
  if (conf_.verbose > 0) t_step.PrintMax("# max ", conf_.max_step);
  if (conf_.verbose > 1) t_init.PrintAll("# ");
  if (conf_.verbose > 0) t_init.PrintMax("# max ");
}

void CartNode::StepForward(int t) {
  using namespace std;
  t_step.Start();
  if (cart_rank == 0) cout << "step " << t << "\n" << flush;
  t_step.Stop();
}

void CartNode::StepBackward(int t) {
  using namespace std;
  t_step.Start();
  if (cart_rank == 0) cout << "step " << t << "\n" << flush;
  t_step.Stop();
}

void CartNode::GenerateParticles() {
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
