#include "cart.hh"
#include "conf.hh"

CartNode::CartNode(Conf &conf): conf_(conf) {
  // setup timer
  t_init.Label("cart init").Comm(conf_.cart_comm).Start();
  t_step.Label("cart step").Comm(conf_.cart_comm);

  MPI_Comm_rank(conf_.cart_comm, &cart_rank);
  MPI_Comm_size(conf_.cart_comm, &cart_size);
  MPI_Cart_coords(conf_.cart_comm, cart_rank, 3, cart_pos);

#if 1
  if (cart_rank == 0) {
    if (conf_.verbose > 0) std::cout << "# cart_size\t" << cart_size << "\n";
  }
#endif

  // topology

  // setup node
  div_min = div_max = conf_.sys_min;
  div_min += v3d(cart_pos    ) * conf_.sys_size / v3d(conf_.cart_num);
  div_max += v3d(cart_pos + 1) * conf_.sys_size / v3d(conf_.cart_num);

#if 0
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
}
