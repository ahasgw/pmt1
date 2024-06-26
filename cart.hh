#ifndef CART_HH_
#define CART_HH_ 1

#include "node.hh"

#include <fstream>
#include "ptcl.hh"
#include "timer.hh"

typedef std::vector<v3r> Force;

class CartNode: public WorkNode {
 private:
  v3r div_size;
  v3r div_size_2;
  v3r div_min;
  v3r div_max;
  v3r sys_size;
  v3r sys_size_2;
  v3r sys_min;
  v3r sys_max;
  real_t dt;
  real_t cutoff2;
  v3i cart_pos;
  v3i boundary;
  v3i div_intr;
  Ptcls ptcls;
  Ptcls recv_buff[26];  // buffers are used in order of arrival
  Force force;
  Timer t_cfrc;
  Timer t_comm;
  Timer t_step;
  Timer t_oput;
  Timer t_init;
  Conf &conf_;
  std::ifstream is;
  std::ofstream os;
  std::ofstream rs;
  MPI_Comm cart_comm;
  int cart_rank;
  int cart_size;
  int steps_to_write;
  int intr_size;
  int use_ringcomm;

  enum { MIDDLE_DIR = 0x00, LOWER_DIR = 0x02, UPPER_DIR = 0x03 };
  struct Connect {
    unsigned dir_tag;
    int send_to;
    int recv_from;
    MPI_Request req;
  };
  typedef std::vector<Connect> Conns;
  Conns conns;  // elements are sorted by dir_tag

  struct Interact {
    int send_to;
    int recv_from;
    MPI_Request req;
  };
  typedef std::vector<Interact> Intrs;
  Intrs intrs;

 public:
  CartNode(Conf &conf);
  ~CartNode();

  void StepForward(int t);

 private:
  void InitConnect();
  int GetCartRankAtOffset(int dx, int dy, int dz);
  void InitInteract();
  int GetCartRankAtOffsetInsideBoundary(int dx, int dy, int dz);
  void GenerateParticles();
  void ExchangeParticles();
  void CalculateForce();

  void CalculateForceCutoffPeriodic();
  void CalculateForceCutoffPeriodic_RingComm();
};

#endif  // CART_HH_
