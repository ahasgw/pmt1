#ifndef CART_HH_
#define CART_HH_ 1

#include "node.hh"
#include "ptcl.hh"
#include "timer.hh"
#include "vec.hh"

class CartNode: public WorkNode {
 private:
  v3d div_min;
  v3d div_max;
  v3i cart_pos;
  Ptcls ptcls;
  Conf &conf_;
  int cart_rank;
  int cart_size;
  Timer t_init;
  Timer t_step;

 public:
  CartNode(Conf &conf);
  ~CartNode();

  void StepForward(int t);
  void StepBackward(int t);

 private:
  void GenerateParticles();
};

#endif  // CART_HH_
