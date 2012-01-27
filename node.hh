#ifndef NODE_HH_
#define NODE_HH_ 1

#include "vec.hh"
#include <vector>

struct Ptcl {
  double chg;
  v3d crd;
};

typedef std::vector<Ptcl> Ptcls;

struct PtclVel {
  double chg;
  v3d crd;
  v3d vel;
};

typedef std::vector<PtclVel> PtclVels;

struct PtclVelIdAttr {
  double chg;
  v3d crd;
  v3d vel;
  int id;
  unsigned attr;
};

typedef std::vector<PtclVelIdAttr> PtclVelIdAttrs;





class Node {
 public:
  Node() {}
  ~Node() {}

  void GenerateParticles(int total_ptcl);
  void StepForward(int t);
};

#endif  // NODE_HH_
