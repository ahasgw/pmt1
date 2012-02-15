#ifndef PTCL_HH_
#define PTCL_HH_ 1

#include "vec.hh"
#include <vector>

struct Ptcl {
  v3d crd;
  v3d vel;
  int id;
  unsigned attr;
};

struct LessId {
  bool operator()(const Ptcl &p0, const Ptcl &p1) const {
    return p0.id < p1.id;
  }
};

struct LessAttr {
  bool operator()(const Ptcl &p0, const Ptcl &p1) const {
    return p0.attr < p1.attr;
  }
};

typedef std::vector<Ptcl> Ptcls;

#endif  // PTCL_HH_
