#ifndef PTCL_HH_
#define PTCL_HH_ 1

#include "type.hh"
#include <vector>

struct Ptcl {
  v3r crd;
  real_t chg;
  v3r vel;
  real_t inv_2mass;
  int id;
  unsigned attr;

  enum AttrMask { ATTR = 0xFFFFFFFF, DIR = 0x3F };
};

struct LessId {
  bool operator()(const Ptcl &p0, const Ptcl &p1) const {
    return p0.id < p1.id;
  }
};

template<unsigned MASK>
struct Less {
  bool operator()(const Ptcl &p0, const Ptcl &p1) const {
    return (p0.attr & MASK) < (p1.attr & MASK);
  }
};

typedef std::vector<Ptcl> Ptcls;

#endif  // PTCL_HH_
