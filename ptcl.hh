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

  friend std::ostream &operator<<(std::ostream &os, const Ptcl &p) {
    os << p.crd << "\t" << p.vel << "\t" <<
        p.chg << "\t" << p.inv_2mass << "\t" << p.id;
    return os;
  }
  friend std::istream &operator>>(std::istream &is, Ptcl &p) {
    p.attr = 0U;
    is >> p.crd >> p.vel >> p.chg >> p.inv_2mass >> p.id;
    return is;
  }
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
