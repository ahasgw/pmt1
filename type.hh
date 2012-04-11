#ifndef TYPE_HH_
#define TYPE_HH_ 1

#include "vec.hh"

//#define REAL_IS_FLOAT 1

#ifdef REAL_IS_FLOAT
typedef float real_t;
#else
typedef double real_t;
#endif

typedef vec<2,int> v2i;
typedef vec<3,int> v3i;
typedef vec<4,int> v4i;

typedef vec<2,float> v2f;
typedef vec<3,float> v3f;
typedef vec<4,float> v4f;

typedef vec<2,double> v2d;
typedef vec<3,double> v3d;
typedef vec<4,double> v4d;

typedef vec<2,real_t> v2r;
typedef vec<3,real_t> v3r;
typedef vec<4,real_t> v4r;

#endif  // TYPE_HH_
