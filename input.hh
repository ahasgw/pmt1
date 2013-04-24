#ifndef INPUT_HH_
#define INPUT_HH_ 1

#include <mpi.h>
#include <iosfwd>
#include "ptcl.hh"

void InputXYZ(std::istream &is,
              Ptcls *ptcls,
              int *total_ptcl,
              v3r &div_min,
              v3r &div_max,
              v3r &sys_min,
              v3r &sys_max,
              v3i &boundary,
              int comm_rank,
              int comm_size,
              MPI_Comm comm);

#endif  // INPUT_HH_
