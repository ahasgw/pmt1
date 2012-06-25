#ifndef OUTPUT_HH_
#define OUTPUT_HH_ 1

#include <mpi.h>
#include <iosfwd>
#include "ptcl.hh"

void OutputXYZ(std::ostream &os,
               std::ostream &rs,
               const char *comment,
               const Ptcls &ptcls,
               int total_ptcl,
               int step,
               int max_step,
               int comm_rank,
               int comm_size,
               MPI_Comm comm);

#endif  // OUTPUT_HH_
