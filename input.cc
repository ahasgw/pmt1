#include "input.hh"

#include <algorithm>
#include <istream>
#include <sstream>

void InputXYZ(std::istream &is,
              Ptcls *ptcls,
              int *total_ptcl,
              v3r &div_min,
              v3r &div_max,
              int comm_rank,
              int comm_size,
              MPI_Comm comm) {
  using namespace std;

  // all particles
  Ptcls all_ptcls;

  if (comm_rank == 0) {

    // input XYZ format
    if (is) {
      enum { E1stLine, E2ndLine, ELaterLine } st = E1stLine;
      string line;
      while (getline(is, line)) {
        istringstream iss(line);
        switch (st) {
          case E1stLine: {
            iss >> *total_ptcl;
            all_ptcls.reserve(*total_ptcl);
            st = E2ndLine;
            break;
          }
          case E2ndLine: {
            string comment = line;
            st = ELaterLine;
            break;
          }
          case ELaterLine:
          default: {
            int type;
            Ptcl ptcl;
            iss >> type >> ptcl;
            all_ptcls.push_back(ptcl);
            break;
          }
        }
      }
    }
    MPI_Bcast(total_ptcl, 1, MPI_INT, 0, comm);
    MPI_Bcast(&all_ptcls[0], *total_ptcl * sizeof(Ptcl), MPI_BYTE, 0, comm);
  }
  else {  // other ranks
    MPI_Bcast(total_ptcl, 1, MPI_INT, 0, comm);
    all_ptcls.resize(*total_ptcl);
    MPI_Bcast(&all_ptcls[0], *total_ptcl * sizeof(Ptcl), MPI_BYTE, 0, comm);
  }

  ptcls->clear();
  for (int n = 0; n < *total_ptcl; ++n) {
    Ptcl &p = all_ptcls[n];
    if (div_min <= p.crd && p.crd < div_max) {
      ptcls->push_back(p);
    }
  }
}
