#include "output.hh"

#include <algorithm>
#include <limits>
#include <ostream>

void OutputXYZ(std::ostream &os,
               std::ostream &rs,
               const char *comment,
               const Ptcls &ptcls,
               int total_ptcl,
               int step,
               int max_step,
               int comm_rank,
               int comm_size,
               MPI_Comm comm) {
  using namespace std;

  // calculate receive count
  int rc = static_cast<int>(ptcls.size() * sizeof(Ptcl));

  if (comm_rank == 0) {
    // gather receive counts
    vector<int> rcounts(comm_size);
    MPI_Gather(&rc, 1, MPI_INT, &rcounts[0], 1, MPI_INT, 0, comm);

    // make displacement list
    vector<int> displs(comm_size);
    displs[0] = 0;
    for (vector<int>::size_type i = 1; i < rcounts.size(); ++i)
      displs[i] = displs[i - 1] + rcounts[i - 1];

    // gather all particles
    Ptcls all_ptcls(total_ptcl);
    MPI_Gatherv(const_cast<Ptcl*>(&ptcls[0]), rc, MPI_BYTE,
                &all_ptcls[0], &rcounts[0], &displs[0], MPI_BYTE, 0, comm);

    // sort particles
    sort(all_ptcls.begin(), all_ptcls.end(), LessId());

    // output XYZ format
    if (os) {
      os << total_ptcl << "\n" << comment << " (step " << step << ")\n";
      for (Ptcls::size_type i = 0; i < all_ptcls.size(); ++i) {
        os << 1 << "\t" << all_ptcls[i].crd << "\n";
      }
      os << flush;
    }

    // output XYZ format for restart save
    if ((step == max_step) && rs) {
      rs << total_ptcl << "\n" << comment << " (step " << step << ")\n";
      streamsize ss = rs.precision(numeric_limits<real_t>::digits10);
      for (Ptcls::size_type i = 0; i < all_ptcls.size(); ++i) {
        rs << 1 << "\t" << all_ptcls[i] << "\n";
      }
      rs.precision(ss);
      rs << flush;
    }
  }
  else {
    // gather receive counts
    MPI_Gather(&rc, 1, MPI_INT, 0, 1, MPI_INT, 0, comm);

    // gather all particles
    MPI_Gatherv(const_cast<Ptcl*>(&ptcls[0]), rc, MPI_BYTE,
                0, 0, 0, MPI_BYTE, 0, comm);
  }
}
