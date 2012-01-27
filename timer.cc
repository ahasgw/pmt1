#include "timer.hh"
#include <iostream>
#include <ostream>

void Timer::Print(const char *hdr, int num_step) {
  using namespace std;
  cout << flush;
  for (int i = 0; i < size_; ++i) {
    MPI_Barrier(comm_);
    if (i == rank_) {
      cout << hdr << "[" << rank_ << "] " << label_ << "\t"
          << total_ << " sec\n";
      if (num_step > 1) {
        cout << hdr << "[" << rank_ << "] " << label_ << "\t"
            << (total_ / num_step) << " sec/step\n";
      }
      cout << flush;
    }
  }
}

void Timer::PrintMax(const char *hdr, int num_step) {
  using namespace std;
  double max = 0.0;
  MPI_Reduce(&total_, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
  cout << flush;
  if (rank_ == 0) {
    cout << hdr << label_ << "\t" << max << " sec\n";
    if (num_step > 1) {
      cout << hdr << label_ << "\t" << (max / num_step) << " sec/step\n";
    }
    cout << flush;
  }
}
