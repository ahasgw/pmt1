#ifndef TIMER_HH_
#define TIMER_HH_ 1

#include <cerrno>
#include <ctime>
#include <mpi.h>
#include <iostream>
#include <ostream>

namespace {

inline double gettime() {
  timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return static_cast<double>(ts.tv_sec) +
      static_cast<double>(ts.tv_nsec) * 1.0e-9;
}

}  // namespace

class Timer {
 public:
  Timer(const char *label = "", MPI_Comm comm = MPI_COMM_WORLD)
      : label_(label), comm_(comm) {
    Clear();
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
  }
  ~Timer() {}

  Timer &Clear() { total_ = 0.0; return *this; }
  Timer &Start() { start_ = gettime(); return *this; }
  Timer &Stop() { total_ += (gettime() - start_); return *this; }

  const char *Label() const { return label_; }
  Timer &Label(const char *label) { label_ = label; return *this; }

  MPI_Comm Comm() const { return comm_; }
  Timer &Comm(MPI_Comm comm) { comm_ = comm; return *this; }

  double Sec() const { return total_; }

  Timer &PrintAll(const char *hdr = "", int num_step = 1,
                  std::ostream &os = std::cout) {
    using namespace std;
    os << flush;
    for (int i = 0; i < size_; ++i) {
      MPI_Barrier(comm_);
      if (i == rank_) {
        os << hdr << "[" << rank_ << "] " << label_ << "\t" << total_ << " sec";
        if (num_step > 1) os << "\t" << (total_ / num_step) << " sec/step";
        os << endl;
      }
    }
    return *this;
  }

  Timer &PrintMax(const char *hdr = "", int num_step = 1,
                  std::ostream &os = std::cout) {
    using namespace std;
    double max = 0.0;
    MPI_Reduce(&total_, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
    os << flush;
    if (rank_ == 0) {
      os << hdr << label_ << "\t" << max << " sec";
      if (num_step > 1) os << "\t" << (max / num_step) << " sec/step";
      os << endl;
    }
    return *this;
  }

 private:
  const char *label_;
  MPI_Comm comm_;
  double total_;
  double start_;
  int rank_;
  int size_;
};

#endif  // TIMER_HH_
