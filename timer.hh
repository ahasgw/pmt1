#ifndef TIMER_HH_
#define TIMER_HH_ 1

#include <cerrno>
#include <ctime>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

#ifdef __MACH__
# include <sys/time.h>
#endif  // __MACH__

namespace {

inline double gettime() {
  timespec ts;
#ifdef __MACH__
  {
    timeval tv;
    gettimeofday(&tv, 0);
    ts.tv_sec = tv.tv_sec;
    ts.tv_nsec = tv.tv_usec * 1000;
  }
#else  // __MACH__
  clock_gettime(CLOCK_REALTIME, &ts);
#endif  // __MACH__
  return static_cast<double>(ts.tv_sec) +
      static_cast<double>(ts.tv_nsec) * 1.0e-9;
}

}  // namespace

class Timer {
 public:
  Timer(const char *label = "", MPI_Comm comm = MPI_COMM_WORLD)
      : label_(label), comm_(comm)
  {
    Clear();
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
  }
  ~Timer() {}

  Timer &Clear() { total_ = 0.0; count_ = 0; return *this; }
  Timer &Start() { start_ = gettime(); ++count_; return *this; }
  Timer &Stop() { total_ += (gettime() - start_); return *this; }

  const char *Label() const { return label_; }
  Timer &Label(const char *label) { label_ = label; return *this; }

  MPI_Comm Comm() const { return comm_; }
  Timer &Comm(MPI_Comm comm) { comm_ = comm; return *this; }

  double Sec() const { return total_; }
  int count() const { return count_; }

  Timer &PrintAll(const char *hdr = "", std::ostream &os = std::cout) {
    using namespace std;
    os << flush;
    for (int i = 0; i < size_; ++i) {
      MPI_Barrier(comm_);
      if (i == rank_) {
        ostringstream oss;
        oss << hdr << "[" << rank_ << "] " << label_;
        os << setw(40) << left << oss.str();
        os << "  " << setw(11) << total_ << " sec";
        if (count_ > 1)
          os << "  " << setw(11) << (total_ / count_) << " sec/step";
        os << endl;
      }
    }
    return *this;
  }

  Timer &PrintMax(const char *hdr = "", std::ostream &os = std::cout) {
    using namespace std;
    double max = 0.0;
    MPI_Reduce(&total_, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
    os << flush;
    if (rank_ == 0) {
      ostringstream oss;
      oss << hdr << label_;
      os << setw(40) << left << oss.str();
      os << "  " << setw(11) << max << " sec";
      if (count_ > 1)
        os << "  " << setw(11) << (max / count_) << " sec/step";
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
  int count_;
};

#endif  // TIMER_HH_
