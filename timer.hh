#ifndef TIMER_HH_
#define TIMER_HH_ 1

#include <cerrno>
#include <ctime>
#include <mpi.h>

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

  void Clear() { total_ = 0.0; }
  void Start() { start_ = gettime(); }
  void Stop() { total_ += (gettime() - start_); }
  const char *Label() const { return label_; }
  void Label(const char *label) { label_ = label; }
  double Sec() const { return total_; }

  void Print(const char *hdr = "", int num_step = 1);
  void PrintMax(const char *hdr = "", int num_step = 1);

 private:
  const char *label_;
  MPI_Comm comm_;
  double total_;
  double start_;
  int rank_;
  int size_;
};

#endif  // TIMER_HH_
