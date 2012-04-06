#ifndef SIGNAL_HH_
#define SIGNAL_HH_ 1

#include <csignal>

class StSignalHandler {
 private:
  void (*prev_handler)(int);
  int signum_;
 public:
  StSignalHandler(int signum, void (*handler)(int)): signum_(signum) {
    prev_handler = std::signal(signum_, handler);
  }
  ~StSignalHandler() {
    std::signal(signum_, prev_handler);
  }
};

#endif  // SIGNAL_HH_
