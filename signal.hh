#ifndef SIGNAL_HH_
#define SIGNAL_HH_ 1

#include <csignal>

class StSignalHandler {
 private:
  sighandler_t prev_handler;
  int signum_;
 public:
  StSignalHandler(int signum, sighandler_t handler): signum_(signum) {
    prev_handler = std::signal(signum_, handler);
  }
  ~StSignalHandler() {
    std::signal(signum_, prev_handler);
  }
};

#endif  // SIGNAL_HH_
