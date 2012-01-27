#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <unistd.h>
#include <vector>
#include "conf.hh"
#include "node.hh"
#include "timer.hh"


volatile std::sig_atomic_t signal_raised = 0;
extern "C" void SignalHandler(int /* signum */) { signal_raised = 1; }


int main(int argc, char *argv[]) {
  using namespace std;
  MPI_Init(&argc, &argv);

  // setup timer
  Timer t_total("total"); t_total.Start();

  // parse arguments
  Conf conf(argc, argv);
  if (conf.verbose > 0) conf.Print();

  // GenerateParticles
  Node node;
  node.GenerateParticles(conf.total_ptcl);

  // Iterate for time steps
  signal(SIGINT, SignalHandler);  // set signal handler
  for (int t = 1; (!signal_raised) && (t <= conf.max_step); ++t) {
    node.StepForward(t);
  }
  signal(SIGINT, SIG_DFL);  // reset signal handler

  // print timer
  t_total.Stop();
  if (conf.verbose > 1) t_total.Print("# ", conf.max_step);
  if (conf.verbose > 0) t_total.PrintMax("# max ", conf.max_step);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
