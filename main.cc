#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "conf.hh"
#include "node.hh"
#include "signal.hh"


// signal handling
volatile std::sig_atomic_t signal_raised = 0;
extern "C" void SetOneIfSignalRaised(int /* signum */) { signal_raised = 1; }


int main(int argc, char *argv[]) {
  using namespace std;
  MPI_Init(&argc, &argv);

  {
    // parse arguments and setup configuration
    Conf conf(argc, argv, MPI_COMM_WORLD);

    // construct node data
    Node node(conf);

    // iterate until maximum time step
    StSignalHandler(SIGINT, SetOneIfSignalRaised);
    for (int t = 1; (!signal_raised) && (t <= conf.max_step); ++t) {
      node.StepForward(t);
    }
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
