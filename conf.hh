#ifndef CONF_HH_
#define CONF_HH_ 1

#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include <unistd.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include "vec.hh"

class Conf {
 public:
  v3d sys_size;
  v3d sys_ofst;
  v3i node_num;
  v3i node_pos;
  const char *oname;
  int max_step;
  int total_ptcl;
  int global_seed;
  int format;
  int verbose;
  int comm_rank;
  int comm_size;
  int argc;
  char **argv;
  MPI_Comm comm;

 public:

  Conf(int argc, char *argv[], MPI_Comm comm = MPI_COMM_WORLD)
      : sys_size(100.0),
        sys_ofst(-50.0),
        node_num(1),
        node_pos(0),
        oname(NULL),
        max_step(1),
        total_ptcl(10000),
        global_seed(1),
        format(0),
        verbose(0),
        argc(argc),
        argv(argv),
        comm(comm) {
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    ParseArguments(argc, argv);
    Setup();
  }

  ~Conf() {}

  void Print(std::ostream &os = std::cout) const {
    os << *this << std::flush;
  }

  friend std::ostream &operator<<(std::ostream &os, const Conf &c) {
    if (c.comm_rank == 0) {
      os << "# command_line\t" << c.argv[0];
      for (char **p = c.argv + 1; (*p); ++p) { os << " " << *p; } os << "\n";
      os << "# sys_size\t" << c.sys_size << "\n";
      os << "# sys_ofst\t" << c.sys_ofst << "\n";
      os << "# node_num\t" << c.node_num << "\n";
      os << "# max_step\t" << c.max_step << "\n";
      os << "# total_ptcl\t" << c.total_ptcl << "\n";
      os << "# global_seed\t" << c.global_seed << "\n";
      os << "# oname\t" << (c.oname ? c.oname : "") << "\n";
      os << "# format\t" << c.format << "\n";
      os << "# verbose\t" << c.verbose << "\n";
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, Conf &c) {
    return is;
  }

 private:

  void ParseArguments(int argc, char *argv[]) {
    using namespace std;
    for (::opterr = 0;;) {
      int opt = ::getopt(argc, argv, ":m:n:S:O:N:s:o:f:dvh");
      if (opt == -1) break;
      switch (opt) {
        case 'm': max_step = abs(atoi(::optarg)); break;
        case 'n': total_ptcl = abs(atoi(::optarg)); break;
        case 'S': { istringstream iss(::optarg); iss >> sys_size; } break;
        case 'O': { istringstream iss(::optarg); iss >> sys_ofst; } break;
        case 'N': { istringstream iss(::optarg); iss >> node_num; } break;
        case 's': global_seed = atoi(::optarg); break;
        case 'o': oname = ::optarg; break;
        case 'f': format = atoi(::optarg); break;
        case 'v': ++verbose; break;
        case 'h': {
          if (comm_rank == 0) {
            cout <<
                "This is pmt0. A particle-moving test program.\n"
                "Usage: pmt0 [options]\n"
                "Options:\n"
                "  -m <n>        maximum number of step\n"
                "  -n <n>        total number of particles\n"
                "  -S <X:Y:Z>    system size in 3d\n"
                "  -O <X:Y:Z>    system offset in 3d\n"
                "  -N <X:Y:Z>    number of nodes in 3d\n"
                "  -s <n>        random seed\n"
                "  -o <name>     output file/directory name\n"
                "  -f <0/1>      output format; 0:XYZ(default), 1:CDV\n"
                "  -v            print message verbosely\n"
                "  -h            show this help message\n"
                << flush;
          }
          MPI_Finalize();
          exit(EXIT_SUCCESS);
        }
        case ':': {  // missing option argument
          if (comm_rank == 0)
            cout << "pmt0: option requires an argument -- '"
                << ::optopt << "'.  try '-h' for help\n" << flush;
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        default: /* case '?': */ {  // unknown option
          if (comm_rank == 0)
            cout << "pmt0: invalid option -- '"
                << ::optopt << "'.  try '-h' for help\n" << flush;
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  void Setup() {
    // set random seed
    std::srand(static_cast<unsigned>(global_seed));
    int num_node = DeterminNumberOfNode();
  }

  int DeterminNumberOfNode() const {
    using namespace std;
    int num_node = 1;
    for (int i = 0; i < 3; ++i) {
      if (node_num[i] > 1)
        num_node *= node_num[i];
    }
    if (num_node > comm_size) {
      if (comm_rank == 0)
        cout << "pmt0: number of nodes exceeds communicator size. abort\n"
            << flush;
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (node_num[0] * node_num[1] * node_num[2] == 0)
      num_node *= (comm_size / num_node);
    return num_node;
  }
};

#endif  // CONF_HH_
