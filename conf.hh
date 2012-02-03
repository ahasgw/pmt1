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
#include "timer.hh"

class Conf {
 public:
  v3d sys_ofst;
  v3d sys_size;
  v3d sys_min;
  v3d sys_max;
  v3i cart_num;
  v3i periods;
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
  MPI_Comm cart_comm;
  Timer t_total;
  Timer t_conf;

  enum { IDLE_NODE = 0, CART_NODE } node_type;

 public:

  Conf(int argc, char *argv[], MPI_Comm comm = MPI_COMM_WORLD)
      : sys_ofst(-50.0),
        sys_size(100.0),
        sys_min(-50.0),
        sys_max(-50.0 + 100.0),
        cart_num(0),
        periods(true),
        oname(NULL),
        max_step(1),
        total_ptcl(10000),
        global_seed(1),
        format(0),
        verbose(0),
        argc(argc),
        argv(argv),
        comm(comm),
        cart_comm(MPI_COMM_NULL),
        node_type(IDLE_NODE)
  {
    // setup timer
    t_total.Label("total").Start();
    t_conf.Label("config").Start();

    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    ParseArguments(argc, argv);

    sys_min = sys_ofst;
    sys_max = sys_ofst + sys_size;

    std::srand(static_cast<unsigned>(global_seed));  // set random seed

    int num_node = DeterminNumberOfNode();

    MPI_Dims_create(num_node, 3, cart_num);
    MPI_Cart_create(comm, 3, cart_num, periods, true, &cart_comm);
    node_type = (cart_comm != MPI_COMM_NULL) ? CART_NODE : IDLE_NODE;

    if (verbose > 0) Print();
    t_conf.Stop();
  }

  ~Conf() {
    if (cart_comm != MPI_COMM_NULL) MPI_Comm_free(&cart_comm);

    t_total.Stop();

    // print timer
    if (verbose > 1) t_conf.PrintAll("# ");
    if (verbose > 0) t_conf.PrintMax("# max ");
    if (verbose > 1) t_total.PrintAll("# ");
    if (verbose > 0) t_total.PrintMax("# max ");
  }

  void Print(std::ostream &os = std::cout) const {
    os << *this << std::flush;
  }

  friend std::ostream &operator<<(std::ostream &os, const Conf &c) {
    if (c.comm_rank == 0) {
      os << "# command_line\t" << c.argv[0];
      for (char **p = c.argv + 1; (*p); ++p) { os << " " << *p; } os << "\n";
      os << "# sys_size\t" << c.sys_size << "\n";
      os << "# sys_ofst\t" << c.sys_ofst << "\n";
      os << "# sys_min\t" << c.sys_min << "\n";
      os << "# sys_max\t" << c.sys_max << "\n";
      os << "# cart_num\t" << c.cart_num << "\n";
      os << "# periods\t" << c.periods << "\n";
      os << "# max_step\t" << c.max_step << "\n";
      os << "# total_ptcl\t" << c.total_ptcl << "\n";
      os << "# global_seed\t" << c.global_seed << "\n";
      os << "# oname\t" << (c.oname ? c.oname : "") << "\n";
      os << "# format\t" << c.format << "\n";
      os << "# verbose\t" << c.verbose << "\n";
      os << "# comm_size\t" << c.comm_size << "\n";
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
        case 'N': { istringstream iss(::optarg); iss >> cart_num; } break;
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
                "  -S <X:Y:Z>    system size\n"
                "  -O <X:Y:Z>    system offset\n"
                "  -N <X:Y:Z>    number of nodes in Cartesian grid\n"
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

  int DeterminNumberOfNode() {
    using namespace std;
    int num_node = 1;
    for (int i = 0; i < 3; ++i) {
      cart_num[i] = abs(cart_num[i]);
      if (cart_num[i] > 1)
        num_node *= cart_num[i];
    }
    if (num_node > comm_size) {
      if (comm_rank == 0)
        cout << "pmt0: number of nodes exceeds communicator size. abort\n"
            << flush;
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (cart_num[0] * cart_num[1] * cart_num[2] == 0) {
      if (cart_num[0] + cart_num[1] + cart_num[2] != 0)
        num_node *= (comm_size / num_node);
      else
        num_node = comm_size;
    }
    return num_node;
  }
};

#endif  // CONF_HH_
