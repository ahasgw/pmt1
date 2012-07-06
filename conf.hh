#ifndef CONF_HH_
#define CONF_HH_ 1

#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include <unistd.h>
#include <iostream>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include "type.hh"
#include "timer.hh"

class Conf {
 public:
  v3r sys_ofst;
  v3r sys_size;
  v3r sys_min;
  v3r sys_max;
  real_t cutoff;
  v3i cart_num;
  v3i periodic;
  int rest_num;
  int max_step;
  int total_ptcl;
  int neutralize;
  int global_seed;
  int write_interval;
  int write_step0;
  int verbose;
  int comm_rank;
  int comm_size;
  int argc;
  char **argv;
  std::string ifname;
  std::string ofname;
  std::string rfname;
  std::string cmd_line;
  MPI_Comm comm;
  MPI_Comm work_comm;
  MPI_Comm cart_comm;
  Timer t_total;
  Timer t_conf;

  enum { IDLE_NODE = 0, CART_NODE, REST_NODE } node_type;

 public:

  Conf(int argc, char *argv[], MPI_Comm comm = MPI_COMM_WORLD)
      : sys_ofst(-50.0),
        sys_size(100.0),
        sys_min(-50.0),
        sys_max(-50.0 + 100.0),
        cutoff(0.0),
        cart_num(0),
        periodic(true),
        rest_num(0),
        max_step(1),
        total_ptcl(10000),
        neutralize(0),
        global_seed(1),
        write_interval(1),
        write_step0(0),
        verbose(0),
        comm_rank(0),
        comm_size(1),
        argc(argc),
        argv(argv),
        ifname(""),
        ofname(""),
        rfname(""),
        cmd_line(""),
        comm(comm),
        work_comm(MPI_COMM_NULL),
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

    if (cutoff > sys_size.min() / 2.0) {
      if (comm_rank == 0) {
        std::cout << "pmt: cutoff radius " << cutoff <<
            " exceeds half of minimum system size " << sys_size.min() / 2.0 <<
            ". abort\n" << std::flush;
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    std::srand(static_cast<unsigned>(global_seed));  // set random seed

    int num_cart_node = DeterminNumberOfCartNode();
    rest_num = comm_size - num_cart_node;
    node_type = ((comm_rank < num_cart_node) ? CART_NODE : REST_NODE);
    MPI_Comm_split(MPI_COMM_WORLD, node_type, 0, &work_comm);

    MPI_Dims_create(num_cart_node, 3, cart_num);
    if (node_type == CART_NODE) {
      MPI_Cart_create(work_comm, 3, cart_num, periodic, true, &cart_comm);
    }

    if (verbose > 0) Print();
    t_conf.Stop();
  }

  ~Conf() {
    if (cart_comm != MPI_COMM_NULL) MPI_Comm_free(&cart_comm);
    if (work_comm != MPI_COMM_NULL) MPI_Comm_free(&work_comm);

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
      os << "# cmd_line\t" << c.cmd_line << "\n";
      os << "# sys_size\t" << c.sys_size << "\n";
      os << "# sys_ofst\t" << c.sys_ofst << "\n";
      os << "# sys_min\t" << c.sys_min << "\n";
      os << "# sys_max\t" << c.sys_max << "\n";
      os << "# cutoff\t" << c.cutoff << "\n";
      os << "# cart_num\t" << c.cart_num << "\n";
      os << "# rest_num\t" << c.rest_num << "\n";
      os << "# periodic\t" << c.periodic << "\n";
      os << "# max_step\t" << c.max_step << "\n";
      os << "# neutralize\t" << c.neutralize << "\n";
      os << "# global_seed\t" << c.global_seed << "\n";
      os << "# ifname\t" << c.ifname << "\n";
      os << "# ofname\t" << c.ofname << "\n";
      os << "# rfname\t" << c.rfname << "\n";
      os << "# write_interval\t" << c.write_interval << "\n";
      os << "# write_step0\t" << c.write_step0 << "\n";
      os << "# verbose\t" << c.verbose << "\n";
      os << "# comm_size\t" << c.comm_size << "\n";
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, Conf &c) {
    return is;
  }

 private:

  template<typename T>
  void Read(T &t) {
    std::istringstream iss(::optarg);
    iss.exceptions(std::ios::badbit | std::ios::failbit);
    iss >> t;
  }

  template<typename T>
  void ReadAbs(T &t) {
    std::istringstream iss(::optarg);
    iss.exceptions(std::ios::badbit | std::ios::failbit);
    iss >> t;
    t = std::abs(t);
  }

  void ParseArguments(int argc, char *argv[]) {
    using namespace std;
    // record command line
    cmd_line += argv[0];
    for (char **p = argv + 1; (*p); ++p) { cmd_line += " "; cmd_line += *p; }
    // set default floating-point number precision
    cout.precision(numeric_limits<double>::digits10);

    for (::opterr = 0;;) {
      int opt = ::getopt(argc, argv, ":m:p:S:O:N:P:c:ns:i:o:r:w:0dvh");
      if (opt == -1) break;
      try {
        switch (opt) {
          case 'm': ReadAbs(max_step); break;
          case 'p': ReadAbs(total_ptcl); break;
          case 'S': Read(sys_size); break;
          case 'O': Read(sys_ofst); break;
          case 'N': Read(cart_num); break;
          case 'P': Read(periodic); break;
          case 'c': Read(cutoff); break;
          case 'n': ++neutralize; break;
          case 's': ReadAbs(global_seed); break;
          case 'i': ifname = ::optarg; break;
          case 'o': ofname = ::optarg; break;
          case 'r': rfname = ::optarg; break;
          case 'w': ReadAbs(write_interval); break;
          case '0': ++write_step0; break;
          case 'v': ++verbose; break;
          case 'h': {
            if (comm_rank == 0) {
              cout <<
                  "This is pmt. A particle-moving test program.\n"
                  "Usage: pmt [options]\n"
                  "Options:\n"
                  "  -m <n>        maximum number of step                 (1)\n"
                  "  -p <n>        total number of particles          (10000)\n"
                  "  -S <X:Y:Z>    system size               (100.:100.:100.)\n"
                  "  -O <X:Y:Z>    system offset             (-50.:-50.:-50.)\n"
                  "  -N <X:Y:Z>    number of nodes in Cartesian grid\n"
                  "  -P <b:b:b>    periodic boundary condition        (1:1:1)\n"
                  "  -c <r>        cutoff radius\n"
                  "  -n            neutralize net charge              (false)\n"
                  "  -s <n>        random seed                            (1)\n"
                  "  -i <name>     XYZ input file name\n"
                  "  -o <name>     XYZ output file name\n"
                  "  -r <name>     XYZ restart save file name\n"
                  "  -w <n>        step interval of XYZ output            (1)\n"
                  "  -0            XYZ output at step 0               (false)\n"
                  "  -v            print message verbosely                (0)\n"
                  "  -h            show this help message\n"
                  << flush;
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
          }
          case ':': {  // missing option argument
            if (comm_rank == 0)
              cout << "pmt: option requires an argument -- '"
                  << static_cast<char>(::optopt)
                  << "'.  try '-h' for help\n" << flush;
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }
          default: /* case '?': */ {  // unknown option
            if (comm_rank == 0)
              cout << "pmt: unknown option -- '"
                  << static_cast<char>(::optopt)
                  << "'.  try '-h' for help\n" << flush;
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }
        }
      }
      catch (...) {  // invalid option argument
        if (comm_rank == 0)
          cout << "pmt: invalid option argument -- '"
              << argv[::optind - 1]
              << "'.  try '-h' for help\n" << flush;
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }

  int DeterminNumberOfCartNode() {
    using namespace std;
    int num_cart_node = 1;
    for (int i = 0; i < 3; ++i) {
      cart_num[i] = abs(cart_num[i]);
      if (cart_num[i] > 1)
        num_cart_node *= cart_num[i];
    }
    if (num_cart_node > comm_size) {
      if (comm_rank == 0)
        cout << "pmt: number of nodes exceeds communicator size. abort\n"
            << flush;
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (cart_num[0] * cart_num[1] * cart_num[2] == 0) {
      if (cart_num[0] + cart_num[1] + cart_num[2] != 0)
        num_cart_node *= (comm_size / num_cart_node);
      else
        num_cart_node = comm_size;
    }
    return num_cart_node;
  }
};

#endif  // CONF_HH_
