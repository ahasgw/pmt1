#ifndef CONF_HH_
#define CONF_HH_ 1

#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include <unistd.h>
#include <iomanip>
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
  v3r sys_size;
  v3r sys_min;
  v3r sys_max;
  real_t delta_t;
  real_t cutoff;
  v3i cart_num;
  v3i boundary;
  int rest_num;
  int max_step;
  int total_ptcl;
  int origin_center;
  int neutralize;
  int uniformize;
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
      : sys_size(100.0),
        delta_t(0.001),
        cutoff(-1.0),
        cart_num(0),
        boundary(1),
        rest_num(0),
        max_step(1),
        total_ptcl(10000),
        origin_center(0),
        neutralize(0),
        uniformize(0),
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

    sys_min = 0.0;
    sys_max = sys_size;
    if (origin_center) {
      sys_min -= 0.5 * sys_size;
      sys_max -= 0.5 * sys_size;
    }

    DeterminCutoff();

    std::srand(static_cast<unsigned>(global_seed));  // set random seed

    int num_cart_node = DeterminNumberOfCartNode();
    rest_num = comm_size - num_cart_node;
    node_type = ((comm_rank < num_cart_node) ? CART_NODE : REST_NODE);
    MPI_Comm_split(MPI_COMM_WORLD, node_type, 0, &work_comm);

    MPI_Dims_create(num_cart_node, 3, cart_num);
    if (node_type == CART_NODE) {
      v3i periodic;
      for (int d = 0; d < 3; ++d) {
        if (boundary[d] < 0) boundary[d] = 1;       // normalize boundary
        periodic[d] = (boundary[d] == 1 ? 1 : 0);
      }
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
    using namespace std;
    if (c.comm_rank == 0) {
      os << setw(24) << left << "# cmd_line" << c.cmd_line << "\n";
      os << setw(24) << left << "# sys_size" << c.sys_size << "\n";
      os << setw(24) << left << "# sys_min" << c.sys_min << "\n";
      os << setw(24) << left << "# sys_max" << c.sys_max << "\n";
      os << setw(24) << left << "# delta_t" << c.delta_t << "\n";
      os << setw(24) << left << "# cutoff" << c.cutoff << "\n";
      os << setw(24) << left << "# cart_num" << c.cart_num << "\n";
      os << setw(24) << left << "# rest_num" << c.rest_num << "\n";
      os << setw(24) << left << "# boundary" << c.boundary << "\n";
      os << setw(24) << left << "# max_step" << c.max_step << "\n";
      os << setw(24) << left << "# origin_center" << c.origin_center << "\n";
      os << setw(24) << left << "# neutralize" << c.neutralize << "\n";
      os << setw(24) << left << "# uniformize" << c.uniformize << "\n";
      os << setw(24) << left << "# global_seed" << c.global_seed << "\n";
      os << setw(24) << left << "# ifname" << c.ifname << "\n";
      os << setw(24) << left << "# ofname" << c.ofname << "\n";
      os << setw(24) << left << "# rfname" << c.rfname << "\n";
      os << setw(24) << left << "# write_interval" << c.write_interval << "\n";
      os << setw(24) << left << "# write_step0" << c.write_step0 << "\n";
      os << setw(24) << left << "# verbose" << c.verbose << "\n";
      os << setw(24) << left << "# comm_size" << c.comm_size << "\n";
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
    cout.precision(numeric_limits<float>::digits10);

    for (::opterr = 0;;) {
      int opt = ::getopt(argc, argv, ":m:p:S:B:N:d:c:Onus:i:o:r:w:0dvh");
      if (opt == -1) break;
      try {
        switch (opt) {
          case 'm': ReadAbs(max_step); break;
          case 'p': ReadAbs(total_ptcl); break;
          case 'S': Read(sys_size); break;
          case 'B': Read(boundary); break;
          case 'N': Read(cart_num); break;
          case 'd': ReadAbs(delta_t); break;
          case 'c': Read(cutoff); break;
          case 'O': ++origin_center; break;
          case 'n': ++neutralize; break;
          case 'u': ++uniformize; break;
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
                  "  -B <n:n:n>    boundary condition                 (1:1:1)\n"
                  "  -N <X:Y:Z>    number of nodes in Cartesian grid\n"
                  "  -d <r>        delta t                            (0.001)\n"
                  "  -c <r>        cutoff radius\n"
                  "  -O            place center of system at origin          \n"
                  "  -n            neutralize net charge              (false)\n"
                  "  -u            uniformize particle charge & mass  (false)\n"
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
    for (int d = 0; d < 3; ++d) {
      cart_num[d] = abs(cart_num[d]);
      if (cart_num[d] > 1)
        num_cart_node *= cart_num[d];
    }
    if (num_cart_node > comm_size) {
      if (comm_rank == 0)
        cout << "pmt: number of nodes exceeds communicator size. abort\n"
            << flush;
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (cart_num.mul() == 0) {
      if (cart_num.sum() != 0)
        num_cart_node *= (comm_size / num_cart_node);
      else
        num_cart_node = comm_size;
    }
    return num_cart_node;
  }

  void DeterminCutoff() {
    using namespace std;
    if (boundary[0] == 1 || boundary[1] == 1 || boundary[2] == 1) {  // periodic
      v3r max_sys_size(sys_size.max());
      v3r periodic_sys_size = sys_size * v3r(boundary == v3i(1));
      periodic_sys_size += max_sys_size * v3r(boundary != v3i(1));
      real_t half_min_periodic_sys_size = 0.5 * periodic_sys_size.min();
      if (cutoff > half_min_periodic_sys_size) {
        if (comm_rank == 0) {
          cout << "pmt: cutoff radius " << cutoff <<
              " exceeds half of minimum periodic system size " <<
              half_min_periodic_sys_size << ". abort\n" << flush;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
      if (cutoff < 0.0) {
        cutoff = half_min_periodic_sys_size;
      }
    } else {  // non-periodic
      if (cutoff < 0.0) {
        if (comm_rank == 0)
          cout << "pmt: cutoff radius is not given. abort\n" << flush;
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
  }
};

#endif  // CONF_HH_
