#include "cart.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <ostream>
#include "conf.hh"
#include "output.hh"
#include "ptcl.hh"
#include "random.hh"

CartNode::CartNode(const Conf &conf): conf_(conf), os(NULL) {
  using namespace std;
  // setup timer
  t_comm.Label("cart comm").Comm(conf_.cart_comm);
  t_step.Label("cart step").Comm(conf_.cart_comm);
  t_init.Label("cart init").Comm(conf_.cart_comm).Start();

  // make local copy
  cart_comm = conf_.cart_comm;
  sys_size = conf_.sys_size;
  sys_min = conf_.sys_min;
  sys_max = conf_.sys_max;

  MPI_Comm_rank(cart_comm, &cart_rank);
  MPI_Comm_size(cart_comm, &cart_size);
  MPI_Cart_coords(cart_comm, cart_rank, 3, cart_pos);

  if (cart_rank == 0) {
    if (conf_.verbose > 0) std::cout << "# cart_size\t" << cart_size << "\n";
    if (!conf_.ofname.empty()) {
      os = new ofstream(conf_.ofname.c_str());  // prepare output stream
      if (!*os) {
        cout << "pmt0: cannot open file '" << conf_.ofname << "'. abort\n"
            << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }
  }

  // setup node
  div_min = div_max = sys_min;
  div_min += v3r(cart_pos    ) * sys_size / v3r(conf_.cart_num);
  div_max += v3r(cart_pos + 1) * sys_size / v3r(conf_.cart_num);

  // topology
  InitConnect();

  // generate particles
  GenerateParticles();

  // write particle coordinates at step 0
  if (!conf_.ofname.empty()) {
    OutputXYZ(*os, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
              0, conf_.max_step, cart_rank, cart_size, cart_comm);
    steps_to_write = conf_.write_interval;
  }

  t_init.Stop();
}

CartNode::~CartNode() {
  if (os) delete os;

  // print timer
  if (conf_.verbose > 1) t_comm.PrintAll("# ", conf_.max_step);
  if (conf_.verbose > 0) t_comm.PrintMax("# max ", conf_.max_step);
  if (conf_.verbose > 1) t_step.PrintAll("# ", conf_.max_step);
  if (conf_.verbose > 0) t_step.PrintMax("# max ", conf_.max_step);
  if (conf_.verbose > 1) t_init.PrintAll("# ");
  if (conf_.verbose > 0) t_init.PrintMax("# max ");
}

void CartNode::StepForward(int t) {
  using namespace std;
  t_step.Start();
  Ptcls::size_type p_size = ptcls.size();
//#pragma omp parallel for
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    Ptcl &ptcl = ptcls[p];
    ptcl.attr = 0;
    ptcl.crd += ptcl.vel;

    // embarkation check & periodic shift
    for (int i = 0; i < 3; ++i) {
      ptcl.attr <<= 2;
      if (ptcl.crd[i] < div_min[i]) {
        if (ptcl.crd[i] < sys_min[i]) { ptcl.crd[i] += sys_size[i]; }
        ptcl.attr |= LOWER_DIR;
      }
      else if (ptcl.crd[i] >= div_max[i]) {
        if (ptcl.crd[i] >= sys_max[i]) { ptcl.crd[i] -= sys_size[i]; }
        ptcl.attr |= UPPER_DIR;
      }
    }
  }
  t_comm.Start();
  ExchangeParticles();
  t_comm.Stop();
  t_step.Stop();

  // write particle coordinates at step t
  if (!conf_.ofname.empty()) {
    if (steps_to_write <= 1) {
      OutputXYZ(*os, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
                t, conf_.max_step, cart_rank, cart_size, cart_comm);
      steps_to_write = conf_.write_interval;
    } else {
      --steps_to_write;
    }
  }
}

void CartNode::InitConnect() {
  const unsigned tag[3] = { MIDDLE_DIR, LOWER_DIR, UPPER_DIR };
  const int offset[3] = { 0, -1, 1 };
  v3i coords;
  conns.reserve(26);
  for (int x = 0; x < 3; ++x) {
    for (int y = 0; y < 3; ++y) {
      for (int z = 0; z < 3; ++z) {
        if (x == 0 && y == 0 && z == 0) continue;
        Connect conn;
        conn.dir_tag = ((tag[x] << 4) | (tag[y] << 2) | tag[z]);
        //
        coords = cart_pos;
        coords[0] += offset[x];
        coords[1] += offset[y];
        coords[2] += offset[z];
        MPI_Cart_rank(cart_comm, coords, &conn.send_to);
        //
        coords = cart_pos;
        coords[0] -= offset[x];
        coords[1] -= offset[y];
        coords[2] -= offset[z];
        MPI_Cart_rank(cart_comm, coords, &conn.recv_from);
        //
        conns.push_back(conn);
      }
    }
  }
}

void CartNode::GenerateParticles() {
  using namespace std;
  double min_sys_size = min(min(sys_size[0], sys_size[1]), sys_size[2]);
  for (int n = 0; n < conf_.total_ptcl; ++n) {
    Ptcl p;
    for (int i = 0; i < 3; ++i) {
      p.crd[i] = Rand() * sys_size[i] + sys_min[i];
      p.vel[i] = Gaussian(0.0, 0.5) * min_sys_size / 1024.0;
    }
    if (div_min <= p.crd && p.crd < div_max) {
      p.id = n;
      ptcls.push_back(p);
    }
  }
}

void CartNode::ExchangeParticles() {
  using namespace std;
  sort(ptcls.begin(), ptcls.end(), Less<Ptcl::DIR>());

  // find lower bound of outgoing ptcl
  Ptcl stayin;
  stayin.attr = 0x00;
  Ptcls::iterator lb_outgoing =
      upper_bound(ptcls.begin(), ptcls.end(), stayin, Less<Ptcl::DIR>());

  // copy outgoing ptcls
  Ptcls send_buff(lb_outgoing, ptcls.end());
  // trim ptcls
  ptcls.erase(lb_outgoing, ptcls.end());

  for (int i = 0; i < 26; ++i) {
    // get range of ptcl whose attr equals a dir_tag
    Ptcl p;
    p.attr = conns[i].dir_tag;
    pair<Ptcls::iterator, Ptcls::iterator> range =
        equal_range(send_buff.begin(), send_buff.end(), p, Less<Ptcl::DIR>());

    // get offsets from ptcls[0] and number of Ptcl to send
    int first = static_cast<int>(distance(send_buff.begin(), range.first));
    int count = static_cast<int>(distance(range.first, range.second));
    int bytes = count * sizeof(Ptcl);

    // send immediately
    MPI_Isend(&send_buff[first], bytes, MPI_BYTE, conns[i].send_to,
              static_cast<int>(p.attr), cart_comm, &conns[i].req);
  }

  for (int i = 0; i < 26; ++i) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, cart_comm, &status);

    // get size of message
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    int count = bytes / sizeof(Ptcl);

    // receive using i-th buffer
    recv_buff[i].resize(count);
    MPI_Recv(&recv_buff[i][0], bytes, MPI_BYTE, status.MPI_SOURCE,
             status.MPI_TAG, cart_comm, MPI_STATUS_IGNORE);

    // append incoming ptcls to ptcls
    ptcls.insert(ptcls.end(), recv_buff[i].begin(), recv_buff[i].end());
  }

  // wait all send
  for (int i = 0; i < 26; ++i) {
    MPI_Wait(&conns[i].req, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(cart_comm);
}
