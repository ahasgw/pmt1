#include "cart.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <ostream>
#include "conf.hh"
#include "input.hh"
#include "output.hh"
#include "ptcl.hh"
#include "random.hh"

CartNode::CartNode(Conf &conf): conf_(conf) {
  using namespace std;
  // setup timer
  t_cfrc.Label("cart force").Comm(conf_.cart_comm);
  t_comm.Label("cart comm").Comm(conf_.cart_comm);
  t_step.Label("cart step").Comm(conf_.cart_comm);
  t_oput.Label("cart output").Comm(conf_.cart_comm);
  t_init.Label("cart init").Comm(conf_.cart_comm).Start();

  // make local copy
  cart_comm = conf_.cart_comm;
  sys_size = conf_.sys_size;
  sys_size_2 = 0.5 * sys_size;
  sys_min = conf_.sys_min;
  sys_max = conf_.sys_max;
  dt = conf_.delta_t;
  cutoff2 = conf_.cutoff * conf_.cutoff;

  MPI_Comm_rank(cart_comm, &cart_rank);
  MPI_Comm_size(cart_comm, &cart_size);
  MPI_Cart_coords(cart_comm, cart_rank, 3, cart_pos);

  boundary = conf_.boundary;

  if (cart_rank == 0) {
    if (conf_.verbose > 0)
      cout << setw(24) << left << "# cart_size" << cart_size << "\n";

    // prepare input stream
    if (!conf_.ifname.empty()) {
      is.open(conf_.ifname.c_str());
      if (!is) {
        cout << "pmt: cannot open input file '" << conf_.ifname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }

    // prepare output stream
    if (!conf_.ofname.empty()) {
      os.open(conf_.ofname.c_str());
      if (!os) {
        cout << "pmt: cannot open output file '" << conf_.ofname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }

    // prepare restart save stream
    if (!conf_.rfname.empty()) {
      rs.open(conf_.rfname.c_str());
      if (!rs) {
        cout << "pmt: cannot open restart save file '" << conf_.rfname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }
  }

  // setup node
  div_size = sys_size / v3r(conf_.cart_num);
  div_size_2 = 0.5 * div_size;
  div_min = div_max = sys_min;
  div_min += v3r(cart_pos    ) * div_size;
  div_max += v3r(cart_pos + 1) * div_size;
  div_intr = v3r(conf_.cutoff) / div_size;
  div_intr += v3r(conf_.cutoff) > (v3r(div_intr) * div_size);
  v3i intr_region = (div_intr * 2) + 1;
  intr_size = intr_region.mul() - 1;
  use_ringcomm = (cart_size - 1 < intr_size);  // <= ?
  //
  if (cart_rank == 0) {
    if (conf_.verbose > 0)
      cout << setw(24) << left << "# intr_size" << intr_size << "\n";
  }

  // topology
  InitConnect();
  if (!use_ringcomm) InitInteract();

  // generate particles
  GenerateParticles();

  // write particle coordinates at step 0
  if (!conf_.ofname.empty() || !conf_.rfname.empty()) {
    t_oput.Start();
    if (conf_.write_step0 || ((conf_.max_step == 0) && rs)) {
      OutputXYZ(os, rs, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
                0, conf_.max_step, cart_rank, cart_size, cart_comm);
    }
    steps_to_write = conf_.write_interval;
    t_oput.Stop();
  }

  // calculate initial force
  if (conf_.max_step > 0) {
    t_cfrc.Start();
    CalculateForce();
    t_cfrc.Stop();
  }

  t_init.Stop();
}

CartNode::~CartNode() {
  // print timer
  if (conf_.verbose > 1) t_cfrc.PrintAll("# ");
  if (conf_.verbose > 0) t_cfrc.PrintMax("# max ");
  if (conf_.verbose > 1) t_comm.PrintAll("# ");
  if (conf_.verbose > 0) t_comm.PrintMax("# max ");
  if (conf_.verbose > 1) t_step.PrintAll("# ");
  if (conf_.verbose > 0) t_step.PrintMax("# max ");
  if (conf_.verbose > 1) t_oput.PrintAll("# ");
  if (conf_.verbose > 0) t_oput.PrintMax("# max ");
  if (conf_.verbose > 1) t_init.PrintAll("# ");
  if (conf_.verbose > 0) t_init.PrintMax("# max ");
}

void CartNode::StepForward(int t) {
  using namespace std;
  t_step.Start();
  Ptcls::size_type p_size = ptcls.size();

  // update positions and velocities
  // v(t+0.5dt) = v(t) + 0.5dt*Fi(t)/mi
  // r(t+dt) = r(t) + dt*v(t+0.5dt)
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    Ptcl &ptcl = ptcls[p];
    ptcl.attr = 0;
    ptcl.crd += dt * (ptcl.vel += (dt * force[p] * ptcl.inv_2mass));

    // embarkation check
    for (int d = 0; d < 3; ++d) {
      ptcl.attr <<= 2;
      if (ptcl.crd[d] < div_min[d]) {
        switch (boundary[d]) {
          case 2: {  // reflecting wall
            if (ptcl.crd[d] < sys_min[d]) {
              ptcl.crd[d] = 2.0 * sys_min[d] - ptcl.crd[d];
              ptcl.vel[d] *= -1.0;
            }
            else ptcl.attr |= LOWER_DIR;
            break;
          }
          case 1: {  // periodic
            if (ptcl.crd[d] < sys_min[d]) ptcl.crd[d] += sys_size[d];
            ptcl.attr |= LOWER_DIR;
            break;
          }
          case 0:
          default: {
            if (ptcl.crd[d] >= sys_min[d]) ptcl.attr |= LOWER_DIR;
          }
        }
      }
      else if (ptcl.crd[d] >= div_max[d]) {
        switch (boundary[d]) {
          case 2: {  // reflecting wall
            if (ptcl.crd[d] >= sys_max[d]) {
              ptcl.crd[d] = 2.0 * sys_max[d] - ptcl.crd[d];
              ptcl.vel[d] *= -1.0;
            }
            else ptcl.attr |= UPPER_DIR;
            break;
          }
          case 1: {  // periodic
            if (ptcl.crd[d] >= sys_max[d]) ptcl.crd[d] -= sys_size[d];
            ptcl.attr |= UPPER_DIR;
            break;
          }
          case 0:
          default: {
            if (ptcl.crd[d] < sys_max[d]) ptcl.attr |= UPPER_DIR;
          }
        }
      }
    }
  }

  // calculate force F(t+0.5dt) using r(t+dt)
  // F(t+0.5dt)
  t_cfrc.Start();
  CalculateForce();
  t_cfrc.Stop();

  // update velocities
  // v(t+dt) = v(t+0.5dt) + 0.5dt*Fi(t+dt)/mi
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    Ptcl &ptcl = ptcls[p];
    ptcl.vel += (dt * force[p] * ptcl.inv_2mass);
  }

  t_comm.Start();
  ExchangeParticles();
  t_comm.Stop();
  t_step.Stop();

  // write particle coordinates at step t
  if (!conf_.ofname.empty() || !conf_.rfname.empty()) {
    t_oput.Start();
    if ((steps_to_write <= 1) || ((conf_.max_step == t) && rs)) {
      OutputXYZ(os, rs, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
                t, conf_.max_step, cart_rank, cart_size, cart_comm);
      steps_to_write = conf_.write_interval;
    } else {
      --steps_to_write;
    }
    t_oput.Stop();
  }
}

void CartNode::InitConnect() {
  const unsigned tag[3] = { MIDDLE_DIR, LOWER_DIR, UPPER_DIR };
  const int ofst[3] = { 0, -1, 1 };
  v3i coords;
  conns.reserve(26);
  for (int x = 0; x < 3; ++x) {
    for (int y = 0; y < 3; ++y) {
      for (int z = 0; z < 3; ++z) {
        if (x == 0 && y == 0 && z == 0) continue;
        Connect conn;
        conn.dir_tag   = ((tag[x] << 4) | (tag[y] << 2) | tag[z]);
        conn.send_to   = GetCartRankAtOffset( ofst[x],  ofst[y],  ofst[z]);
        conn.recv_from = GetCartRankAtOffset(-ofst[x], -ofst[y], -ofst[z]);
        conns.push_back(conn);
      }
    }
  }
}

int CartNode::GetCartRankAtOffset(int dx, int dy, int dz) {
  v3i coords = cart_pos;
  coords[0] += dx;
  coords[1] += dy;
  coords[2] += dz;
  // adjust coords
  for (int d = 0; d < 3; ++d) {
    if (boundary[d] != 1) {  // non-periodic
      if (coords[d] < 0) {
        coords[d] = 0;
      } else if (coords[d] >= conf_.cart_num[d]) {
        coords[d] = conf_.cart_num[d] - 1;
      }
    }
  }
  // get cart rank at the coords
  int rank = MPI_PROC_NULL;
  MPI_Cart_rank(cart_comm, coords, &rank);
  return rank;
}

void CartNode::InitInteract() {
  using namespace std;
  intrs.reserve(intr_size);
  for (int r = 1; r <= div_intr.max(); ++r) {
    int mx = min(r, div_intr[0]);
    int my = min(r, div_intr[1]);
    int mz = min(r, div_intr[2]);
    for (int x = -mx; x <= mx; ++x) {
      for (int y = -my; y <= my; ++y) {
        for (int z = -mz; z <= mz; ++z) {
          if ((-mx < x) && (x < mx) &&
              (-my < y) && (y < my) &&
              (-mz < z) && (z < mz))
            continue;
          Interact intr;
          intr.send_to   = GetCartRankAtOffsetInsideBoundary( x,  y,  z);
          intr.recv_from = GetCartRankAtOffsetInsideBoundary(-x, -y, -z);
          intrs.push_back(intr);
        }
      }
    }
  }
}

int CartNode::GetCartRankAtOffsetInsideBoundary(int dx, int dy, int dz) {
  v3i coords = cart_pos;
  coords[0] += dx;
  coords[1] += dy;
  coords[2] += dz;
  // check if coords is out of bounds
  bool outofbounds = false;
  for (int d = 0; d < 3; ++d) {
    if (boundary[d] != 1) {  // non-periodic
      if (coords[d] < 0 || coords[d] >= conf_.cart_num[d]) {
        outofbounds = true;
      }
    }
  }
  // get cart rank at the coords
  int rank = MPI_PROC_NULL;
  if (!outofbounds) {
    MPI_Cart_rank(cart_comm, coords, &rank);
  }
  return rank;
}

void CartNode::GenerateParticles() {
  using namespace std;
  if (!conf_.ifname.empty()) {
    // read from XYZ input file
    InputXYZ(is, &ptcls, &conf_.total_ptcl, div_min, div_max,
             sys_min, sys_max, boundary,
             cart_rank, cart_size, cart_comm);
  } else {
    // generate at random
    real_t chg = 0.0;
    const int last_id = conf_.total_ptcl - 1;
    for (int n = 0; n < conf_.total_ptcl; ++n) {
      Ptcl p;
      for (int d = 0; d < 3; ++d) {
        p.crd[d] = Rand() * sys_size[d] + sys_min[d];
        p.vel[d] = Gaussian(0.0, 0.5);
      }

      if (conf_.uniformize) {
        if (conf_.neutralize) {
          chg = ((n & 1) == 0)
              ? ((n < last_id) ? 20.0 : 0.0)
              : -chg;
        } else {
          chg = 20.0;
        }
        p.chg = chg;
        p.inv_2mass = 0.05;
      } else {
        if (conf_.neutralize) {
          chg = ((n & 1) == 0)
              ? ((n < last_id) ? (Gaussian(0.0, 0.2) * 100) : 0.0)
              : -chg;
        } else {
          chg = Gaussian(0.0, 0.2) * 100;
        }
        p.chg = chg;
        p.inv_2mass = 0.5 / (Rand() * 15.0 + 1.0);
      }

      if ((div_min <= p.crd).mul() * (p.crd < div_max).mul()) {
        p.id = n;
        ptcls.push_back(p);
      }
    }
  }
  force.resize(ptcls.size());

  if (cart_rank == 0) {
    if (conf_.verbose > 0)
      cout << setw(24) << left << "# total_ptcl" << conf_.total_ptcl << "\n";
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

  force.resize(ptcls.size());
}

void CartNode::CalculateForce() {
#if 1
  if (!use_ringcomm) {
    CalculateForceCutoffPeriodic();
  } else {
    CalculateForceCutoffPeriodic_RingComm();
  }
#else
  force.assign(ptcls.size(), 0.0);
#endif
}

// n-neighber communication (use Intrs)
void CartNode::CalculateForceCutoffPeriodic() {
  static const double kInv4Pi = 0.25 / acos(-1.0);
  Ptcls::size_type p_size = ptcls.size();

  // clear force
  force.assign(p_size, 0.0);

  // prepare send buffer
  std::vector<v4r> send_buffer(conf_.total_ptcl);

  // pack my particles in send_buffer
#pragma omp parallel for schedule(static)
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    send_buffer[p][0] = ptcls[p].crd[0];
    send_buffer[p][1] = ptcls[p].crd[1];
    send_buffer[p][2] = ptcls[p].crd[2];
    send_buffer[p][3] = ptcls[p].chg;
  }

  Ptcls::size_type j_size = p_size;

  // calculate force from particles of local node
#pragma omp parallel for schedule(dynamic, 8)
  for (Ptcls::size_type i = 0; i < p_size; ++i) {
    for (Ptcls::size_type j = 0; j < i; ++j) {
      v3r r;
      for (int d = 0; d < 3; ++d) {
        r[d] = ptcls[i].crd[d] - send_buffer[j][d];
        if (boundary[d] == 1) {
          if (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
          if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
        }
      }
      double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      if (r2 < cutoff2) {
        double r6 = r2 * r2 * r2;
        double _r3 = 1.0 / sqrt(r6);  // divzero ?
        double chgs = ptcls[i].chg * send_buffer[j][3];
        r *= (kInv4Pi * chgs) * _r3;
        force[i] += r;
        force[j] -= r;
      }
    }
  }

  // prepare recv buffer
  std::vector<v4r> recv_buffer(conf_.total_ptcl);

  // calculate force from particles of external node
  for (Intrs::size_type i = 0; i < intrs.size(); ++i) {
    // send
    MPI_Request req;
    MPI_Isend(&send_buffer[0], j_size * sizeof(v4r), MPI_BYTE,
              intrs[i].send_to, cart_rank, cart_comm, &req);

    // recv
    MPI_Status status;
    MPI_Recv(&recv_buffer[0], conf_.total_ptcl * sizeof(v4r), MPI_BYTE,
             intrs[i].recv_from, intrs[i].recv_from, cart_comm, &status);

    // wait
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    // get size of receive message
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    j_size = bytes / sizeof(v4r);

    // calculate force
    if (intrs[i].recv_from != MPI_PROC_NULL) {
#pragma omp parallel for collapse(2), schedule(dynamic, 4)
      for (Ptcls::size_type i = 0; i < p_size; ++i) {
        for (Ptcls::size_type j = 0; j < j_size; ++j) {
          v3r r;
          for (int d = 0; d < 3; ++d) {
            r[d] = ptcls[i].crd[d] - recv_buffer[j][d];
            if (boundary[d] == 1) {
              if (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
              if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
            }
          }
          double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
          if (r2 < cutoff2) {
            double r6 = r2 * r2 * r2;
            double _r3 = 1.0 / sqrt(r6);  // divzero ?
            double chgs = ptcls[i].chg * recv_buffer[j][3];
            r *= (kInv4Pi * chgs) * _r3;
            force[i] += r;
          }
        }
      }
    }
  }
}

// ring communication (do not use Intrs)
void CartNode::CalculateForceCutoffPeriodic_RingComm() {
  static const double kInv4Pi = 0.25 / acos(-1.0);
  Ptcls::size_type p_size = ptcls.size();

  // clear force
  force.assign(p_size, 0.0);

  // prepare recieve buffer
  std::vector<v4r> j_ptcls(conf_.total_ptcl);

  // pack my particles in j_ptcls
#pragma omp parallel for schedule(static)
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    j_ptcls[p][0] = ptcls[p].crd[0];
    j_ptcls[p][1] = ptcls[p].crd[1];
    j_ptcls[p][2] = ptcls[p].crd[2];
    j_ptcls[p][3] = ptcls[p].chg;
  }

  Ptcls::size_type j_size = p_size;

  // calculate force from particles of local node
#pragma omp parallel for schedule(dynamic, 8)
  for (Ptcls::size_type i = 0; i < p_size; ++i) {
    for (Ptcls::size_type j = 0; j < i; ++j) {
      v3r r;
      for (int d = 0; d < 3; ++d) {
        r[d] = ptcls[i].crd[d] - j_ptcls[j][d];
        if (boundary[d] == 1) {
          if (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
          if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
        }
      }
      double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      if (r2 < cutoff2) {
        double r6 = r2 * r2 * r2;
        double _r3 = 1.0 / sqrt(r6);  // divzero ?
        double chgs = ptcls[i].chg * j_ptcls[j][3];
        r *= (kInv4Pi * chgs) * _r3;
        force[i] += r;
        force[j] -= r;
      }
    }
  }

  // prepare send buffer
  std::vector<v4r> send_buffer;

  // calculate force from particles of external node
  for (int i = 1; i < cart_size; ++i) {
    // copy j_ptcls to send_buffer
    send_buffer.resize(j_size);
    for (Ptcls::size_type j = 0; j < j_size; ++j) {
      send_buffer[j] = j_ptcls[j];
    }

    // send
    int dest = (cart_rank + 1) % cart_size;
    MPI_Request req;
    MPI_Isend(&send_buffer[0], j_size * sizeof(v4r), MPI_BYTE,
              dest, cart_rank, cart_comm, &req);

    // recv
    int source = (cart_rank + cart_size - 1) % cart_size;
    MPI_Status status;
    MPI_Recv(&j_ptcls[0], conf_.total_ptcl * sizeof(v4r), MPI_BYTE,
             source, source, cart_comm, &status);

    // wait
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    // get size of receive message
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    j_size = bytes / sizeof(v4r);

    // calculate force
#pragma omp parallel for collapse(2), schedule(dynamic, 4)
    for (Ptcls::size_type i = 0; i < p_size; ++i) {
      for (Ptcls::size_type j = 0; j < j_size; ++j) {
        v3r r;
        for (int d = 0; d < 3; ++d) {
          r[d] = ptcls[i].crd[d] - j_ptcls[j][d];
          if (boundary[d] == 1) {
            if (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
            if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
          }
        }
        double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        if (r2 < cutoff2) {
          double r6 = r2 * r2 * r2;
          double _r3 = 1.0 / sqrt(r6);  // divzero ?
          double chgs = ptcls[i].chg * j_ptcls[j][3];
          r *= (kInv4Pi * chgs) * _r3;
          force[i] += r;
        }
      }
    }
  }
}
